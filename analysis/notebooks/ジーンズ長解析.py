# -*- coding: utf-8 -*-
# %% cell1: ジーンズ長の診断
"""等温・Cartesian・static SMRのAthena++ VTKを元のセル幅で診断する。

実行例:
python3 analysis/notebooks/ジーンズ長解析.py /work/beta/aida/results/Toyouchi-test72

依存: numpy, matplotlib, pyvista。スクリプト冒頭の設定でも実行可能。
シンク判定・任意の解析箱への所属はセル中心で判定する。
VTKはghost_zones=false、重複のない全leaf MeshBlock一式を使う。
音速・重力定数・シンク半径・単位系は冒頭で手動設定する。
領域はVTKから推定するため、全時刻で同じ外側領域が欠けている場合は検出できない。
"""

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

# IDEから実行するときの設定。対象runの値に変更する。
DATADIR = "~/athena-project/results/〇〇"
ISO_SOUND_SPEED_CODE = 0.707
GRAV_CONST_CODE = 1.0
SINK_RADIUS_AU = 1000.0  # 対象runに合わせる。シンクなしの場合は0。
OUTPUT_DIR = "./jeans_diagnostics"
THRESHOLDS = (4, 8, 16, 32)
# 全領域を診断。中心のみなら例: ((-5000, 5000),) * 3 [au]
ANALYSIS_BOUNDS_AU = None
SLICE_HALF_WIDTH_AU = 5000.0
L_UNIT_CGS = 7.03e15
M_UNIT_CGS = 4.0e33
T_UNIT_CGS = 3.61e10
AU_CGS = 1.495978707e13
L_AU = L_UNIT_CGS / AU_CGS
T_KYR = T_UNIT_CGS / (365.25 * 86400 * 1000)


def vtk_time(path):
    with path.open("rb") as handle:
        handle.readline()
        title = handle.readline().decode("ascii", errors="replace")
    match = re.search(r"time\s*=\s*([+\-\d.eE]+)", title)
    if not match:
        raise ValueError(f"VTK timeがありません: {path}")
    return float(match.group(1))


def group_files(directory, output_number):
    pattern = re.compile(rf"\.block(\d+)\.out{output_number}\.(\d+)\.vtk$")
    groups = defaultdict(dict)
    for path in sorted(directory.rglob(f"*.block*.out{output_number}.*.vtk")):
        match = pattern.search(path.name)
        if match:
            block, step = map(int, match.groups())
            if block in groups[step]:
                raise ValueError(f"重複block {block}, step {step}: {path}")
            groups[step][block] = path
    if not groups:
        raise FileNotFoundError(f"MeshBlock VTKがありません: {directory}")
    expected = set(groups[min(groups)])
    for step, files in groups.items():
        if set(files) != expected:
            raise ValueError(f"static SMRのblock集合が一致しません: step={step}")
    return groups


def block_arrays(grid):
    """VTKセル順序（x最速）を保持。非一様な直交格子にも対応。"""
    if not isinstance(grid, pv.RectilinearGrid):
        raise ValueError("Athena++ Cartesian RECTILINEAR_GRID VTKが必要です")
    edges = [np.asarray(getattr(grid, axis)) for axis in "xyz"]
    widths = [np.diff(edge) for edge in edges]
    if any(np.any(w <= 0) or not np.all(np.isfinite(w)) for w in widths):
        raise ValueError("不正なセル幅です")
    centers = [(edge[:-1] + edge[1:]) / 2 for edge in edges]
    xyz = np.stack([a.ravel(order="F") for a in
                    np.meshgrid(*centers, indexing="ij")], axis=1)
    delta = np.stack([a.ravel(order="F") for a in
                      np.meshgrid(*widths, indexing="ij")], axis=1)
    names = {name.lower(): name for name in grid.cell_data}
    key = next((names[name] for name in ("rho", "dens", "density", "prim_dens")
                if name in names), None)
    if key is None:
        raise KeyError(f"密度配列がありません: {list(grid.cell_data)}")
    rho = np.asarray(grid.cell_data[key]).reshape(-1)
    if len(rho) != len(xyz) or not np.all(np.isfinite(rho)) or np.any(rho <= 0):
        raise ValueError("密度のサイズ不一致、非有限値または非正値があります")
    return xyz, delta, rho


def draw_slice(ax, xyz, delta, nj, selected, args, norm=None, cmap="viridis"):
    normal = "xyz".index(args.normal)
    u, v = [axis for axis in range(3) if axis != normal]
    lower, upper = xyz - delta / 2, xyz + delta / 2
    plane = args.slice_au / L_AU
    # 面上では正側のセルだけを選び、ブロック境界での二重描画を避ける。
    mask = selected & (lower[:, normal] <= plane) & (plane < upper[:, normal])
    half = args.half_width_au / L_AU
    mask &= ((lower[:, u] < half) & (upper[:, u] > -half)
             & (lower[:, v] < half) & (upper[:, v] > -half))
    ids = np.flatnonzero(mask)
    if not len(ids):
        return 0
    polygons = [np.array([[lower[i, u], lower[i, v]],
                          [upper[i, u], lower[i, v]],
                          [upper[i, u], upper[i, v]],
                          [lower[i, u], upper[i, v]]]) * L_AU for i in ids]
    # 元のセルをそのまま描画する。補間は行わない。
    collection = PolyCollection(polygons, array=nj[ids], cmap=cmap,
                                norm=LogNorm(1, 64) if norm is None else norm,
                                edgecolors="none")
    ax.add_collection(collection)
    bounds_lo, bounds_hi = lower.min(axis=0), upper.max(axis=0)
    ax.add_patch(Rectangle((bounds_lo[u]*L_AU, bounds_lo[v]*L_AU),
                          (bounds_hi[u]-bounds_lo[u])*L_AU,
                          (bounds_hi[v]-bounds_lo[v])*L_AU,
                          fill=False, edgecolor="black", linewidth=0.45, alpha=0.6))
    return len(ids)


def analyze(args):
    cs = float(ISO_SOUND_SPEED_CODE)
    grav = float(GRAV_CONST_CODE)
    sink = float(SINK_RADIUS_AU)
    if not (np.isfinite(cs) and np.isfinite(grav) and cs > 0 and grav > 0
            and np.isfinite(sink) and sink >= 0):
        raise ValueError("音速・重力定数・シンク半径を確認してください")
    region = None if ANALYSIS_BOUNDS_AU is None else np.asarray(ANALYSIS_BOUNDS_AU) / L_AU
    if region is not None and (region.shape != (3, 2) or
                              not np.all(np.isfinite(region)) or
                              np.any(region[:, 1] <= region[:, 0])):
        raise ValueError("ANALYSIS_BOUNDS_AUを確認してください")
    groups = group_files(Path(args.datadir).expanduser(), args.output_number)
    # 最初の出力の全MeshBlockから直方体領域を推定する。
    domain = np.array([[np.inf, -np.inf]] * 3)
    for path in groups[min(groups)].values():
        grid = pv.read(path)
        bounds = np.asarray(grid.bounds).reshape(3, 2)
        domain[:, 0] = np.minimum(domain[:, 0], bounds[:, 0])
        domain[:, 1] = np.maximum(domain[:, 1], bounds[:, 1])
        del grid
    if not np.all(np.isfinite(domain)) or np.any(domain[:, 1] <= domain[:, 0]):
        raise ValueError("VTKから有効な3次元領域を取得できません")
    expected_volume = np.prod(domain[:, 1] - domain[:, 0])
    out = Path(args.output).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    (out / "settings.json").write_text(json.dumps(dict(
        arguments=vars(args), sound_speed_code=cs, G_code=grav, sink_radius_au=sink,
        thresholds=THRESHOLDS, analysis_bounds_au=ANALYSIS_BOUNDS_AU,
        inferred_domain_code=domain.tolist(),
        L_unit_cgs=L_UNIT_CGS, M_unit_cgs=M_UNIT_CGS, T_unit_cgs=T_UNIT_CGS,
        selection="cell centers; exclude r < sink radius; all leaf blocks; no ghosts"
    ), indent=2, ensure_ascii=False))
    rows = []
    for index, (step, files) in enumerate(sorted(groups.items())):
        time = vtk_time(next(iter(files.values())))
        do_slice = index % args.slice_every == 0
        if do_slice:
            fig, ax = plt.subplots(figsize=(8, 7))
        mass = volume = 0.0
        bad_mass = np.zeros(len(THRESHOLDS))
        minimum = np.inf
        worst = None
        count = slice_count = 0
        for path in files.values():
            if not np.isclose(vtk_time(path), time, rtol=1e-10, atol=1e-12):
                raise ValueError(f"同じstepの時刻が一致しません: {path}")
            grid = pv.read(path)
            xyz, delta, rho = block_arrays(grid)
            dv = delta.prod(axis=1)
            volume += dv.sum()
            if np.any(xyz-delta/2 < domain[:, 0]-1e-5) or np.any(xyz+delta/2 > domain[:, 1]+1e-5):
                raise ValueError(f"最初のVTKで推定した領域外のセルがあります: {path}")
            selected = np.linalg.norm(xyz, axis=1) >= sink / L_AU
            if region is not None:
                selected &= np.all((xyz >= region[:, 0]) & (xyz < region[:, 1]), axis=1)
            nj = np.sqrt(np.pi * cs**2 / (grav * rho)) / delta.max(axis=1)
            cell_mass = rho * dv
            mass += cell_mass[selected].sum()
            count += int(selected.sum())
            for k, threshold in enumerate(THRESHOLDS):
                bad_mass[k] += cell_mass[selected & (nj < threshold)].sum()
            if np.any(selected):
                i = np.argmin(np.where(selected, nj, np.inf))
                if nj[i] < minimum:
                    minimum = float(nj[i])
                    worst = (*xyz[i]*L_AU, rho[i], delta[i].max()*L_AU, path.name)
            if do_slice:
                slice_count += draw_slice(ax, xyz, delta, nj, selected, args)
            del grid
        # 推定した直方体領域と全セル体積を比較し、欠損・重複をチェック。
        if not np.isclose(volume, expected_volume, rtol=1e-5):
            raise ValueError(f"全セル体積がVTKの推定領域と不一致（欠損・重複・ghostを確認）: step={step}, {volume}/{expected_volume}")
        if mass <= 0 or worst is None:
            raise ValueError(f"解析領域にシンク外のガスがありません: step={step}")
        row = dict(step=step, time_code=time, time_kyr=time*T_KYR, min_NJ=minimum,
                   mass_code=mass, selected_cells=count, meshblocks=len(files),
                   min_x_au=worst[0], min_y_au=worst[1], min_z_au=worst[2],
                   min_rho_code=worst[3], min_dx_au=worst[4], min_block=worst[5])
        row.update({f"mass_fraction_NJ_lt_{n}": bad_mass[k]/mass
                    for k, n in enumerate(THRESHOLDS)})
        rows.append(row)
        if do_slice:
            axes = [a for a in "xyz" if a != args.normal]
            ax.set(xlim=(-args.half_width_au, args.half_width_au),
                   ylim=(-args.half_width_au, args.half_width_au),
                   xlabel=f"{axes[0]} [au]", ylabel=f"{axes[1]} [au]",
                   title=f"t={time*T_KYR:.3f} kyr; {args.normal}={args.slice_au:g} au\n"
                         "Native cells; black: MeshBlock edges (including SMR boundaries)")
            ax.set_aspect("equal")
            ax.set_facecolor("0.85")
            radius2 = sink**2 - args.slice_au**2
            if radius2 > 0:
                ax.add_patch(Circle((0, 0), np.sqrt(radius2), fill=False,
                                    color="red", linestyle="--", label="Sink radius"))
                ax.legend(loc="upper right")
            if not slice_count:
                ax.text(.5, .5, "No selected cells on this plane", transform=ax.transAxes, ha="center")
            colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=LogNorm(1, 64), cmap="viridis"),
                                   ax=ax, extend="both", ticks=[1, 2, 4, 8, 16, 32, 64])
            colorbar.set_label(r"$N_J=\lambda_J/\max(\Delta x_i)$")
            fig.savefig(out / f"jeans_slice_{step:05d}.png", dpi=180, bbox_inches="tight")
            plt.close(fig)
        print(f"[{index+1}/{len(groups)}] step={step} t={time*T_KYR:.3f} kyr "
              f"min NJ={minimum:.5g}, mass(NJ<4)/mass={bad_mass[0]/mass:.5g}")
    rows.sort(key=lambda row: row["time_code"])
    with (out / "jeans_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    times = [row["time_kyr"] for row in rows]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(times, [row["min_NJ"] for row in rows], "o-", label="Minimum outside sink")
    for n in THRESHOLDS:
        ax.axhline(n, linestyle="--", alpha=.6, label=f"NJ={n}")
    ax.set(xlabel="Time [kyr]", ylabel="Minimum Jeans cells", yscale="log")
    ax.legend()
    ax.grid(alpha=.3)
    fig.savefig(out / "jeans_min_vs_time.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5))
    for n in THRESHOLDS:
        ax.plot(times, [row[f"mass_fraction_NJ_lt_{n}"] for row in rows], "o-", label=f"NJ < {n}")
    ax.set(xlabel="Time [kyr]", ylabel="Underresolved gas mass / selected gas mass", ylim=(0, 1))
    ax.legend()
    ax.grid(alpha=.3)
    fig.savefig(out / "jeans_mass_fraction_vs_time.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out.resolve()}")


def parse_arguments(argv=None):
    # セルに貼り付けた実行では、sys.argvはカーネル起動用の引数。
    # %runではargv[0]が解析スクリプトになるため、通常どおり引数を読む。
    if argv is None and Path(sys.argv[0]).name in ("ipykernel_launcher.py", "ipykernel_launcher"):
        argv = []
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datadir", nargs="?", default=DATADIR)
    parser.add_argument("--output", default=OUTPUT_DIR)
    parser.add_argument("--output-number", type=int, default=2)
    parser.add_argument("--normal", choices=list("xyz"), default="z", help="断面の法線軸")
    parser.add_argument("--slice-au", type=float, default=0.0)
    parser.add_argument("--half-width-au", type=float, default=SLICE_HALF_WIDTH_AU)
    parser.add_argument("--slice-every", type=int, default=1, help="断面画像を何出力ごとに保存するか")
    args = parser.parse_args(argv)
    if args.slice_every < 1 or args.half_width_au <= 0:
        parser.error("slice-everyとhalf-width-auは正値が必要です")
    return args


def main(argv=None):
    args = parse_arguments(argv)
    global np, pv, plt, PolyCollection, LogNorm, Rectangle, Circle
    import numpy as np
    import pyvista as pv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle, Circle
    analyze(args)


if __name__ == "__main__":
    main()


# %% cell2: 現在の密度から必要な追加細分化レベルを見積もる
# cell1を先に実行する。DATADIR・音速・G・シンク半径・単位系はcell1を共有。
# 追加段数は各セルの「現在のレベルからの増分」であり、絶対SMRレベルではない。
# 再計算時の密度変化、MeshBlock単位の細分化、隣接レベル制約は含まない。
TARGET_JEANS_CELLS = 4
REFINEMENT_BOUNDS_AU = ((-10000, 10000),) * 3  # セル中心で選択
REFINEMENT_HALF_WIDTH_AU = 10000.0
REFINEMENT_SLICE_EVERY = 1
REFINEMENT_COLOR_MAX = 8  # 8以上は最上位の色。CSV値には上限をかけない。


def draw_refinement_slice(ax, xyz, delta, nj, selected, args, norm=None, cmap="viridis"):
    normal = "xyz".index(args.normal)
    u, v = [axis for axis in range(3) if axis != normal]
    lower, upper = xyz - delta / 2, xyz + delta / 2
    plane = args.slice_au / L_AU
    # 面上では正側のセルだけを選び、ブロック境界での二重描画を避ける。
    mask = selected & (lower[:, normal] <= plane) & (plane < upper[:, normal])
    half = args.half_width_au / L_AU
    mask &= ((lower[:, u] < half) & (upper[:, u] > -half)
             & (lower[:, v] < half) & (upper[:, v] > -half))
    ids = np.flatnonzero(mask)
    if not len(ids):
        return 0
    polygons = [np.array([[lower[i, u], lower[i, v]],
                          [upper[i, u], lower[i, v]],
                          [upper[i, u], upper[i, v]],
                          [lower[i, u], upper[i, v]]]) * L_AU for i in ids]
    # 元のセルをそのまま描画する。補間は行わない。
    collection = PolyCollection(polygons, array=nj[ids], cmap=cmap,
                                norm=LogNorm(1, 64) if norm is None else norm,
                                edgecolors="none")
    ax.add_collection(collection)
    bounds_lo, bounds_hi = lower.min(axis=0), upper.max(axis=0)
    ax.add_patch(Rectangle((bounds_lo[u]*L_AU, bounds_lo[v]*L_AU),
                          (bounds_hi[u]-bounds_lo[u])*L_AU,
                          (bounds_hi[v]-bounds_lo[v])*L_AU,
                          fill=False, edgecolor="black", linewidth=0.45, alpha=0.6))
    return len(ids)


def estimate_refinement(argv=None):
    from matplotlib.colors import BoundaryNorm

    args = parse_arguments(argv)
    args.half_width_au = REFINEMENT_HALF_WIDTH_AU
    args.slice_every = REFINEMENT_SLICE_EVERY
    target = float(TARGET_JEANS_CELLS)
    if not np.isfinite(target) or target <= 0:
        raise ValueError("TARGET_JEANS_CELLSは正の有限値が必要です")
    if args.slice_every < 1 or args.half_width_au <= 0 or REFINEMENT_COLOR_MAX < 1:
        raise ValueError("cell2の描画設定は正値が必要です")
    if not (np.isfinite(ISO_SOUND_SPEED_CODE) and ISO_SOUND_SPEED_CODE > 0
            and np.isfinite(GRAV_CONST_CODE) and GRAV_CONST_CODE > 0
            and np.isfinite(SINK_RADIUS_AU) and SINK_RADIUS_AU >= 0):
        raise ValueError("音速・重力定数・シンク半径を確認してください")
    region = np.asarray(REFINEMENT_BOUNDS_AU, dtype=float) / L_AU
    if region.shape != (3, 2) or not np.all(np.isfinite(region)) or np.any(region[:, 1] <= region[:, 0]):
        raise ValueError("REFINEMENT_BOUNDS_AUを確認してください")
    groups = group_files(Path(args.datadir).expanduser(), args.output_number)
    out = Path(args.output).expanduser() / "refinement_estimate"
    out.mkdir(parents=True, exist_ok=True)
    (out / "settings.json").write_text(json.dumps(dict(
        arguments=vars(args), target_jeans_cells=target,
        bounds_au=REFINEMENT_BOUNDS_AU, sound_speed_code=ISO_SOUND_SPEED_CODE,
        G_code=GRAV_CONST_CODE, sink_radius_au=SINK_RADIUS_AU,
        L_unit_cgs=L_UNIT_CGS, T_unit_cgs=T_UNIT_CGS,
        assumption="fixed density; additional levels relative to each native cell"
    ), indent=2))
    cmap = plt.get_cmap("viridis", REFINEMENT_COLOR_MAX + 1)
    norm = BoundaryNorm(np.arange(REFINEMENT_COLOR_MAX + 2) - .5, cmap.N)
    history, spatial = [], []
    for index, (step, files) in enumerate(sorted(groups.items())):
        time = vtk_time(next(iter(files.values())))
        draw = index % args.slice_every == 0
        if draw:
            fig, ax = plt.subplots(figsize=(8, 7))
        mass_by_level = defaultdict(float)
        # 各時刻・各追加段数について、該当セルの外縁を囲む箱を記録する。
        boxes = {}
        min_allowed_dx = np.inf
        slice_count = 0
        for path in files.values():
            if not np.isclose(vtk_time(path), time, rtol=1e-10, atol=1e-12):
                raise ValueError(f"時刻不一致: {path}")
            grid = pv.read(path)
            xyz, delta, rho = block_arrays(grid)
            selected = (np.linalg.norm(xyz, axis=1) >= SINK_RADIUS_AU / L_AU)
            selected &= np.all((xyz >= region[:, 0]) & (xyz < region[:, 1]), axis=1)
            jeans = np.sqrt(np.pi * ISO_SOUND_SPEED_CODE**2 /
                            (GRAV_CONST_CODE * rho.astype(float)))
            nj = jeans / delta.max(axis=1)
            extra = np.maximum(0, np.ceil(np.log2(target / nj))).astype(int)
            masses = rho * delta.prod(axis=1)
            if np.any(selected):
                min_allowed_dx = min(min_allowed_dx, float((jeans[selected] / target).min()) * L_AU)
            for level in np.unique(extra[selected]):
                mask = selected & (extra == level)
                mass_by_level[int(level)] += float(masses[mask].sum())
                lo = (xyz[mask] - delta[mask]/2).min(axis=0) * L_AU
                hi = (xyz[mask] + delta[mask]/2).max(axis=0) * L_AU
                if level in boxes:
                    lo = np.minimum(lo, boxes[level][0])
                    hi = np.maximum(hi, boxes[level][1])
                boxes[level] = (lo, hi)
            if draw:
                slice_count += draw_refinement_slice(ax, xyz, delta,
                    np.minimum(extra, REFINEMENT_COLOR_MAX), selected, args, norm, cmap)
            del grid
        total = sum(mass_by_level.values())
        if total <= 0:
            raise ValueError(f"対象領域にシンク外のガスがありません: step={step}")
        history.append(dict(step=step, time_kyr=time*T_KYR,
                            max_extra_levels=max(mass_by_level),
                            min_allowed_dx_au=min_allowed_dx, mass_code=total,
                            fractions={k: v/total for k, v in mass_by_level.items()}))
        for level, (lo, hi) in sorted(boxes.items()):
            spatial.append(dict(step=step, time_kyr=time*T_KYR, extra_levels=int(level),
                                mass_fraction=mass_by_level[int(level)]/total,
                                **{f"{axis}min_au": float(lo[i]) for i, axis in enumerate("xyz")},
                                **{f"{axis}max_au": float(hi[i]) for i, axis in enumerate("xyz")}))
        if draw:
            axes = [a for a in "xyz" if a != args.normal]
            ax.set(xlim=(-args.half_width_au, args.half_width_au),
                   ylim=(-args.half_width_au, args.half_width_au),
                   xlabel=f"{axes[0]} [au]", ylabel=f"{axes[1]} [au]",
                   title=f"Target NJ={target:g}; t={time*T_KYR:.3f} kyr\n"
                         f"{args.normal}={args.slice_au:g} au; black: MeshBlock edges")
            ax.set_aspect("equal")
            ax.set_facecolor("0.85")
            radius2 = SINK_RADIUS_AU**2 - args.slice_au**2
            if radius2 > 0:
                ax.add_patch(Circle((0, 0), np.sqrt(radius2), fill=False, color="red", linestyle="--"))
            if not slice_count:
                ax.text(.5, .5, "No selected cells on this plane", transform=ax.transAxes, ha="center")
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                              ticks=range(REFINEMENT_COLOR_MAX+1))
            cb.set_ticklabels([str(i) for i in range(REFINEMENT_COLOR_MAX)] + [f"{REFINEMENT_COLOR_MAX}+"])
            cb.set_label("Additional refinement levels (fixed density)")
            fig.savefig(out / f"extra_levels_slice_{step:05d}.png", dpi=180, bbox_inches="tight")
            plt.close(fig)
        print(f"[cell2 {index+1}/{len(groups)}] step={step}: "
              f"max additional levels={max(mass_by_level)}, min allowed dx={min_allowed_dx:.3g} au")
    history.sort(key=lambda r: r["time_kyr"])
    max_level = max(r["max_extra_levels"] for r in history)
    rows = [{**{k: v for k, v in r.items() if k != "fractions"},
             **{f"mass_fraction_extra_{k}": r["fractions"].get(k, 0.0)
                for k in range(max_level+1)}} for r in history]
    for name, data in (("extra_levels_history.csv", rows), ("extra_levels_bounds.csv", spatial)):
        with (out / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    times = [r["time_kyr"] for r in history]
    axes[0].step(times, [r["max_extra_levels"] for r in history], where="mid")
    axes[0].set(ylabel="Maximum additional levels", yticks=range(max_level+1),
                title=f"Target NJ={target:g}; fixed-density estimate")
    axes[1].stackplot(times, *[[r["fractions"].get(k, 0.0) for r in history]
                              for k in range(max_level+1)],
                      labels=[f"+{k}" for k in range(max_level+1)])
    axes[1].set(xlabel="Time [kyr]", ylabel="Gas mass fraction by additional levels", ylim=(0, 1))
    axes[1].legend(ncol=4)
    fig.savefig(out / "extra_levels_vs_time.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out.resolve()}")


if __name__ == "__main__":
    estimate_refinement()

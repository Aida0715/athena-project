# -*- coding: utf-8 -*-
"""SMR VTKから円筒円環平均したv_phiプロファイルを作る。

重力成分iの円筒動径加速度g_R,iから
    v_c,i^2 = -R <g_R,i>
を計算する。合計は速度の線形和ではなく、
    v_c,total^2 = -R <g_R,self + g_R,star + g_R,nfw>
である。実際のガス速度から求めた<v_phi>も比較用に描画する。
黒線は内向きのガス自己重力＋中心星重力＋NFW重力と釣り合うために必要なVφである。
つまり紫の点線〜黒線なら概ね回転平衡、黒線＞紫点線なら回転支持不足の可能性と読むことができる。
なお中心星重力曲線は質量加重平均をされるガスの鉛直分布が変化すると上下することがあるが、これは中心星質量が上昇したわけではないので注意。
"""

from collections import defaultdict
from pathlib import Path
import csv
import gc
import json
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# ============================== ユーザー設定 ==============================
DATADIR = Path("~/athena-project/results/○○").expanduser().resolve()
VTK_GLOB = "Toyouchi.block*.out2.*.vtk"
STEPS = None                         # None: 全VTK時刻
SKIP_FIRST_SNAPSHOT = True           # 初期診断量が未更新の最初のVTKを除外
OUTPUT_DIR = DATADIR / "vphi_profiles"
Z_HALF_WIDTH_AU = 2500.0             # Noneでz方向全域
R_MIN_AU = None                      # None: 最小セル幅から自動決定
R_MAX_AU = None                      # None: 完全な円環が入る最大半径
N_RADIAL_BINS = 80
MIN_CELLS_PER_BIN = 8
WEIGHT = "mass"                      # "mass" / "volume"
FIGSIZE = (8.5, 6.2)
DPI = 180
SHOW_IMAGES = False

AU_CGS = 1.495978707e13
KYR_CGS = 1.0e3 * 365.25 * 24.0 * 3600.0


def load_metadata(directory):
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(directory.glob("*.metadata.json"))
    if len(candidates) != 1:
        raise FileNotFoundError(f"Expected exactly one metadata JSON in {directory}; found {len(candidates)}")
    metadata = json.loads(candidates[0].read_text(encoding="utf-8"))
    if metadata.get("format") != "Toyouchi runtime metadata":
        raise ValueError(f"Unsupported metadata: {candidates[0]}")
    metadata["_path"] = str(candidates[0])
    return metadata


def collect_steps():
    groups = defaultdict(list)
    for path in DATADIR.glob(VTK_GLOB):
        match = re.search(r"\.(\d+)\.vtk$", path.name)
        if match:
            step = int(match.group(1))
            if STEPS is None or step in STEPS:
                groups[step].append(path)
    if not groups:
        raise FileNotFoundError(f"No VTK files: {DATADIR / VTK_GLOB}")
    return {step: sorted(paths) for step, paths in sorted(groups.items())}


def read_vtk_time(path):
    with path.open("rb") as handle:
        header = handle.read(4096).decode("ascii", errors="ignore")
    match = re.search(r"\btime\s*=\s*([-+0-9.eEdD]+)", header, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"VTK time was not found: {path}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def block_geometry(grid, path):
    bounds = np.asarray(grid.bounds, dtype=float).reshape(3, 2)
    dimensions = np.asarray(grid.dimensions, dtype=int) - 1
    if np.any(dimensions <= 0) or int(np.prod(dimensions)) != grid.n_cells:
        raise ValueError(f"Unexpected dimensions: {path}; dimensions={grid.dimensions}, n_cells={grid.n_cells}")
    spacing = (bounds[:, 1] - bounds[:, 0]) / dimensions
    if np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError(f"Invalid cell spacing: {path}; {spacing}")
    return bounds, spacing


def make_radial_edges(first_paths, full_bounds, length_au):
    finest_xy = np.inf
    for path in first_paths:
        grid = pv.read(path)
        try:
            _, spacing = block_geometry(grid, path)
            finest_xy = min(finest_xy, spacing[0], spacing[1])
        finally:
            del grid
    auto_rmax = min(full_bounds[0, 1], -full_bounds[0, 0], full_bounds[1, 1], -full_bounds[1, 0])
    # 最内側は1最小セル幅まで表示候補に入れる。ただし、
    # MIN_CELLS_PER_BINを満たさない円環は後段で自動的に無効にする。
    rmin = float(R_MIN_AU) / length_au if R_MIN_AU is not None else finest_xy
    rmax = float(R_MAX_AU) / length_au if R_MAX_AU is not None else auto_rmax
    if not (np.isfinite(rmin) and np.isfinite(rmax) and 0.0 < rmin < rmax):
        raise ValueError(f"Invalid radial range: {rmin=} {rmax=} [code]")
    return np.geomspace(rmin, rmax, int(N_RADIAL_BINS) + 1)


def accumulate_step(paths, edges, z_half_width_code, epsilon_code):
    required = ("rho", "vel", "gr_star", "gr_nfw", "gself_x", "gself_y", "gself_z")
    nbins = edges.size - 1
    denominator = np.zeros(nbins)
    counts = np.zeros(nbins, dtype=np.int64)
    totals = {key: np.zeros(nbins) for key in ("vphi_gas", "gR_self", "gR_star", "gR_nfw")}

    for block_number, path in enumerate(paths, start=1):
        grid = pv.read(path)
        try:
            missing = [name for name in required if name not in grid.cell_data]
            if missing:
                raise KeyError(f"{path}: missing cell arrays {missing}")
            data = grid.cell_data
            rho = np.asarray(data["rho"], dtype=float).reshape(-1)
            velocity = np.asarray(data["vel"], dtype=float)
            xyz = np.asarray(grid.cell_centers().points, dtype=float)
            _, spacing = block_geometry(grid, path)
            dx, dy, dz = spacing
            x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
            radius = np.hypot(x, y)
            spherical_radius = np.sqrt(x*x + y*y + z*z)
            softened_radius = np.sqrt(spherical_radius**2 + epsilon_code**2)

            if z_half_width_code is None:
                dz_overlap = np.full(rho.shape, dz)
            else:
                zlo = np.maximum(z - 0.5 * dz, -z_half_width_code)
                zhi = np.minimum(z + 0.5 * dz, z_half_width_code)
                dz_overlap = np.maximum(zhi - zlo, 0.0)

            vphi = np.divide(-y*velocity[:, 0] + x*velocity[:, 1], radius,
                             out=np.zeros_like(radius), where=radius > 0.0)
            gR_self = np.divide(x*np.asarray(data["gself_x"]) + y*np.asarray(data["gself_y"]), radius,
                                out=np.zeros_like(radius), where=radius > 0.0)
            # gr_starは非softened球半径方向のスカラー値。
            gR_star = np.asarray(data["gr_star"], dtype=float) * np.divide(
                radius, spherical_radius, out=np.zeros_like(radius), where=spherical_radius > 0.0)
            # gr_nfwはsoftened radius方向のスカラー値。
            gR_nfw = np.asarray(data["gr_nfw"], dtype=float) * np.divide(
                radius, softened_radius, out=np.zeros_like(radius), where=softened_radius > 0.0)
            values = (vphi, gR_self, gR_star, gR_nfw)

            finite = (np.isfinite(rho) & (rho > 0.0) & (radius > 0.0)
                      & (dz_overlap > 0.0) & np.all(np.isfinite(velocity), axis=1))
            for value in values:
                finite &= np.isfinite(value)
            ibin = np.searchsorted(edges, radius, side="right") - 1
            use = finite & (ibin >= 0) & (ibin < nbins)
            if not np.any(use):
                continue
            volume = dx * dy * dz_overlap[use]
            weights = volume * rho[use] if WEIGHT == "mass" else volume
            bins = ibin[use]
            denominator += np.bincount(bins, weights=weights, minlength=nbins)
            counts += np.bincount(bins, minlength=nbins)
            for key, value in zip(totals, values):
                totals[key] += np.bincount(bins, weights=weights*value[use], minlength=nbins)
        finally:
            del grid
        if block_number % 128 == 0:
            gc.collect()

    valid_bin = (denominator > 0.0) & (counts >= MIN_CELLS_PER_BIN)
    profiles = {key: np.divide(total, denominator, out=np.full(nbins, np.nan), where=valid_bin)
                for key, total in totals.items()}
    return profiles, counts, denominator


def circular_speed_from_acceleration(radius_code, radial_acceleration):
    """内向き加速度のみ実数の円運動速度へ変換する。"""
    speed_squared = -radius_code * radial_acceleration
    return np.sqrt(speed_squared, out=np.full_like(speed_squared, np.nan),
                   where=np.isfinite(speed_squared) & (speed_squared >= 0.0))


def save_step(step, time_kyr, edges, profiles, counts, weights, units, z_label):
    radius_code = np.sqrt(edges[:-1] * edges[1:])
    radius_au = radius_code * units["length_au"]
    gR_total = profiles["gR_self"] + profiles["gR_star"] + profiles["gR_nfw"]
    curves_code = {
        "self": circular_speed_from_acceleration(radius_code, profiles["gR_self"]),
        "star": circular_speed_from_acceleration(radius_code, profiles["gR_star"]),
        "nfw": circular_speed_from_acceleration(radius_code, profiles["gR_nfw"]),
        "total": circular_speed_from_acceleration(radius_code, gR_total),
        "gas": profiles["vphi_gas"],
    }
    curves = {key: value * units["velocity_kms"] for key, value in curves_code.items()}
    png = OUTPUT_DIR / f"vphi_profile_{step:05d}.png"
    csv_path = OUTPUT_DIR / f"vphi_profile_{step:05d}.csv"

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(radius_au, curves["self"], lw=1.5, label="Gas self-gravity")
    ax.plot(radius_au, curves["star"], lw=1.5, label="Central star")
    ax.plot(radius_au, curves["nfw"], lw=1.5, label="NFW halo")
    ax.plot(radius_au, curves["total"], color="black", lw=2.2,
            label="Total gravitational circular speed")
    ax.plot(radius_au, curves["gas"], color="tab:purple", ls="--", lw=1.7,
            label=r"Measured gas $\langle v_\phi\rangle$")
    ax.set_xscale("log")
    ax.set_xlabel("Cylindrical radius R [AU]")
    ax.set_ylabel(r"$v_\phi$ [km s$^{-1}$]")
    ax.set_title(
        rf"Azimuthal-velocity contributions, step={step:05d}, "
        rf"$t={time_kyr:.1f}$ kyr" + "\n"
        + f"{WEIGHT}-weighted, {z_label}"
    )
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=DPI)
    if SHOW_IMAGES:
        plt.show()
    plt.close(fig)

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["radius_AU", "cell_count", "weight_sum_code",
                         "vphi_self_km_s", "vphi_star_km_s", "vphi_nfw_km_s",
                         "vphi_total_gravity_km_s", "vphi_gas_mass_weighted_km_s",
                         "gR_self_code", "gR_star_code", "gR_nfw_code", "gR_total_code"])
        for i in range(radius_au.size):
            writer.writerow([f"{radius_au[i]:.10e}", int(counts[i]), f"{weights[i]:.10e}",
                             f"{curves['self'][i]:.10e}", f"{curves['star'][i]:.10e}",
                             f"{curves['nfw'][i]:.10e}", f"{curves['total'][i]:.10e}",
                             f"{curves['gas'][i]:.10e}", f"{profiles['gR_self'][i]:.10e}",
                             f"{profiles['gR_star'][i]:.10e}", f"{profiles['gR_nfw'][i]:.10e}",
                             f"{gR_total[i]:.10e}"])
    return png, csv_path


def main():
    if not DATADIR.is_dir():
        raise FileNotFoundError(DATADIR)
    if WEIGHT not in ("mass", "volume") or N_RADIAL_BINS < 2:
        raise ValueError("Invalid WEIGHT or N_RADIAL_BINS")
    metadata = load_metadata(DATADIR)
    unit_data = metadata["units"]
    length_au = float(unit_data["length_cgs"]) / AU_CGS
    units = {"length_au": length_au,
             "velocity_kms": float(unit_data["velocity_cgs"]) / 1.0e5,
             "time_cgs": float(unit_data["time_cgs"])}
    epsilon_code = float(metadata["gravity"]["softening_code"])
    mesh = metadata["mesh"]
    full_bounds = np.array([[mesh["x1min_code"], mesh["x1max_code"]],
                            [mesh["x2min_code"], mesh["x2max_code"]],
                            [mesh["x3min_code"], mesh["x3max_code"]]], dtype=float)
    z_half_width_code = None if Z_HALF_WIDTH_AU is None else float(Z_HALF_WIDTH_AU)/length_au
    z_label = "full z domain" if Z_HALF_WIDTH_AU is None else rf"$|z|\leq {Z_HALF_WIDTH_AU:g}$ AU"
    groups = collect_steps()
    if SKIP_FIRST_SNAPSHOT:
        first_step = next(iter(groups))
        print(
            f"[INFO] Skip     : first snapshot step={first_step:05d} "
            "(gravity diagnostics may be uninitialized)"
        )
        groups = dict(list(groups.items())[1:])
        if not groups:
            raise RuntimeError(
                "No VTK snapshots remain after SKIP_FIRST_SNAPSHOT=True"
            )
    edges = make_radial_edges(next(iter(groups.values())), full_bounds, length_au)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Metadata : {metadata['_path']}")
    print(f"[INFO] Input    : {DATADIR}")
    print(f"[INFO] Output   : {OUTPUT_DIR}")
    print(f"[INFO] Steps    : {len(groups)}")
    print(f"[INFO] z range  : {z_label}")

    for number, (step, paths) in enumerate(groups.items(), start=1):
        time_code = read_vtk_time(paths[0])
        time_kyr = time_code * units["time_cgs"] / KYR_CGS
        print(f"[RUN {number:3d}/{len(groups):3d}] step={step:05d}, t={time_kyr:.1f} kyr, blocks={len(paths)}")
        profiles, counts, weights = accumulate_step(paths, edges, z_half_width_code, epsilon_code)
        if not np.any(np.isfinite(profiles["gR_self"])):
            warnings.warn(f"step={step:05d}: no valid gravity profile")
        png, csv_path = save_step(step, time_kyr, edges, profiles, counts, weights, units, z_label)
        print(f"      [SAVED] {png.name}, {csv_path.name}")
        del profiles, counts, weights
        gc.collect()
    print(f"[DONE] {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

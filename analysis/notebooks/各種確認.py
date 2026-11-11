# -*- coding: utf-8 -*-
"""Toyouchi: シンクなしの現行計算向け HST/VTK 健全性確認。"""

# %% cell1: HST/VTK健全性確認
from pathlib import Path
from collections import defaultdict
import csv, gc, re
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
try:
    from toyouchi_metadata import load_toyouchi_metadata, run_units
except ModuleNotFoundError:
    import json
    def load_toyouchi_metadata(data_path):
        directory = Path(data_path).expanduser().resolve()
        directory = directory if directory.is_dir() else directory.parent
        preferred = directory / "Toyouchi.metadata.json"
        candidates = [preferred] if preferred.is_file() else sorted(directory.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError(f"Expected one metadata JSON in {directory}; found {len(candidates)}")
        metadata = json.loads(candidates[0].read_text(encoding="utf-8")); metadata["_path"] = str(candidates[0]); return metadata
    def run_units(metadata):
        keys = ("mass_cgs", "length_cgs", "time_cgs", "velocity_cgs", "density_cgs", "pressure_cgs", "magnetic_field_gauss", "acceleration_cgs", "gravity_constant_code")
        return {key: float(metadata["units"][key]) for key in keys}

# ======================== ユーザー設定 ========================
DATA_DIR = Path("~/athena-project/results/○○").expanduser()
HST_FILENAME = "Toyouchi.hst"
VTK_GLOB = "**/*.block*.out2.*.vtk"
OUTPUT_DIR = DATA_DIR / "checks"
SHOW_PLOTS = True

# 現行の基準単位
RUN_UNITS = run_units(load_toyouchi_metadata(DATA_DIR))
M_UNIT_G = RUN_UNITS["mass_cgs"]
L_UNIT_CM = RUN_UNITS["length_cgs"]
T_UNIT_S = RUN_UNITS["time_cgs"]
MSUN_G = 1.98847e33
TIME_KYR = T_UNIT_S / (1.e3 * 365.25 * 86400.)
MASS_MSUN = M_UNIT_G / MSUN_G
V_UNIT = L_UNIT_CM / T_UNIT_S
RHO_UNIT = M_UNIT_G / L_UNIT_CM**3
B_UNIT = RUN_UNITS["magnetic_field_gauss"]
E_UNIT = M_UNIT_G * V_UNIT**2
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def hst_columns(path):
    result = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.lstrip().startswith("#"):
                for number, name in re.findall(r"\[(\d+)\]\s*=\s*([^\s]+)", line):
                    result[name] = int(number)-1
    return result

def load_hst(path):
    columns = hst_columns(path)
    if "time" not in columns:
        raise KeyError("HST header has no time column")
    data = np.loadtxt(path, comments="#", ndmin=2)
    data = data[np.isfinite(data[:, columns["time"]])]
    data = data[np.argsort(data[:, columns["time"]], kind="stable")]
    _, reverse_index = np.unique(data[::-1, columns["time"]], return_index=True)
    return columns, data[np.sort(len(data)-1-reverse_index)]

def plot_hst(columns, data):
    time = data[:, columns["time"]] * TIME_KYR
    series = []
    if "mass" in columns:
        series.append(("Gas mass", data[:, columns["mass"]]*MASS_MSUN, r"$M_{\rm gas}\ [M_\odot]$"))
    if "Mstar" in columns:  # 固定中心星は現行コードに存在する
        series.append(("Central-star mass", data[:, columns["Mstar"]]*MASS_MSUN, r"$M_\star\ [M_\odot]$"))
    if not series:
        return
    fig, axes = plt.subplots(len(series), 1, figsize=(8, 3.2*len(series)), sharex=True, squeeze=False)
    for ax, (title, values, ylabel) in zip(axes[:, 0], series):
        ax.plot(time, values, lw=1.5); ax.set(title=title, ylabel=ylabel); ax.grid(alpha=.3)
    axes[-1, 0].set_xlabel("Time [kyr]")
    fig.tight_layout(); fig.savefig(OUTPUT_DIR/"history_overview.png", dpi=180)
    if SHOW_PLOTS: plt.show()
    plt.close(fig)

def vtk_step(path):
    match = re.search(r"\.out\d+\.(\d+)\.vtk$", path.name)
    return int(match.group(1)) if match else None

def grids(dataset):
    if isinstance(dataset, pv.MultiBlock):
        for block in dataset:
            if block is not None: yield from grids(block)
    else:
        yield dataset

def analyze_step(step, paths):
    out = defaultdict(float)
    out.update(step=step, blocks=0, cells=0, nonfinite=0, nonpositive_rho=0,
               rho_min=np.inf, rho_max=-np.inf, b_max=0., time=np.nan)
    times = []
    for path in paths:
        dataset = pv.read(path)
        for grid in grids(dataset):
            sized = grid.compute_cell_sizes(length=False, area=False, volume=True)
            vol = np.asarray(sized.cell_data["Volume"], dtype=float)
            rho = np.asarray(grid.cell_data["rho"], dtype=float).reshape(-1)
            out["blocks"] += 1; out["cells"] += rho.size
            out["nonfinite"] += np.count_nonzero(~np.isfinite(rho))
            out["nonpositive_rho"] += np.count_nonzero(rho <= 0.)
            good = np.isfinite(rho) & np.isfinite(vol) & (vol > 0.)
            if np.any(good):
                out["mass"] += np.sum(rho[good]*vol[good]); out["volume"] += np.sum(vol[good])
                out["rho_min"] = min(out["rho_min"], np.min(rho[good])); out["rho_max"] = max(out["rho_max"], np.max(rho[good]))
            if "Bcc" in grid.cell_data:
                b = np.asarray(grid.cell_data["Bcc"], dtype=float).reshape(-1, 3)
                b2 = np.einsum("ij,ij->i", b, b); bgood = good & np.isfinite(b2)
                if np.any(bgood):
                    out["magnetic_energy"] += np.sum(.5*b2[bgood]*vol[bgood]); out["b_max"] = max(out["b_max"], np.sqrt(np.max(b2[bgood])))
            if "TIME" in grid.field_data: times.extend(np.ravel(grid.field_data["TIME"]).astype(float))
            del sized
        del dataset; gc.collect()
    if times: out["time"] = float(np.median(times))
    return dict(out)

def save_vtk(rows):
    if not rows:
        print("[WARNING] VTK files not found; HST check only."); return
    fields = ["step","time_kyr","blocks","cells","volume_cm3","mass_msun","magnetic_energy_erg","rho_min_gcm3","rho_max_gcm3","b_max_gauss","nonfinite","nonpositive_rho"]
    with (OUTPUT_DIR/"vtk_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for r in rows:
            writer.writerow(dict(step=r["step"], time_kyr=r["time"]*TIME_KYR, blocks=int(r["blocks"]), cells=int(r["cells"]), volume_cm3=r["volume"]*L_UNIT_CM**3, mass_msun=r["mass"]*MASS_MSUN, magnetic_energy_erg=r["magnetic_energy"]*E_UNIT, rho_min_gcm3=r["rho_min"]*RHO_UNIT, rho_max_gcm3=r["rho_max"]*RHO_UNIT, b_max_gauss=r["b_max"]*B_UNIT, nonfinite=int(r["nonfinite"]), nonpositive_rho=int(r["nonpositive_rho"])))
    t = np.array([r["time"] for r in rows])*TIME_KYR
    values = [(np.array([r["mass"] for r in rows])*MASS_MSUN, r"Gas mass [$M_\odot$]", False), (np.array([r["magnetic_energy"] for r in rows])*E_UNIT, "Magnetic energy [erg]", True), (np.array([r["rho_min"] for r in rows])*RHO_UNIT, r"Minimum density [g cm$^{-3}$]", True), (np.array([r["rho_max"] for r in rows])*RHO_UNIT, r"Maximum density [g cm$^{-3}$]", True)]
    fig, axes = plt.subplots(2,2,figsize=(10,7),sharex=True)
    for ax,(y,label,log) in zip(axes.flat, values):
        ax.plot(t,y,marker="o",ms=2); ax.set_ylabel(label); ax.grid(alpha=.3)
        if log and np.any(y>0): ax.set_yscale("log")
    axes[1,0].set_xlabel("Time [kyr]"); axes[1,1].set_xlabel("Time [kyr]")
    fig.tight_layout(); fig.savefig(OUTPUT_DIR/"vtk_overview.png",dpi=180)
    if SHOW_PLOTS: plt.show()
    plt.close(fig)

hst_path = DATA_DIR/HST_FILENAME
if hst_path.is_file():
    cmap, table = load_hst(hst_path); plot_hst(cmap, table); print(f"[OK] HST rows: {len(table)}")
else: print(f"[WARNING] HST file not found: {hst_path}")
groups = defaultdict(list)
for path in sorted(DATA_DIR.glob(VTK_GLOB)):
    step = vtk_step(path)
    if step is not None: groups[step].append(path)
rows = []
for index,(step,paths) in enumerate(sorted(groups.items()),1):
    print(f"[VTK {index}/{len(groups)}] step={step:05d}, files={len(paths)}")
    rows.append(analyze_step(step,paths))
save_vtk(rows)
print(f"[SAVED] {OUTPUT_DIR}")

# %% cell2: 全VTK時刻のxy速度カラーマップ（full / Zoom）
# このセルはcell1を実行せず単独で実行できます。
# 速度は面内成分ではなく sqrt(vx^2+vy^2+vz^2) [km/s]。
# RectilinearGridのセル値をそのまま描画し、SMR境界で補間しません。
from pathlib import Path
from collections import defaultdict
import json
import re
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# -------- このセル専用のユーザー設定 --------
SPEED_DATA_DIR = Path("~/athena-project/results/○○").expanduser()
SPEED_VTK_GLOB = "**/*.block*.out2.*.vtk"
SPEED_OUTPUT_DIR = SPEED_DATA_DIR / "checks" / "speed_maps"
SPEED_ZOOM_FRACTION = 0.1  # データ全域の幅の1/10、原点中心
SPEED_SLICE_AU = 0.0             # xy: z位置 / xz: y位置
# fullは代表時刻で共通、Zoomは各時刻の表示領域から色範囲を決める。
SPEED_PERCENTILES = (1.0, 99.0)
SPEED_SCALE_SAMPLE_COUNT = 12
SPEED_CMAP = "magma"
SPEED_DPI = 180
SPEED_SHOW = False              # 大量出力時はFalse推奨


def speed_map_groups():
    groups = defaultdict(list)
    seen = set()
    for path in sorted(SPEED_DATA_DIR.glob(SPEED_VTK_GLOB)):
        match = re.fullmatch(r"(.+)\.block(\d+)\.out(\d+)\.(\d+)\.vtk", path.name)
        if not match:
            continue
        prefix, block, stream, step = match.groups()
        key = (prefix, int(stream), int(step))
        identity = key + (int(block),)
        if identity in seen:
            raise ValueError(f"同じ出力時刻・ブロックが重複しています: {path}")
        seen.add(identity)
        groups[key].append(path)
    if not groups:
        raise FileNotFoundError(f"VTKがありません: {SPEED_DATA_DIR}")
    return sorted(groups.items())


def speed_map_time(path):
    with path.open("rb") as f:
        f.readline()
        header = f.readline().decode("ascii", errors="replace")
    match = re.search(r"time\s*=\s*([-+\d.eE]+)", header)
    if not match:
        raise ValueError(f"VTKヘッダーに時刻がありません: {path}")
    return float(match.group(1))


def speed_map_tiles(paths, plane, length_au, velocity_kms):
    """各ブロックの交差セル層を返す。境界上では正側を選び二重描画を防ぐ。"""
    normal = 2 if plane == "xy" else 1
    cut = SPEED_SLICE_AU / length_au
    for path in paths:
        grid = pv.read(path)
        if not isinstance(grid, pv.RectilinearGrid):
            raise TypeError(f"RectilinearGrid以外は未対応: {path}")
        edges = [np.asarray(a, dtype=float) for a in (grid.x, grid.y, grid.z)]
        axis = edges[normal]
        if not axis[0] <= cut < axis[-1]:
            continue
        if any(np.any(np.diff(a) <= 0) for a in edges):
            raise ValueError(f"座標が単調増加ではありません: {path}")
        layer = int(np.searchsorted(axis, cut, side="right") - 1)
        nx, ny, nz = (len(a)-1 for a in edges)
        vel = np.asarray(grid.cell_data["vel"], dtype=float)
        if vel.shape != (nx*ny*nz, 3):
            raise ValueError(f"velの形状が不正です: {path}: {vel.shape}")
        speed = np.linalg.norm(vel, axis=1).reshape(nz, ny, nx) * velocity_kms
        values = speed[layer, :, :] if plane == "xy" else speed[:, layer, :]
        if not np.isfinite(values).all():
            print(f"[WARNING] 非有限速度を白抜き: {path.name}")
        yield (edges[0]*length_au, edges[1 if plane == "xy" else 2]*length_au,
               np.ma.masked_invalid(values))


def speed_map_norm(values, fallback=None):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return fallback if fallback is not None else Normalize(0, 1)
    lo, hi = np.percentile(values, SPEED_PERCENTILES)
    if hi <= lo:
        pad = max(abs(float(lo))*0.01, 1.e-6)
        lo, hi = max(0.0, lo-pad), hi+pad
    return Normalize(float(lo), float(hi), clip=False)


def render_speed_maps(plane):
    if plane not in ("xy", "xz"):
        raise ValueError(plane)
    if not 0 < SPEED_ZOOM_FRACTION <= 1:
        raise ValueError("Zoom倍率は0より大きく1以下にしてください")
    if not 0 <= SPEED_PERCENTILES[0] < SPEED_PERCENTILES[1] <= 100:
        raise ValueError("パーセンタイル範囲が不正です")
    if SPEED_SCALE_SAMPLE_COUNT < 1:
        raise ValueError("色範囲のサンプル時刻数は1以上にしてください")
    metadata_path = SPEED_DATA_DIR / "Toyouchi.metadata.json"
    if not metadata_path.is_file():
        candidates = sorted(SPEED_DATA_DIR.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError("DATA_DIRに単位系のmetadata JSONを置いてください")
        metadata_path = candidates[0]
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    units = metadata["units"]
    length_au = float(units["length_cgs"]) / 1.495978707e13
    velocity_kms = float(units["velocity_cgs"]) / 1.e5
    time_kyr = float(units["time_cgs"]) / (365.25*86400*1000)
    groups = speed_map_groups()
    samples = []
    indices = np.unique(np.linspace(0, len(groups)-1,
                        min(SPEED_SCALE_SAMPLE_COUNT, len(groups))).astype(int))
    for index in indices:
        _, paths = groups[index]
        print(f"[{plane} 色範囲サンプル] {index+1}/{len(groups)}")
        for _, _, values in speed_map_tiles(paths, plane, length_au, velocity_kms):
            samples.append(values.compressed())
    norm = speed_map_norm(np.concatenate(samples) if samples else [])
    del samples
    cmap = plt.get_cmap(SPEED_CMAP).copy()
    cmap.set_bad("white")
    out = SPEED_OUTPUT_DIR / plane
    for view in ("full", "Zoom"):
        (out / view).mkdir(parents=True, exist_ok=True)
    for number, ((prefix, stream, step), paths) in enumerate(groups, 1):
        times = [speed_map_time(p) for p in paths]
        if not np.allclose(times, times[0], rtol=1.e-7, atol=1.e-12):
            raise ValueError(f"同じ番号のVTKで時刻が不一致: {step}")
        figures = [plt.subplots(figsize=(8, 7)) for _ in range(2)]
        bounds = [np.inf, -np.inf, np.inf, -np.inf]
        tiles = 0
        drawn_tiles = []
        zoom_meshes = []
        raw_max = 0.0
        try:
            for x, vertical, values in speed_map_tiles(paths, plane, length_au, velocity_kms):
                tiles += 1
                drawn_tiles.append((x, vertical, values))
                if values.count():
                    raw_max = max(raw_max, float(values.max()))
                bounds = [min(bounds[0], x[0]), max(bounds[1], x[-1]),
                          min(bounds[2], vertical[0]), max(bounds[3], vertical[-1])]
                for view_index, (_, ax) in enumerate(figures):
                    mesh = ax.pcolormesh(x, vertical, values, shading="flat", cmap=cmap,
                                        norm=norm, rasterized=True, edgecolors="none")
                    if view_index == 1:
                        zoom_meshes.append(mesh)
            if not tiles:
                raise ValueError(f"指定断面に交差するブロックなし: {plane}, step={step}")
            xhalf = (bounds[1]-bounds[0])*0.5*SPEED_ZOOM_FRACTION
            yhalf = (bounds[3]-bounds[2])*0.5*SPEED_ZOOM_FRACTION
            zoom_values = []
            for x, vertical, values in drawn_tiles:
                # 表示範囲に重なるセルを採用（セル中心が外でも端が見える場合を含む）。
                inside = ((vertical[:-1, None] < yhalf)
                          & (vertical[1:, None] > -yhalf)
                          & (x[None, :-1] < xhalf) & (x[None, 1:] > -xhalf))
                zoom_values.append(values[inside].compressed())
            zoom_norm = speed_map_norm(np.concatenate(zoom_values), fallback=norm)
            for mesh in zoom_meshes:
                mesh.set_norm(zoom_norm)
            for view, (fig, ax) in zip(("full", "Zoom"), figures):
                if view == "full":
                    ax.set_xlim(bounds[:2]); ax.set_ylim(bounds[2:])
                else:
                    ax.set_xlim(-xhalf, xhalf); ax.set_ylim(-yhalf, yhalf)
                ax.set_aspect("equal")
                ax.set_xlabel("x [au]"); ax.set_ylabel(f"{plane[1]} [au]")
                ax.set_title(f"{plane} |v| — {times[0]*time_kyr:.6g} kyr ({view})\n"
                             f"{'z' if plane == 'xy' else 'y'} = {SPEED_SLICE_AU:g} au")
                view_norm = norm if view == "full" else zoom_norm
                fig.colorbar(plt.cm.ScalarMappable(norm=view_norm, cmap=cmap), ax=ax,
                             extend="both", label="Speed |v| [km/s]")
                fig.tight_layout()
                fig.savefig(out / view / f"{prefix}.out{stream}.{step:05d}.png", dpi=SPEED_DPI)
                if SPEED_SHOW:
                    plt.show()
            print(f"[{plane} {number}/{len(groups)}] step={step:05d}, slice blocks={tiles}, "
                  f"実最大速度={raw_max:.6g} km/s, "
                  f"Zoom色範囲={zoom_norm.vmin:.6g}〜{zoom_norm.vmax:.6g} km/s")
        finally:
            for fig, _ in figures:
                plt.close(fig)
    print(f"[SAVED] {out} (vmin={norm.vmin:g}, vmax={norm.vmax:g} km/s)")


render_speed_maps("xy")

# %% cell3: 全VTK時刻のxz速度カラーマップ（full / Zoom）
# このセルはcell1・cell2を実行せず単独で実行できます。
# 速度は面内成分ではなく sqrt(vx^2+vy^2+vz^2) [km/s]。
# RectilinearGridのセル値をそのまま描画し、SMR境界で補間しません。
from pathlib import Path
from collections import defaultdict
import json
import re
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# -------- このセル専用のユーザー設定 --------
XZ_SPEED_DATA_DIR = Path("~/athena-project/results/○○").expanduser()
XZ_SPEED_VTK_GLOB = "**/*.block*.out2.*.vtk"
XZ_SPEED_OUTPUT_DIR = XZ_SPEED_DATA_DIR / "checks" / "speed_maps"
XZ_SPEED_ZOOM_FRACTION = 0.1  # データ全域の幅の1/10、原点中心
XZ_SPEED_SLICE_AU = 0.0             # xy: z位置 / xz: y位置
# fullは代表時刻で共通、Zoomは各時刻の表示領域から色範囲を決める。
XZ_SPEED_PERCENTILES = (1.0, 99.0)
XZ_SPEED_SCALE_SAMPLE_COUNT = 12
XZ_SPEED_CMAP = "magma"
XZ_SPEED_DPI = 180
XZ_SPEED_SHOW = False              # 大量出力時はFalse推奨


def xz_speed_map_groups():
    groups = defaultdict(list)
    seen = set()
    for path in sorted(XZ_SPEED_DATA_DIR.glob(XZ_SPEED_VTK_GLOB)):
        match = re.fullmatch(r"(.+)\.block(\d+)\.out(\d+)\.(\d+)\.vtk", path.name)
        if not match:
            continue
        prefix, block, stream, step = match.groups()
        key = (prefix, int(stream), int(step))
        identity = key + (int(block),)
        if identity in seen:
            raise ValueError(f"同じ出力時刻・ブロックが重複しています: {path}")
        seen.add(identity)
        groups[key].append(path)
    if not groups:
        raise FileNotFoundError(f"VTKがありません: {XZ_SPEED_DATA_DIR}")
    return sorted(groups.items())


def xz_speed_map_time(path):
    with path.open("rb") as f:
        f.readline()
        header = f.readline().decode("ascii", errors="replace")
    match = re.search(r"time\s*=\s*([-+\d.eE]+)", header)
    if not match:
        raise ValueError(f"VTKヘッダーに時刻がありません: {path}")
    return float(match.group(1))


def xz_speed_map_tiles(paths, plane, length_au, velocity_kms):
    """各ブロックの交差セル層を返す。境界上では正側を選び二重描画を防ぐ。"""
    normal = 2 if plane == "xy" else 1
    cut = XZ_SPEED_SLICE_AU / length_au
    for path in paths:
        grid = pv.read(path)
        if not isinstance(grid, pv.RectilinearGrid):
            raise TypeError(f"RectilinearGrid以外は未対応: {path}")
        edges = [np.asarray(a, dtype=float) for a in (grid.x, grid.y, grid.z)]
        axis = edges[normal]
        if not axis[0] <= cut < axis[-1]:
            continue
        if any(np.any(np.diff(a) <= 0) for a in edges):
            raise ValueError(f"座標が単調増加ではありません: {path}")
        layer = int(np.searchsorted(axis, cut, side="right") - 1)
        nx, ny, nz = (len(a)-1 for a in edges)
        vel = np.asarray(grid.cell_data["vel"], dtype=float)
        if vel.shape != (nx*ny*nz, 3):
            raise ValueError(f"velの形状が不正です: {path}: {vel.shape}")
        speed = np.linalg.norm(vel, axis=1).reshape(nz, ny, nx) * velocity_kms
        values = speed[layer, :, :] if plane == "xy" else speed[:, layer, :]
        if not np.isfinite(values).all():
            print(f"[WARNING] 非有限速度を白抜き: {path.name}")
        yield (edges[0]*length_au, edges[1 if plane == "xy" else 2]*length_au,
               np.ma.masked_invalid(values))


def xz_speed_map_norm(values, fallback=None):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return fallback if fallback is not None else Normalize(0, 1)
    lo, hi = np.percentile(values, XZ_SPEED_PERCENTILES)
    if hi <= lo:
        pad = max(abs(float(lo))*0.01, 1.e-6)
        lo, hi = max(0.0, lo-pad), hi+pad
    return Normalize(float(lo), float(hi), clip=False)


def render_xz_speed_maps(plane):
    if plane not in ("xy", "xz"):
        raise ValueError(plane)
    if not 0 < XZ_SPEED_ZOOM_FRACTION <= 1:
        raise ValueError("Zoom倍率は0より大きく1以下にしてください")
    if not 0 <= XZ_SPEED_PERCENTILES[0] < XZ_SPEED_PERCENTILES[1] <= 100:
        raise ValueError("パーセンタイル範囲が不正です")
    if XZ_SPEED_SCALE_SAMPLE_COUNT < 1:
        raise ValueError("色範囲のサンプル時刻数は1以上にしてください")
    metadata_path = XZ_SPEED_DATA_DIR / "Toyouchi.metadata.json"
    if not metadata_path.is_file():
        candidates = sorted(XZ_SPEED_DATA_DIR.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError("DATA_DIRに単位系のmetadata JSONを置いてください")
        metadata_path = candidates[0]
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    units = metadata["units"]
    length_au = float(units["length_cgs"]) / 1.495978707e13
    velocity_kms = float(units["velocity_cgs"]) / 1.e5
    time_kyr = float(units["time_cgs"]) / (365.25*86400*1000)
    groups = xz_speed_map_groups()
    samples = []
    indices = np.unique(np.linspace(0, len(groups)-1,
                        min(XZ_SPEED_SCALE_SAMPLE_COUNT, len(groups))).astype(int))
    for index in indices:
        _, paths = groups[index]
        print(f"[{plane} 色範囲サンプル] {index+1}/{len(groups)}")
        for _, _, values in xz_speed_map_tiles(paths, plane, length_au, velocity_kms):
            samples.append(values.compressed())
    norm = xz_speed_map_norm(np.concatenate(samples) if samples else [])
    del samples
    cmap = plt.get_cmap(XZ_SPEED_CMAP).copy()
    cmap.set_bad("white")
    out = XZ_SPEED_OUTPUT_DIR / plane
    for view in ("full", "Zoom"):
        (out / view).mkdir(parents=True, exist_ok=True)
    for number, ((prefix, stream, step), paths) in enumerate(groups, 1):
        times = [xz_speed_map_time(p) for p in paths]
        if not np.allclose(times, times[0], rtol=1.e-7, atol=1.e-12):
            raise ValueError(f"同じ番号のVTKで時刻が不一致: {step}")
        figures = [plt.subplots(figsize=(8, 7)) for _ in range(2)]
        bounds = [np.inf, -np.inf, np.inf, -np.inf]
        tiles = 0
        drawn_tiles = []
        zoom_meshes = []
        raw_max = 0.0
        try:
            for x, vertical, values in xz_speed_map_tiles(paths, plane, length_au, velocity_kms):
                tiles += 1
                drawn_tiles.append((x, vertical, values))
                if values.count():
                    raw_max = max(raw_max, float(values.max()))
                bounds = [min(bounds[0], x[0]), max(bounds[1], x[-1]),
                          min(bounds[2], vertical[0]), max(bounds[3], vertical[-1])]
                for view_index, (_, ax) in enumerate(figures):
                    mesh = ax.pcolormesh(x, vertical, values, shading="flat", cmap=cmap,
                                        norm=norm, rasterized=True, edgecolors="none")
                    if view_index == 1:
                        zoom_meshes.append(mesh)
            if not tiles:
                raise ValueError(f"指定断面に交差するブロックなし: {plane}, step={step}")
            xhalf = (bounds[1]-bounds[0])*0.5*XZ_SPEED_ZOOM_FRACTION
            yhalf = (bounds[3]-bounds[2])*0.5*XZ_SPEED_ZOOM_FRACTION
            zoom_values = []
            for x, vertical, values in drawn_tiles:
                # 表示範囲に重なるセルを採用（セル中心が外でも端が見える場合を含む）。
                inside = ((vertical[:-1, None] < yhalf)
                          & (vertical[1:, None] > -yhalf)
                          & (x[None, :-1] < xhalf) & (x[None, 1:] > -xhalf))
                zoom_values.append(values[inside].compressed())
            zoom_norm = xz_speed_map_norm(np.concatenate(zoom_values), fallback=norm)
            for mesh in zoom_meshes:
                mesh.set_norm(zoom_norm)
            for view, (fig, ax) in zip(("full", "Zoom"), figures):
                if view == "full":
                    ax.set_xlim(bounds[:2]); ax.set_ylim(bounds[2:])
                else:
                    ax.set_xlim(-xhalf, xhalf); ax.set_ylim(-yhalf, yhalf)
                ax.set_aspect("equal")
                ax.set_xlabel("x [au]"); ax.set_ylabel(f"{plane[1]} [au]")
                ax.set_title(f"{plane} |v| — {times[0]*time_kyr:.6g} kyr ({view})\n"
                             f"{'z' if plane == 'xy' else 'y'} = {XZ_SPEED_SLICE_AU:g} au")
                view_norm = norm if view == "full" else zoom_norm
                fig.colorbar(plt.cm.ScalarMappable(norm=view_norm, cmap=cmap), ax=ax,
                             extend="both", label="Speed |v| [km/s]")
                fig.tight_layout()
                fig.savefig(out / view / f"{prefix}.out{stream}.{step:05d}.png", dpi=XZ_SPEED_DPI)
                if XZ_SPEED_SHOW:
                    plt.show()
            print(f"[{plane} {number}/{len(groups)}] step={step:05d}, slice blocks={tiles}, "
                  f"実最大速度={raw_max:.6g} km/s, "
                  f"Zoom色範囲={zoom_norm.vmin:.6g}〜{zoom_norm.vmax:.6g} km/s")
        finally:
            for fig, _ in figures:
                plt.close(fig)
    print(f"[SAVED] {out} (vmin={norm.vmin:g}, vmax={norm.vmax:g} km/s)")


render_xz_speed_maps("xz")


# %% cell4: xy最大速度投影（full / Zoom、独立実行）
# 中央断面ではなく視線方向の最大|v|。3次元速度の大きさを使用。
# ピクセルと重なるセルの最大値を採用：微小な高速セルも平均で消さない。
# そのためピクセル未満の構造は画面上で広がる場合があります。
from pathlib import Path
from collections import defaultdict
import csv
import json
import re
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap

MAX_XY_DATA_DIR = Path("~/athena-project/results/○○").expanduser()
MAX_XY_VTK_GLOB = "**/*.block*.out2.*.vtk"
MAX_XY_OUTPUT_DIR = MAX_XY_DATA_DIR / "checks" / "speed_max_projection" / "xy"
MAX_XY_THRESHOLD_KMS = 100.0  # この値以上を赤で強調（手動変更可）
MAX_XY_ZOOM_FRACTION = 0.1    # 原点中心、全域幅の1/10
MAX_XY_PIXELS = 800          # full / Zoomそれぞれの一辺の画素数
MAX_XY_DPI = 180
MAX_XY_SHOW = False
MAX_XY_COLOR_MIN_KMS = None  # None: 表示画像の正の速度の1パーセンタイル
MAX_XY_COLOR_MAX_KMS = None  # None: 各時刻の全領域最大速度（外れ値も含む）


def max_xy_raster_max(tiles, bounds, pixels):
    """セル矩形の値を重なる全ピクセルへmaxで集約。欠損はNaN。"""
    xmin, xmax, ymin, ymax = bounds
    xe = np.linspace(xmin, xmax, pixels+1)
    ye = np.linspace(ymin, ymax, pixels+1)
    image = np.full((pixels, pixels), np.nan)
    for x, y, values in tiles:
        for j in range(len(y)-1):
            if y[j+1] <= ymin or y[j] >= ymax:
                continue
            j0 = max(0, np.searchsorted(ye, y[j], side="right")-1)
            j1 = min(pixels, np.searchsorted(ye, y[j+1], side="left"))
            for i in range(len(x)-1):
                value = values[j, i]
                if not np.isfinite(value) or x[i+1] <= xmin or x[i] >= xmax:
                    continue
                i0 = max(0, np.searchsorted(xe, x[i], side="right")-1)
                i1 = min(pixels, np.searchsorted(xe, x[i+1], side="left"))
                view = image[j0:j1, i0:i1]
                np.fmax(view, value, out=view)
    return image


def max_xy_run():
    if MAX_XY_THRESHOLD_KMS <= 0 or not np.isfinite(MAX_XY_THRESHOLD_KMS):
        raise ValueError("閾値は有限の正数にしてください")
    if not 0 < MAX_XY_ZOOM_FRACTION <= 1 or MAX_XY_PIXELS < 2:
        raise ValueError("Zoom倍率またはピクセル数が不正です")
    metadata = MAX_XY_DATA_DIR / "Toyouchi.metadata.json"
    if not metadata.is_file():
        candidates = sorted(MAX_XY_DATA_DIR.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError("入力フォルダーにmetadata JSONを置いてください")
        metadata = candidates[0]
    units = json.loads(metadata.read_text(encoding="utf-8"))["units"]
    length_au = float(units["length_cgs"]) / 1.495978707e13
    velocity_kms = float(units["velocity_cgs"]) / 1.e5
    time_kyr = float(units["time_cgs"]) / (365.25*86400*1000)
    groups = defaultdict(list)
    seen = set()
    for path in sorted(MAX_XY_DATA_DIR.glob(MAX_XY_VTK_GLOB)):
        m = re.fullmatch(r"(.+)\.block(\d+)\.out(\d+)\.(\d+)\.vtk", path.name)
        if m is None:
            continue
        prefix, block, stream, step = m.groups()
        key = (prefix, int(stream), int(step))
        identity = key + (int(block),)
        if identity in seen:
            raise ValueError(f"同じ出力・ブロックが重複: {path}")
        seen.add(identity)
        groups[key].append(path)
    if not groups:
        raise FileNotFoundError("対象VTKがありません")
    for name in ("full", "Zoom"):
        (MAX_XY_OUTPUT_DIR / name).mkdir(parents=True, exist_ok=True)
    summary = []
    vertical_axis = 1
    for (prefix, stream, step), paths in sorted(groups.items()):
        tiles, times = [], []
        bounds = [np.inf, -np.inf, np.inf, -np.inf]
        peak = None
        invalid = high_count = 0
        for path in paths:
            with path.open("rb") as f:
                f.readline()
                header = f.readline().decode("ascii", errors="replace")
            match = re.search(r"time\s*=\s*([-+\d.eE]+)", header)
            if match is None:
                raise ValueError(f"時刻が読めません: {path}")
            times.append(float(match.group(1)))
            g = pv.read(path)
            if not isinstance(g, pv.RectilinearGrid):
                raise TypeError(f"RectilinearGrid以外は未対応: {path}")
            edges = [np.asarray(a, dtype=float)*length_au for a in (g.x, g.y, g.z)]
            if any(len(a)<2 or not np.all(np.isfinite(a)) or np.any(np.diff(a)<=0)
                   for a in edges):
                raise ValueError(f"格子座標が不正: {path}")
            nx, ny, nz = [len(a)-1 for a in edges]
            velocity = np.asarray(g.cell_data["vel"], dtype=float)
            if velocity.shape != (nx*ny*nz, 3):
                raise ValueError(f"velの形状が不正: {path}")
            speed = np.linalg.norm(velocity, axis=1)*velocity_kms
            good = np.isfinite(speed)
            invalid += int((~good).sum())
            high_count += int(np.count_nonzero(good & (speed >= MAX_XY_THRESHOLD_KMS)))
            safe = np.where(good, speed, -np.inf)
            if good.any():
                n = int(np.argmax(safe))
                if peak is None or safe[n] > peak[0]:
                    k, j, i = np.unravel_index(n, (nz, ny, nx))
                    position = np.array([(edges[0][i]+edges[0][i+1])/2,
                                         (edges[1][j]+edges[1][j+1])/2,
                                         (edges[2][k]+edges[2][k+1])/2])
                    peak = (float(safe[n]), position, path.name, i, j, k)
            projected = safe.reshape(nz, ny, nx).max(axis=0)
            projected[~np.isfinite(projected)] = np.nan
            x, y = edges[0], edges[vertical_axis]
            tiles.append((x, y, projected))
            bounds = [min(bounds[0], x[0]), max(bounds[1], x[-1]),
                      min(bounds[2], y[0]), max(bounds[3], y[-1])]
            del g
        if not np.allclose(times, times[0], rtol=1.e-7, atol=1.e-12):
            raise ValueError(f"出力時刻が不一致: step={step}")
        if peak is None:
            raise ValueError(f"有限の速度がありません: step={step}")
        xhalf = (bounds[1]-bounds[0])*MAX_XY_ZOOM_FRACTION/2
        yhalf = (bounds[3]-bounds[2])*MAX_XY_ZOOM_FRACTION/2
        for name, extent in (("full", bounds), ("Zoom", [-xhalf,xhalf,-yhalf,yhalf])):
            values = max_xy_raster_max(tiles, extent, MAX_XY_PIXELS)
            positive = values[np.isfinite(values) & (values>0)]
            vmin = MAX_XY_COLOR_MIN_KMS
            if vmin is None:
                vmin = float(np.percentile(positive, 1)) if positive.size else 1.e-6
            vmax = MAX_XY_COLOR_MAX_KMS
            if vmax is None:
                vmax = max(peak[0], vmin*1.01)
            if not (np.isfinite(vmin) and np.isfinite(vmax) and 0 < vmin < vmax):
                raise ValueError("対数色範囲は0 < min < maxの有限値にしてください")
            fig, ax = plt.subplots(figsize=(9, 8))
            try:
                # 真のゼロは最低色、欠損は白。線形補間をせずセルのピークを保持。
                background = np.ma.masked_invalid(np.maximum(values, vmin))
                cmap = plt.get_cmap("viridis").copy()
                cmap.set_bad("white")
                im = ax.imshow(background, origin="lower", extent=extent,
                               interpolation="nearest", cmap=cmap, norm=LogNorm(vmin,vmax))
                overlay = np.ma.masked_where(~np.isfinite(values) | (values<MAX_XY_THRESHOLD_KMS),
                                            np.ones_like(values))
                ax.imshow(overlay, origin="lower", extent=extent, interpolation="nearest",
                          cmap=ListedColormap(["red"]), vmin=0, vmax=1, alpha=0.65)
                ax.plot([], [], color="red", linewidth=6, alpha=0.65,
                        label=f"Projected speed >= {MAX_XY_THRESHOLD_KMS:g} km/s")
                px, py = peak[1][0], peak[1][vertical_axis]
                if extent[0] <= px <= extent[1] and extent[2] <= py <= extent[3]:
                    ax.scatter([px], [py], marker="*", s=170, facecolors="cyan",
                               edgecolors="black", linewidths=0.9,
                               label=f"Global max: {peak[0]:.5g} km/s", zorder=5)
                else:
                    ax.text(0.02,0.02,"Global maximum is outside this view",
                            transform=ax.transAxes, bbox=dict(facecolor="white",alpha=0.8))
                ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="x [au]", ylabel="y [au]")
                ax.set_aspect("equal")
                ax.set_title(f"xy max-|v| projection along z ({name})\n"
                             f"t={times[0]*time_kyr:.6g} kyr; global max={peak[0]:.6g} km/s\n"
                             f"peak (x,y,z)={tuple(np.round(peak[1],3))} au")
                ax.legend(loc="upper right", fontsize=8)
                fig.colorbar(im,ax=ax,label="Maximum speed along sightline [km/s]",extend="both")
                fig.tight_layout()
                fig.savefig(MAX_XY_OUTPUT_DIR/name/f"{prefix}.out{stream}.{step:05d}.png",dpi=MAX_XY_DPI)
                if MAX_XY_SHOW:
                    plt.show()
            finally:
                plt.close(fig)
        summary.append(dict(prefix=prefix,stream=stream,step=step,time_kyr=times[0]*time_kyr,
                            max_speed_kms=peak[0],x_au=peak[1][0],y_au=peak[1][1],z_au=peak[1][2],
                            block_file=peak[2],i=peak[3],j=peak[4],k=peak[5],
                            threshold_kms=MAX_XY_THRESHOLD_KMS,cells_above_threshold=high_count,
                            nonfinite_cells=invalid))
        print(f"[xy step={step}] max={peak[0]:.6g} km/s, "
              f">=閾値: {high_count} cells, 非有限値: {invalid}")
    with (MAX_XY_OUTPUT_DIR/"maximum_speed_cells.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)
    print(f"[SAVED] {MAX_XY_OUTPUT_DIR}")


max_xy_run()


# %% cell5: xz最大速度投影（full / Zoom、独立実行）
# 中央断面ではなく視線方向の最大|v|。3次元速度の大きさを使用。
# ピクセルと重なるセルの最大値を採用：微小な高速セルも平均で消さない。
# そのためピクセル未満の構造は画面上で広がる場合があります。
from pathlib import Path
from collections import defaultdict
import csv
import json
import re
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap

MAX_XZ_DATA_DIR = Path("~/athena-project/results/○○").expanduser()
MAX_XZ_VTK_GLOB = "**/*.block*.out2.*.vtk"
MAX_XZ_OUTPUT_DIR = MAX_XZ_DATA_DIR / "checks" / "speed_max_projection" / "xz"
MAX_XZ_THRESHOLD_KMS = 100.0  # この値以上を赤で強調（手動変更可）
MAX_XZ_ZOOM_FRACTION = 0.1    # 原点中心、全域幅の1/10
MAX_XZ_PIXELS = 800          # full / Zoomそれぞれの一辺の画素数
MAX_XZ_DPI = 180
MAX_XZ_SHOW = False
MAX_XZ_COLOR_MIN_KMS = None  # None: 表示画像の正の速度の1パーセンタイル
MAX_XZ_COLOR_MAX_KMS = None  # None: 各時刻の全領域最大速度（外れ値も含む）


def max_xz_raster_max(tiles, bounds, pixels):
    """セル矩形の値を重なる全ピクセルへmaxで集約。欠損はNaN。"""
    xmin, xmax, ymin, ymax = bounds
    xe = np.linspace(xmin, xmax, pixels+1)
    ye = np.linspace(ymin, ymax, pixels+1)
    image = np.full((pixels, pixels), np.nan)
    for x, y, values in tiles:
        for j in range(len(y)-1):
            if y[j+1] <= ymin or y[j] >= ymax:
                continue
            j0 = max(0, np.searchsorted(ye, y[j], side="right")-1)
            j1 = min(pixels, np.searchsorted(ye, y[j+1], side="left"))
            for i in range(len(x)-1):
                value = values[j, i]
                if not np.isfinite(value) or x[i+1] <= xmin or x[i] >= xmax:
                    continue
                i0 = max(0, np.searchsorted(xe, x[i], side="right")-1)
                i1 = min(pixels, np.searchsorted(xe, x[i+1], side="left"))
                view = image[j0:j1, i0:i1]
                np.fmax(view, value, out=view)
    return image


def max_xz_run():
    if MAX_XZ_THRESHOLD_KMS <= 0 or not np.isfinite(MAX_XZ_THRESHOLD_KMS):
        raise ValueError("閾値は有限の正数にしてください")
    if not 0 < MAX_XZ_ZOOM_FRACTION <= 1 or MAX_XZ_PIXELS < 2:
        raise ValueError("Zoom倍率またはピクセル数が不正です")
    metadata = MAX_XZ_DATA_DIR / "Toyouchi.metadata.json"
    if not metadata.is_file():
        candidates = sorted(MAX_XZ_DATA_DIR.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError("入力フォルダーにmetadata JSONを置いてください")
        metadata = candidates[0]
    units = json.loads(metadata.read_text(encoding="utf-8"))["units"]
    length_au = float(units["length_cgs"]) / 1.495978707e13
    velocity_kms = float(units["velocity_cgs"]) / 1.e5
    time_kyr = float(units["time_cgs"]) / (365.25*86400*1000)
    groups = defaultdict(list)
    seen = set()
    for path in sorted(MAX_XZ_DATA_DIR.glob(MAX_XZ_VTK_GLOB)):
        m = re.fullmatch(r"(.+)\.block(\d+)\.out(\d+)\.(\d+)\.vtk", path.name)
        if m is None:
            continue
        prefix, block, stream, step = m.groups()
        key = (prefix, int(stream), int(step))
        identity = key + (int(block),)
        if identity in seen:
            raise ValueError(f"同じ出力・ブロックが重複: {path}")
        seen.add(identity)
        groups[key].append(path)
    if not groups:
        raise FileNotFoundError("対象VTKがありません")
    for name in ("full", "Zoom"):
        (MAX_XZ_OUTPUT_DIR / name).mkdir(parents=True, exist_ok=True)
    summary = []
    vertical_axis = 2
    for (prefix, stream, step), paths in sorted(groups.items()):
        tiles, times = [], []
        bounds = [np.inf, -np.inf, np.inf, -np.inf]
        peak = None
        invalid = high_count = 0
        for path in paths:
            with path.open("rb") as f:
                f.readline()
                header = f.readline().decode("ascii", errors="replace")
            match = re.search(r"time\s*=\s*([-+\d.eE]+)", header)
            if match is None:
                raise ValueError(f"時刻が読めません: {path}")
            times.append(float(match.group(1)))
            g = pv.read(path)
            if not isinstance(g, pv.RectilinearGrid):
                raise TypeError(f"RectilinearGrid以外は未対応: {path}")
            edges = [np.asarray(a, dtype=float)*length_au for a in (g.x, g.y, g.z)]
            if any(len(a)<2 or not np.all(np.isfinite(a)) or np.any(np.diff(a)<=0)
                   for a in edges):
                raise ValueError(f"格子座標が不正: {path}")
            nx, ny, nz = [len(a)-1 for a in edges]
            velocity = np.asarray(g.cell_data["vel"], dtype=float)
            if velocity.shape != (nx*ny*nz, 3):
                raise ValueError(f"velの形状が不正: {path}")
            speed = np.linalg.norm(velocity, axis=1)*velocity_kms
            good = np.isfinite(speed)
            invalid += int((~good).sum())
            high_count += int(np.count_nonzero(good & (speed >= MAX_XZ_THRESHOLD_KMS)))
            safe = np.where(good, speed, -np.inf)
            if good.any():
                n = int(np.argmax(safe))
                if peak is None or safe[n] > peak[0]:
                    k, j, i = np.unravel_index(n, (nz, ny, nx))
                    position = np.array([(edges[0][i]+edges[0][i+1])/2,
                                         (edges[1][j]+edges[1][j+1])/2,
                                         (edges[2][k]+edges[2][k+1])/2])
                    peak = (float(safe[n]), position, path.name, i, j, k)
            projected = safe.reshape(nz, ny, nx).max(axis=1)
            projected[~np.isfinite(projected)] = np.nan
            x, y = edges[0], edges[vertical_axis]
            tiles.append((x, y, projected))
            bounds = [min(bounds[0], x[0]), max(bounds[1], x[-1]),
                      min(bounds[2], y[0]), max(bounds[3], y[-1])]
            del g
        if not np.allclose(times, times[0], rtol=1.e-7, atol=1.e-12):
            raise ValueError(f"出力時刻が不一致: step={step}")
        if peak is None:
            raise ValueError(f"有限の速度がありません: step={step}")
        xhalf = (bounds[1]-bounds[0])*MAX_XZ_ZOOM_FRACTION/2
        yhalf = (bounds[3]-bounds[2])*MAX_XZ_ZOOM_FRACTION/2
        for name, extent in (("full", bounds), ("Zoom", [-xhalf,xhalf,-yhalf,yhalf])):
            values = max_xz_raster_max(tiles, extent, MAX_XZ_PIXELS)
            positive = values[np.isfinite(values) & (values>0)]
            vmin = MAX_XZ_COLOR_MIN_KMS
            if vmin is None:
                vmin = float(np.percentile(positive, 1)) if positive.size else 1.e-6
            vmax = MAX_XZ_COLOR_MAX_KMS
            if vmax is None:
                vmax = max(peak[0], vmin*1.01)
            if not (np.isfinite(vmin) and np.isfinite(vmax) and 0 < vmin < vmax):
                raise ValueError("対数色範囲は0 < min < maxの有限値にしてください")
            fig, ax = plt.subplots(figsize=(9, 8))
            try:
                # 真のゼロは最低色、欠損は白。線形補間をせずセルのピークを保持。
                background = np.ma.masked_invalid(np.maximum(values, vmin))
                cmap = plt.get_cmap("viridis").copy()
                cmap.set_bad("white")
                im = ax.imshow(background, origin="lower", extent=extent,
                               interpolation="nearest", cmap=cmap, norm=LogNorm(vmin,vmax))
                overlay = np.ma.masked_where(~np.isfinite(values) | (values<MAX_XZ_THRESHOLD_KMS),
                                            np.ones_like(values))
                ax.imshow(overlay, origin="lower", extent=extent, interpolation="nearest",
                          cmap=ListedColormap(["red"]), vmin=0, vmax=1, alpha=0.65)
                ax.plot([], [], color="red", linewidth=6, alpha=0.65,
                        label=f"Projected speed >= {MAX_XZ_THRESHOLD_KMS:g} km/s")
                px, py = peak[1][0], peak[1][vertical_axis]
                if extent[0] <= px <= extent[1] and extent[2] <= py <= extent[3]:
                    ax.scatter([px], [py], marker="*", s=170, facecolors="cyan",
                               edgecolors="black", linewidths=0.9,
                               label=f"Global max: {peak[0]:.5g} km/s", zorder=5)
                else:
                    ax.text(0.02,0.02,"Global maximum is outside this view",
                            transform=ax.transAxes, bbox=dict(facecolor="white",alpha=0.8))
                ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="x [au]", ylabel="z [au]")
                ax.set_aspect("equal")
                ax.set_title(f"xz max-|v| projection along y ({name})\n"
                             f"t={times[0]*time_kyr:.6g} kyr; global max={peak[0]:.6g} km/s\n"
                             f"peak (x,y,z)={tuple(np.round(peak[1],3))} au")
                ax.legend(loc="upper right", fontsize=8)
                fig.colorbar(im,ax=ax,label="Maximum speed along sightline [km/s]",extend="both")
                fig.tight_layout()
                fig.savefig(MAX_XZ_OUTPUT_DIR/name/f"{prefix}.out{stream}.{step:05d}.png",dpi=MAX_XZ_DPI)
                if MAX_XZ_SHOW:
                    plt.show()
            finally:
                plt.close(fig)
        summary.append(dict(prefix=prefix,stream=stream,step=step,time_kyr=times[0]*time_kyr,
                            max_speed_kms=peak[0],x_au=peak[1][0],y_au=peak[1][1],z_au=peak[1][2],
                            block_file=peak[2],i=peak[3],j=peak[4],k=peak[5],
                            threshold_kms=MAX_XZ_THRESHOLD_KMS,cells_above_threshold=high_count,
                            nonfinite_cells=invalid))
        print(f"[xz step={step}] max={peak[0]:.6g} km/s, "
              f">=閾値: {high_count} cells, 非有限値: {invalid}")
    with (MAX_XZ_OUTPUT_DIR/"maximum_speed_cells.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)
    print(f"[SAVED] {MAX_XZ_OUTPUT_DIR}")


max_xz_run()

# %% cell6: 全領域の最大速度・質量重み付き平均速度の時間変化（独立実行）
# mean = sum(rho * cell_volume * |v|) / sum(rho * cell_volume)
# |平均速度ベクトル|ではなく、速度の大きさの平均。
# SMRのセル体積を含める。欠落ブロック・重複・不正値は黙って除外しない。
from pathlib import Path
from collections import defaultdict
import csv
import json
import re
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

HISTORY_SPEED_DATA_DIR = Path("~/athena-project/results/○○").expanduser()
HISTORY_SPEED_VTK_GLOB = "**/*.block*.out2.*.vtk"
HISTORY_SPEED_OUTPUT_DIR = HISTORY_SPEED_DATA_DIR / "checks" / "speed_history"
HISTORY_SPEED_EXPECTED_BLOCKS = None  # SMR: Noneなら最初の時刻のブロック集合と照合
HISTORY_SPEED_LOG_Y = False          # 桁差が大きい場合はTrue（対数軸）
HISTORY_SPEED_SHOW = True
HISTORY_SPEED_DPI = 180


def history_speed_run():
    meta = HISTORY_SPEED_DATA_DIR / "Toyouchi.metadata.json"
    if not meta.is_file():
        candidates = sorted(HISTORY_SPEED_DATA_DIR.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError("入力フォルダーに単位系のmetadata JSONを置いてください")
        meta = candidates[0]
    units = json.loads(meta.read_text(encoding="utf-8"))["units"]
    velocity_kms = float(units["velocity_cgs"]) / 1.e5
    time_kyr = float(units["time_cgs"]) / (365.25 * 86400 * 1000)
    # 質量の共通換算係数は加重平均の分子・分母で相殺される。
    groups = defaultdict(dict)
    for path in sorted(HISTORY_SPEED_DATA_DIR.glob(HISTORY_SPEED_VTK_GLOB)):
        match = re.fullmatch(r"(.+)\.block(\d+)\.out(\d+)\.(\d+)\.vtk", path.name)
        if match is None:
            continue
        prefix, block, stream, step = match.groups()
        key = (prefix, int(stream), int(step))
        block = int(block)
        if block in groups[key]:
            raise ValueError(f"同じ時刻・ブロックが重複しています: {path}")
        groups[key][block] = path
    if not groups:
        raise FileNotFoundError("対象VTKがありません")
    streams = {(key[0], key[1]) for key in groups}
    if len(streams) != 1:
        raise ValueError("複数のproblem_id/outストリームがあります。VTK_GLOBで1つに絞ってください")
    reference_blocks = set(groups[min(groups)])
    if HISTORY_SPEED_EXPECTED_BLOCKS is not None and len(reference_blocks) != HISTORY_SPEED_EXPECTED_BLOCKS:
        raise ValueError("最初の時刻のブロック数がEXPECTED_BLOCKSと一致しません")
    print(f"[INFO] SMRブロック数: {len(reference_blocks)}（全時刻で同じ集合を確認）")
    rows = []
    for number, ((prefix, stream, step), blocks) in enumerate(sorted(groups.items()), 1):
        if set(blocks) != reference_blocks:
            raise ValueError(f"ブロック集合が不一致: step={step}。欠落またはAMRを確認してください")
        total_mass = weighted_speed = 0.0
        peak = -np.inf
        peak_position = None
        peak_file = None
        cell_count = 0
        times = []
        for block, path in sorted(blocks.items()):
            with path.open("rb") as f:
                f.readline()
                header = f.readline().decode("ascii", errors="replace")
            match = re.search(r"time\s*=\s*([-+\d.eE]+)", header)
            if match is None:
                raise ValueError(f"時刻を取得できません: {path}")
            times.append(float(match.group(1)))
            grid = pv.read(path)
            if not isinstance(grid, pv.RectilinearGrid):
                raise TypeError(f"RectilinearGrid以外は未対応: {path}")
            edges = [np.asarray(a, dtype=float) for a in (grid.x, grid.y, grid.z)]
            widths = [np.diff(a) for a in edges]
            if any(len(a) == 0 or not np.all(np.isfinite(a)) or np.any(a <= 0) for a in widths):
                raise ValueError(f"セル幅が不正: {path}")
            dx, dy, dz = widths
            volume = (dz[:, None, None] * dy[None, :, None] * dx[None, None, :]).ravel()
            rho = np.asarray(grid.cell_data["rho"], dtype=float)
            vel = np.asarray(grid.cell_data["vel"], dtype=float)
            if rho.shape != volume.shape or vel.shape != (volume.size, 3):
                raise ValueError(f"rho/velの配列形状が不正: {path}")
            if not np.all(np.isfinite(rho)) or np.any(rho <= 0) or not np.all(np.isfinite(vel)):
                raise ValueError(f"密度<=0または非有限値があります: {path}")
            speed = np.linalg.norm(vel, axis=1) * velocity_kms
            mass = rho * volume
            if not np.all(np.isfinite(speed)) or not np.all(np.isfinite(mass)):
                raise ValueError(f"速度・セル質量の計算が非有限値になりました: {path}")
            total_mass += float(np.sum(mass))
            weighted_speed += float(np.dot(mass, speed))
            cell_count += volume.size
            n = int(np.argmax(speed))
            if speed[n] > peak:
                peak = float(speed[n])
                k, j, i = np.unravel_index(n, (len(dz), len(dy), len(dx)))
                peak_position = [(edges[0][i]+edges[0][i+1])/2,
                                 (edges[1][j]+edges[1][j+1])/2,
                                 (edges[2][k]+edges[2][k+1])/2]
                peak_file = path.name
            del grid
        if not np.allclose(times, times[0], rtol=1.e-7, atol=1.e-12):
            raise ValueError(f"ブロック間の時刻が不一致: step={step}")
        if total_mass <= 0 or not np.isfinite(total_mass) or not np.isfinite(weighted_speed):
            raise ValueError(f"総質量または加重和が不正: step={step}")
        position_au = np.asarray(peak_position)*float(units["length_cgs"])/1.495978707e13
        rows.append(dict(step=step, time_kyr=times[0]*time_kyr,
                         maximum_speed_kms=peak, mass_weighted_mean_speed_kms=weighted_speed/total_mass,
                         gas_mass_code=total_mass, blocks=len(blocks), cells=cell_count,
                         maximum_block_file=peak_file, maximum_x_au=position_au[0],
                         maximum_y_au=position_au[1], maximum_z_au=position_au[2]))
        print(f"[{number}/{len(groups)}] t={rows[-1]['time_kyr']:.6g} kyr, "
              f"max={peak:.6g}, mass-weighted mean={weighted_speed/total_mass:.6g} km/s")
    rows.sort(key=lambda row: row["time_kyr"])
    HISTORY_SPEED_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (HISTORY_SPEED_OUTPUT_DIR / "speed_history.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(10, 6))
    try:
        t = np.array([row["time_kyr"] for row in rows])
        for key, label, color in (
            ("maximum_speed_kms", "Maximum |v| (all cells)", "tab:red"),
            ("mass_weighted_mean_speed_kms", "Mass-weighted mean |v|", "tab:blue"),
        ):
            values = np.array([row[key] for row in rows])
            if HISTORY_SPEED_LOG_Y:
                if np.any(values <= 0):
                    print(f"[WARNING] 対数軸では0の点を表示しません: {label}")
                values = np.where(values > 0, values, np.nan)
            ax.plot(t, values, marker="o", markersize=3, linewidth=1.2, color=color, label=label)
        if HISTORY_SPEED_LOG_Y:
            ax.set_yscale("log")
        else:
            ax.set_ylim(bottom=0)
        ax.set(xlabel="Time [kyr]", ylabel="Speed |v| [km/s]",
               title="Whole-domain speed history (VTK snapshots)")
        ax.grid(alpha=0.3); ax.legend(); fig.tight_layout()
        suffix = "log" if HISTORY_SPEED_LOG_Y else "linear"
        fig.savefig(HISTORY_SPEED_OUTPUT_DIR/f"speed_history_{suffix}.png", dpi=HISTORY_SPEED_DPI)
        if HISTORY_SPEED_SHOW:
            plt.show()
    finally:
        plt.close(fig)
    print(f"[SAVED] {HISTORY_SPEED_OUTPUT_DIR}")


history_speed_run()

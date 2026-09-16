# -*- coding: utf-8 -*-
# Converted from 密度map.ipynb

# %% cell 1
# ============================================================
# x-z 平面密度マップ（AMR グリッドなし）
#   ・計算領域を VTK の bounds から自動取得
#   ・全域、1/10 を別ディレクトリへ保存
#   ・シンク半径は設定欄で手動指定
#   ・シンク内部を白抜き
#   ・時間を kyr で表示（TIME_DECIMALS 桁に丸める）
#
# このファイル全体を Jupyter Notebook の1セルに貼り付けて実行できます。
# ============================================================

import copy
import gc
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, LogFormatterSciNotation
from scipy.interpolate import griddata


# ------------------------- 入出力 ----------------------------
vtk_dir = Path(
    os.path.expanduser("~/athena-project/results/〇〇")
).resolve()

output_root = Path("./xz_density_maps_auto").resolve()
output_dirs = {
    "full": output_root / "full_domain",
    "zoom_1over10": output_root / "zoom_1over10",
}
for directory in output_dirs.values():
    directory.mkdir(parents=True, exist_ok=True)

if not vtk_dir.is_dir():
    raise FileNotFoundError(f"VTK directory not found: {vtk_dir}")


# --------------------------- 単位 ----------------------------
M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15  # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0

LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
MASS_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS
DENSITY_UNIT = MASS_UNIT_MSUN / LENGTH_UNIT_AU**3
TIME_UNIT_KYR = T_UNIT_CGS / YEAR_CGS / 1.0e3


# -------------------------- 描画設定 --------------------------
FULL_RESOLUTION = 800
ZOOM_RESOLUTION = 800

FULL_PERCENTILES = (1.0, 99.0)
ZOOM_PERCENTILES = (1.0, 99.0)

# 元コードと同じカラーマップ
FULL_CMAP = "inferno"
ZOOM_CMAP = "turbo"
ZOOM_INTERPOLATION = "bilinear"

DPI = 200
FIGSIZE = (8.6, 7.2)
TIME_DECIMALS = 0       # 0: 整数 kyr。例: 10 kyr
SHOW_IMAGES = False     # 全画像を Notebook に表示するなら True
GLOBAL_SCALE_SAMPLE_COUNT = 10

# シンク半径を手動指定 [AU]
SINK_RADIUS_AU = 1000.0


# ============================================================
# 小関数
# ============================================================
def read_vtk_time_code(filename):
    with open(filename, "rb") as handle:
        header = handle.read(2048).decode("ascii", errors="ignore")
    match = re.search(
        r"time\s*=\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        header,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTK header time not found: {filename}")
    return float(match.group(1))


def get_step_number(path):
    match = re.search(r"(?:prim\.)?out2\.(\d+)", Path(path).name)
    return int(match.group(1)) if match else None


def find_density_name(grid):
    for name in ("dens", "density", "rho", "prim_dens", "prim_density"):
        if name in grid.cell_data:
            return name
    raise KeyError(
        "Density field not found. "
        f"cell_data={list(grid.cell_data.keys())}"
    )


def log_limits(values, percentiles, fallback=None):
    values = np.asarray(values).ravel()
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        if fallback is not None:
            return fallback
        raise RuntimeError("No positive finite density values for color scale.")
    vmin, vmax = np.percentile(values, percentiles)
    if not np.isfinite(vmin) or vmin <= 0.0:
        vmin = np.min(values)
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = np.max(values)
    if vmax <= vmin:
        vmin = max(vmin * 0.5, np.finfo(float).tiny)
        vmax = vmax * 2.0
    return float(vmin), float(vmax)


def extract_xz_midplane(grid, density_name):
    """y=0 と交差する block から、y=0 に最も近いセル中心層を抽出。"""
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    tolerance = 1.0e-12 * max(abs(ymin), abs(ymax), 1.0)
    if not (ymin - tolerance <= 0.0 <= ymax + tolerance):
        return None, None

    centers = grid.cell_centers().points
    if centers.size == 0:
        return None, None
    rho_code = np.asarray(grid.cell_data[density_name]).reshape(-1)

    y_values = centers[:, 1]
    unique_y = np.unique(y_values)
    y_nearest = unique_y[np.argmin(np.abs(unique_y))]
    if unique_y.size >= 2:
        dy = np.min(np.diff(np.sort(unique_y)))
        layer_tolerance = max(0.25 * abs(dy), 1.0e-12)
    else:
        layer_tolerance = 1.0e-10
    mask = np.abs(y_values - y_nearest) <= layer_tolerance

    points_au = centers[mask] * LENGTH_UNIT_AU
    rho = rho_code[mask] * DENSITY_UNIT
    return points_au, rho


def domain_from_step(files):
    """実際の VTK block bounds から x-z 計算領域を得る。"""
    xmins, xmaxs, zmins, zmaxs = [], [], [], []
    for filename in files:
        grid = pv.read(filename)
        xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
        tol = 1.0e-12 * max(abs(ymin), abs(ymax), 1.0)
        if ymin - tol <= 0.0 <= ymax + tol:
            xmins.append(xmin * LENGTH_UNIT_AU)
            xmaxs.append(xmax * LENGTH_UNIT_AU)
            zmins.append(zmin * LENGTH_UNIT_AU)
            zmaxs.append(zmax * LENGTH_UNIT_AU)
        del grid
    if not xmins:
        raise RuntimeError("No blocks intersect the y=0 plane.")
    return min(xmins), max(xmaxs), min(zmins), max(zmaxs)


def fractional_bounds(full_bounds, fraction):
    """原点を中心に、全域の各辺を fraction 倍した表示範囲を返す。"""
    xmin, xmax, zmin, zmax = full_bounds
    x_half = 0.5 * (xmax - xmin) * fraction
    z_half = 0.5 * (zmax - zmin) * fraction
    return -x_half, x_half, -z_half, z_half


def load_step_slice(files, density_name):
    point_arrays, density_arrays = [], []
    for filename in files:
        try:
            grid = pv.read(filename)
            points, density = extract_xz_midplane(grid, density_name)
            del grid
            if points is not None and points.size:
                point_arrays.append(points)
                density_arrays.append(density)
        except Exception as error:
            print(f"[WARNING] {filename}: {error}")
    if not point_arrays:
        raise RuntimeError("No x-z midplane cells were extracted.")
    points = np.vstack(point_arrays)
    density = np.concatenate(density_arrays)
    valid = (
        np.all(np.isfinite(points), axis=1)
        & np.isfinite(density)
        & (density > 0.0)
    )
    return points[valid], density[valid]


def interpolate_density(points, density, bounds, resolution):
    xmin, xmax, zmin, zmax = bounds
    x_axis = np.linspace(xmin, xmax, resolution)
    z_axis = np.linspace(zmin, zmax, resolution)
    X, Z = np.meshgrid(x_axis, z_axis)

    # 表示域近傍だけを補間に使用し、ズーム時のメモリを節約する。
    xpad = 0.05 * (xmax - xmin)
    zpad = 0.05 * (zmax - zmin)
    use = (
        (points[:, 0] >= xmin - xpad)
        & (points[:, 0] <= xmax + xpad)
        & (points[:, 2] >= zmin - zpad)
        & (points[:, 2] <= zmax + zpad)
    )
    source_points = points[use][:, (0, 2)]
    source_density = density[use]
    if source_points.shape[0] < 4:
        source_points = points[:, (0, 2)]
        source_density = density

    linear = griddata(source_points, source_density, (X, Z), method="linear")
    if np.any(~np.isfinite(linear)):
        nearest = griddata(source_points, source_density, (X, Z), method="nearest")
        density_map = np.where(np.isfinite(linear), linear, nearest)
        del nearest
    else:
        density_map = linear
    density_map[(~np.isfinite(density_map)) | (density_map <= 0.0)] = np.nan
    return X, Z, density_map, int(source_points.shape[0])


def scientific_axis_scale(bounds):
    largest = max(abs(value) for value in bounds)
    if not np.isfinite(largest) or largest <= 0.0:
        return 1.0, 0
    exponent = int(np.floor(np.log10(largest)))
    scale = 10.0**exponent
    return scale, exponent


def plot_and_save(
    X, Z, density_map, bounds, norm, cmap_name, interpolation,
    step, time_kyr, sink_radius_au, source_cells, label, output_path,
):
    # シンクは補間後に幾何学的にマスクするため、常に同じ真円になる。
    if sink_radius_au > 0.0:
        sink_mask = np.hypot(X, Z) <= sink_radius_au
        density_plot = np.ma.array(density_map, mask=sink_mask)
    else:
        density_plot = np.ma.masked_invalid(density_map)

    cmap = copy.copy(plt.get_cmap(cmap_name))
    cmap.set_bad("white")

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if label == "full":
        image = ax.pcolormesh(
            X, Z, density_plot, cmap=cmap, norm=norm,
            shading="auto", rasterized=True,
        )
    else:
        image = ax.imshow(
            density_plot,
            origin="lower",
            extent=[bounds[0], bounds[1], bounds[2], bounds[3]],
            cmap=cmap,
            norm=norm,
            interpolation=interpolation,
            aspect="equal",
        )

    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect("equal")
    ax.axhline(0.0, color="gray", lw=0.6, alpha=0.55)
    ax.axvline(0.0, color="gray", lw=0.6, alpha=0.55)

    scale, exponent = scientific_axis_scale(bounds)
    formatter = FuncFormatter(lambda value, position: f"{value / scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\,{{\rm AU}}]$")
    ax.set_ylabel(rf"$z\ [10^{{{exponent}}}\,{{\rm AU}}]$")

    time_text = f"{time_kyr:.{TIME_DECIMALS}f}"
    title_suffix = {
        "full": "full domain",
        "zoom_1over10": "1/10 zoom",
    }[label]
    ax.set_title(
        "Density map on the x-z midplane\n"
        f"step={step:05d}, t={time_text} kyr, {title_suffix}",
        fontsize=14,
    )

    colorbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    colorbar.ax.yaxis.set_major_formatter(LogFormatterSciNotation())
    colorbar.set_label(r"$\rho\ [M_\odot\, {\rm AU}^{-3}]$")

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    if SHOW_IMAGES:
        plt.show()
    plt.close(fig)
    del fig, density_plot, image


# ============================================================
# VTK ファイルの整理
# ============================================================
files_by_step = defaultdict(list)
for path in vtk_dir.glob("*.vtk"):
    if "Toyouchi.block" not in path.name:
        continue
    step = get_step_number(path)
    if step is not None:
        files_by_step[step].append(str(path))

steps = sorted(files_by_step)
for step in steps:
    files_by_step[step].sort()
if not steps:
    raise FileNotFoundError(f"No Athena++ block VTK files found in {vtk_dir}")

test_grid = pv.read(files_by_step[steps[0]][0])
density_name = find_density_name(test_grid)
del test_grid

full_bounds = domain_from_step(files_by_step[steps[0]])
plot_bounds = {
    "full": full_bounds,
    "zoom_1over10": fractional_bounds(full_bounds, 1.0 / 10.0),
}

sink_radius_au = float(SINK_RADIUS_AU)
if not np.isfinite(sink_radius_au) or sink_radius_au < 0.0:
    raise ValueError("SINK_RADIUS_AU must be a finite value >= 0 [AU].")

print(f"[INFO] Input          : {vtk_dir}")
print(f"[INFO] Steps          : {len(steps)}")
print(f"[INFO] Density field  : {density_name}")
print(f"[INFO] Full bounds AU : {full_bounds}")
for key, directory in output_dirs.items():
    print(f"[INFO] Output {key:13s}: {directory}")


# ============================================================
# 全体図用の共通カラースケールをサンプル時刻から決定
# シンク内部は統計から除外する。
# ============================================================
sample_indices = np.unique(
    np.linspace(0, len(steps) - 1, min(GLOBAL_SCALE_SAMPLE_COUNT, len(steps))).astype(int)
)
sample_density = []
for sample_index in sample_indices:
    step = steps[sample_index]
    points, density = load_step_slice(files_by_step[step], density_name)
    outside_sink = np.hypot(points[:, 0], points[:, 2]) > sink_radius_au
    if np.any(outside_sink):
        sample_density.append(density[outside_sink])
    del points, density, outside_sink
    gc.collect()

full_vmin, full_vmax = log_limits(
    np.concatenate(sample_density), FULL_PERCENTILES
)
full_norm = LogNorm(vmin=full_vmin, vmax=full_vmax)
del sample_density
gc.collect()

print(f"[INFO] Full color range: {full_vmin:.6e} -- {full_vmax:.6e}")


# ============================================================
# 各 step を1つずつ読み込み、描画後すぐ解放
# ============================================================
for index, step in enumerate(steps, start=1):
    files = files_by_step[step]
    points, density = load_step_slice(files, density_name)
    time_code = read_vtk_time_code(files[0])
    time_kyr = time_code * TIME_UNIT_KYR

    radius_xz = np.hypot(points[:, 0], points[:, 2])
    outside_sink = radius_xz > sink_radius_au

    print(
        f"[{index:4d}/{len(steps):4d}] step={step:05d}, "
        f"t={time_kyr:.{TIME_DECIMALS}f} kyr, cells={len(density)}"
    )

    for label in ("full", "zoom_1over10"):
        bounds = plot_bounds[label]
        resolution = FULL_RESOLUTION if label == "full" else ZOOM_RESOLUTION
        X, Z, density_map, source_cells = interpolate_density(
            points, density, bounds, resolution
        )

        if label == "full":
            norm = full_norm
            cmap_name = FULL_CMAP
            interpolation = "none"
        else:
            xmin, xmax, zmin, zmax = bounds
            in_view = (
                (points[:, 0] >= xmin) & (points[:, 0] <= xmax)
                & (points[:, 2] >= zmin) & (points[:, 2] <= zmax)
                & outside_sink
            )
            zoom_vmin, zoom_vmax = log_limits(
                density[in_view],
                ZOOM_PERCENTILES,
                fallback=(full_vmin, full_vmax),
            )
            norm = LogNorm(vmin=zoom_vmin, vmax=zoom_vmax)
            cmap_name = ZOOM_CMAP
            interpolation = ZOOM_INTERPOLATION

        output_path = output_dirs[label] / (
            f"density_xz_step_{step:05d}_{label}.png"
        )
        plot_and_save(
            X=X,
            Z=Z,
            density_map=density_map,
            bounds=bounds,
            norm=norm,
            cmap_name=cmap_name,
            interpolation=interpolation,
            step=step,
            time_kyr=time_kyr,
            sink_radius_au=sink_radius_au,
            source_cells=source_cells,
            label=label,
            output_path=output_path,
        )
        print(f"    saved: {output_path}")
        del X, Z, density_map, norm
        gc.collect()

    del points, density, radius_xz, outside_sink
    gc.collect()

print("[DONE] All density maps were saved.")

# %% cell 2
# ============================================================
# x-y 平面密度マップ（AMR グリッドなし）
#   ・計算領域を VTK の bounds から自動取得
#   ・全域、1/10 を別ディレクトリへ保存
#   ・シンク半径は設定欄で手動指定
#   ・シンク内部を白抜き
#   ・時間を kyr で表示（TIME_DECIMALS 桁に丸める）
#
# このファイル全体を Jupyter Notebook の1セルに貼り付けて実行できます。
# ============================================================

import copy
import gc
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, LogFormatterSciNotation
from scipy.interpolate import griddata


# ------------------------- 入出力 ----------------------------
vtk_dir = Path(
    os.path.expanduser("~/athena-project/results/〇〇")
).resolve()

output_root = Path("./xy_density_maps_auto").resolve()
output_dirs = {
    "full": output_root / "full_domain",
    "zoom_1over10": output_root / "zoom_1over10",
}
for directory in output_dirs.values():
    directory.mkdir(parents=True, exist_ok=True)

if not vtk_dir.is_dir():
    raise FileNotFoundError(f"VTK directory not found: {vtk_dir}")


# --------------------------- 単位 ----------------------------
M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15  # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0

LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
MASS_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS
DENSITY_UNIT = MASS_UNIT_MSUN / LENGTH_UNIT_AU**3
TIME_UNIT_KYR = T_UNIT_CGS / YEAR_CGS / 1.0e3


# -------------------------- 描画設定 --------------------------
FULL_RESOLUTION = 800
ZOOM_RESOLUTION = 800

FULL_PERCENTILES = (1.0, 99.0)
ZOOM_PERCENTILES = (1.0, 99.0)

# 元コードと同じカラーマップ
FULL_CMAP = "inferno"
ZOOM_CMAP = "turbo"
ZOOM_INTERPOLATION = "bilinear"

DPI = 200
FIGSIZE = (8.6, 7.2)
TIME_DECIMALS = 0       # 0: 整数 kyr。例: 10 kyr
SHOW_IMAGES = False     # 全画像を Notebook に表示するなら True
GLOBAL_SCALE_SAMPLE_COUNT = 10

# シンク半径を手動指定 [AU]
SINK_RADIUS_AU = 1000.0


# ============================================================
# 小関数
# ============================================================
def read_vtk_time_code(filename):
    with open(filename, "rb") as handle:
        header = handle.read(2048).decode("ascii", errors="ignore")
    match = re.search(
        r"time\s*=\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        header,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTK header time not found: {filename}")
    return float(match.group(1))


def get_step_number(path):
    match = re.search(r"(?:prim\.)?out2\.(\d+)", Path(path).name)
    return int(match.group(1)) if match else None


def find_density_name(grid):
    for name in ("dens", "density", "rho", "prim_dens", "prim_density"):
        if name in grid.cell_data:
            return name
    raise KeyError(
        "Density field not found. "
        f"cell_data={list(grid.cell_data.keys())}"
    )


def log_limits(values, percentiles, fallback=None):
    values = np.asarray(values).ravel()
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        if fallback is not None:
            return fallback
        raise RuntimeError("No positive finite density values for color scale.")
    vmin, vmax = np.percentile(values, percentiles)
    if not np.isfinite(vmin) or vmin <= 0.0:
        vmin = np.min(values)
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = np.max(values)
    if vmax <= vmin:
        vmin = max(vmin * 0.5, np.finfo(float).tiny)
        vmax = vmax * 2.0
    return float(vmin), float(vmax)


def extract_xy_midplane(grid, density_name):
    """z=0 と交差する block から、z=0 に最も近いセル中心層を抽出。"""
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    tolerance = 1.0e-12 * max(abs(zmin), abs(zmax), 1.0)
    if not (zmin - tolerance <= 0.0 <= zmax + tolerance):
        return None, None

    centers = grid.cell_centers().points
    if centers.size == 0:
        return None, None
    rho_code = np.asarray(grid.cell_data[density_name]).reshape(-1)

    z_values = centers[:, 2]
    unique_z = np.unique(z_values)
    z_nearest = unique_z[np.argmin(np.abs(unique_z))]
    if unique_z.size >= 2:
        dz = np.min(np.diff(np.sort(unique_z)))
        layer_tolerance = max(0.25 * abs(dz), 1.0e-12)
    else:
        layer_tolerance = 1.0e-10
    mask = np.abs(z_values - z_nearest) <= layer_tolerance

    points_au = centers[mask] * LENGTH_UNIT_AU
    rho = rho_code[mask] * DENSITY_UNIT
    return points_au, rho


def domain_from_step(files):
    """実際の VTK block bounds から x-y 計算領域を得る。"""
    xmins, xmaxs, ymins, ymaxs = [], [], [], []
    for filename in files:
        grid = pv.read(filename)
        xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
        tol = 1.0e-12 * max(abs(zmin), abs(zmax), 1.0)
        if zmin - tol <= 0.0 <= zmax + tol:
            xmins.append(xmin * LENGTH_UNIT_AU)
            xmaxs.append(xmax * LENGTH_UNIT_AU)
            ymins.append(ymin * LENGTH_UNIT_AU)
            ymaxs.append(ymax * LENGTH_UNIT_AU)
        del grid
    if not xmins:
        raise RuntimeError("No blocks intersect the z=0 plane.")
    return min(xmins), max(xmaxs), min(ymins), max(ymaxs)


def fractional_bounds(full_bounds, fraction):
    """原点を中心に、全域の各辺を fraction 倍した表示範囲を返す。"""
    xmin, xmax, ymin, ymax = full_bounds
    x_half = 0.5 * (xmax - xmin) * fraction
    y_half = 0.5 * (ymax - ymin) * fraction
    return -x_half, x_half, -y_half, y_half


def load_step_slice(files, density_name):
    point_arrays, density_arrays = [], []
    for filename in files:
        try:
            grid = pv.read(filename)
            points, density = extract_xy_midplane(grid, density_name)
            del grid
            if points is not None and points.size:
                point_arrays.append(points)
                density_arrays.append(density)
        except Exception as error:
            print(f"[WARNING] {filename}: {error}")
    if not point_arrays:
        raise RuntimeError("No x-y midplane cells were extracted.")
    points = np.vstack(point_arrays)
    density = np.concatenate(density_arrays)
    valid = (
        np.all(np.isfinite(points), axis=1)
        & np.isfinite(density)
        & (density > 0.0)
    )
    return points[valid], density[valid]


def interpolate_density(points, density, bounds, resolution):
    xmin, xmax, ymin, ymax = bounds
    x_axis = np.linspace(xmin, xmax, resolution)
    y_axis = np.linspace(ymin, ymax, resolution)
    X, Y = np.meshgrid(x_axis, y_axis)

    # 表示域近傍だけを補間に使用し、ズーム時のメモリを節約する。
    xpad = 0.05 * (xmax - xmin)
    ypad = 0.05 * (ymax - ymin)
    use = (
        (points[:, 0] >= xmin - xpad)
        & (points[:, 0] <= xmax + xpad)
        & (points[:, 1] >= ymin - ypad)
        & (points[:, 1] <= ymax + ypad)
    )
    source_points = points[use][:, (0, 1)]
    source_density = density[use]
    if source_points.shape[0] < 4:
        source_points = points[:, (0, 1)]
        source_density = density

    linear = griddata(source_points, source_density, (X, Y), method="linear")
    if np.any(~np.isfinite(linear)):
        nearest = griddata(source_points, source_density, (X, Y), method="nearest")
        density_map = np.where(np.isfinite(linear), linear, nearest)
        del nearest
    else:
        density_map = linear
    density_map[(~np.isfinite(density_map)) | (density_map <= 0.0)] = np.nan
    return X, Y, density_map, int(source_points.shape[0])


def scientific_axis_scale(bounds):
    largest = max(abs(value) for value in bounds)
    if not np.isfinite(largest) or largest <= 0.0:
        return 1.0, 0
    exponent = int(np.floor(np.log10(largest)))
    scale = 10.0**exponent
    return scale, exponent


def plot_and_save(
    X, Y, density_map, bounds, norm, cmap_name, interpolation,
    step, time_kyr, sink_radius_au, source_cells, label, output_path,
):
    # シンクは補間後に幾何学的にマスクするため、常に同じ真円になる。
    if sink_radius_au > 0.0:
        sink_mask = np.hypot(X, Y) <= sink_radius_au
        density_plot = np.ma.array(density_map, mask=sink_mask)
    else:
        density_plot = np.ma.masked_invalid(density_map)

    cmap = copy.copy(plt.get_cmap(cmap_name))
    cmap.set_bad("white")

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if label == "full":
        image = ax.pcolormesh(
            X, Y, density_plot, cmap=cmap, norm=norm,
            shading="auto", rasterized=True,
        )
    else:
        image = ax.imshow(
            density_plot,
            origin="lower",
            extent=[bounds[0], bounds[1], bounds[2], bounds[3]],
            cmap=cmap,
            norm=norm,
            interpolation=interpolation,
            aspect="equal",
        )

    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect("equal")
    ax.axhline(0.0, color="gray", lw=0.6, alpha=0.55)
    ax.axvline(0.0, color="gray", lw=0.6, alpha=0.55)

    scale, exponent = scientific_axis_scale(bounds)
    formatter = FuncFormatter(lambda value, position: f"{value / scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\,{{\rm AU}}]$")
    ax.set_ylabel(rf"$y\ [10^{{{exponent}}}\,{{\rm AU}}]$")

    time_text = f"{time_kyr:.{TIME_DECIMALS}f}"
    title_suffix = {
        "full": "full domain",
        "zoom_1over10": "1/10 zoom",
    }[label]
    ax.set_title(
        "Density map on the x-y midplane\n"
        f"step={step:05d}, t={time_text} kyr, {title_suffix}",
        fontsize=14,
    )

    colorbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    colorbar.ax.yaxis.set_major_formatter(LogFormatterSciNotation())
    colorbar.set_label(r"$\rho\ [M_\odot\, {\rm AU}^{-3}]$")

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    if SHOW_IMAGES:
        plt.show()
    plt.close(fig)
    del fig, density_plot, image


# ============================================================
# VTK ファイルの整理
# ============================================================
files_by_step = defaultdict(list)
for path in vtk_dir.glob("*.vtk"):
    if "Toyouchi.block" not in path.name:
        continue
    step = get_step_number(path)
    if step is not None:
        files_by_step[step].append(str(path))

steps = sorted(files_by_step)
for step in steps:
    files_by_step[step].sort()
if not steps:
    raise FileNotFoundError(f"No Athena++ block VTK files found in {vtk_dir}")

test_grid = pv.read(files_by_step[steps[0]][0])
density_name = find_density_name(test_grid)
del test_grid

full_bounds = domain_from_step(files_by_step[steps[0]])
plot_bounds = {
    "full": full_bounds,
    "zoom_1over10": fractional_bounds(full_bounds, 1.0 / 10.0),
}

sink_radius_au = float(SINK_RADIUS_AU)
if not np.isfinite(sink_radius_au) or sink_radius_au < 0.0:
    raise ValueError("SINK_RADIUS_AU must be a finite value >= 0 [AU].")

print(f"[INFO] Input          : {vtk_dir}")
print(f"[INFO] Steps          : {len(steps)}")
print(f"[INFO] Density field  : {density_name}")
print(f"[INFO] Full bounds AU : {full_bounds}")
for key, directory in output_dirs.items():
    print(f"[INFO] Output {key:13s}: {directory}")


# ============================================================
# 全体図用の共通カラースケールをサンプル時刻から決定
# シンク内部は統計から除外する。
# ============================================================
sample_indices = np.unique(
    np.linspace(0, len(steps) - 1, min(GLOBAL_SCALE_SAMPLE_COUNT, len(steps))).astype(int)
)
sample_density = []
for sample_index in sample_indices:
    step = steps[sample_index]
    points, density = load_step_slice(files_by_step[step], density_name)
    outside_sink = np.hypot(points[:, 0], points[:, 1]) > sink_radius_au
    if np.any(outside_sink):
        sample_density.append(density[outside_sink])
    del points, density, outside_sink
    gc.collect()

full_vmin, full_vmax = log_limits(
    np.concatenate(sample_density), FULL_PERCENTILES
)
full_norm = LogNorm(vmin=full_vmin, vmax=full_vmax)
del sample_density
gc.collect()

print(f"[INFO] Full color range: {full_vmin:.6e} -- {full_vmax:.6e}")


# ============================================================
# 各 step を1つずつ読み込み、描画後すぐ解放
# ============================================================
for index, step in enumerate(steps, start=1):
    files = files_by_step[step]
    points, density = load_step_slice(files, density_name)
    time_code = read_vtk_time_code(files[0])
    time_kyr = time_code * TIME_UNIT_KYR

    radius_xy = np.hypot(points[:, 0], points[:, 1])
    outside_sink = radius_xy > sink_radius_au

    print(
        f"[{index:4d}/{len(steps):4d}] step={step:05d}, "
        f"t={time_kyr:.{TIME_DECIMALS}f} kyr, cells={len(density)}"
    )

    for label in ("full", "zoom_1over10"):
        bounds = plot_bounds[label]
        resolution = FULL_RESOLUTION if label == "full" else ZOOM_RESOLUTION
        X, Y, density_map, source_cells = interpolate_density(
            points, density, bounds, resolution
        )

        if label == "full":
            norm = full_norm
            cmap_name = FULL_CMAP
            interpolation = "none"
        else:
            xmin, xmax, ymin, ymax = bounds
            in_view = (
                (points[:, 0] >= xmin) & (points[:, 0] <= xmax)
                & (points[:, 1] >= ymin) & (points[:, 1] <= ymax)
                & outside_sink
            )
            zoom_vmin, zoom_vmax = log_limits(
                density[in_view],
                ZOOM_PERCENTILES,
                fallback=(full_vmin, full_vmax),
            )
            norm = LogNorm(vmin=zoom_vmin, vmax=zoom_vmax)
            cmap_name = ZOOM_CMAP
            interpolation = ZOOM_INTERPOLATION

        output_path = output_dirs[label] / (
            f"density_xy_step_{step:05d}_{label}.png"
        )
        plot_and_save(
            X=X,
            Y=Y,
            density_map=density_map,
            bounds=bounds,
            norm=norm,
            cmap_name=cmap_name,
            interpolation=interpolation,
            step=step,
            time_kyr=time_kyr,
            sink_radius_au=sink_radius_au,
            source_cells=source_cells,
            label=label,
            output_path=output_path,
        )
        print(f"    saved: {output_path}")
        del X, Y, density_map, norm
        gc.collect()

    del points, density, radius_xy, outside_sink
    gc.collect()

print("[DONE] All density maps were saved.")


# %% cell 3
# ============================================================
# x-y 平面密度マップ（実セルグリッド線あり）
#   ・計算領域を VTK の bounds から自動取得
#   ・全域、1/10 を別ディレクトリへ保存
#   ・シンク半径は設定欄で手動指定
#   ・実際のセル境界を取得して黒線で重ねる
#   ・シンク内部を白抜き
#   ・時間を kyr で表示（TIME_DECIMALS 桁に丸める）
#
# このファイル全体を Jupyter Notebook の1セルに貼り付けて実行できます。
# ============================================================


import copy
import gc
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from matplotlib.ticker import FuncFormatter, LogFormatterSciNotation
from scipy.interpolate import griddata


# ------------------------- 入出力 ----------------------------
vtk_dir = Path(
    os.path.expanduser("~/athena-project/results/〇〇")
).resolve()

output_root = Path("./xy_density_maps_with_cell_grid").resolve()
output_dirs = {
    "full": output_root / "full_domain",
    "zoom_1over10": output_root / "zoom_1over10",
}
for directory in output_dirs.values():
    directory.mkdir(parents=True, exist_ok=True)

if not vtk_dir.is_dir():
    raise FileNotFoundError(f"VTK directory not found: {vtk_dir}")


# --------------------------- 単位 ----------------------------
M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15  # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0

LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
MASS_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS
DENSITY_UNIT = MASS_UNIT_MSUN / LENGTH_UNIT_AU**3
TIME_UNIT_KYR = T_UNIT_CGS / YEAR_CGS / 1.0e3


# -------------------------- 描画設定 --------------------------
FULL_RESOLUTION = 800
ZOOM_RESOLUTION = 800

FULL_PERCENTILES = (1.0, 99.0)
ZOOM_PERCENTILES = (1.0, 99.0)

# 元コードと同じカラーマップ
FULL_CMAP = "inferno"
ZOOM_CMAP = "turbo"
ZOOM_INTERPOLATION = "bilinear"

DPI = 200
FIGSIZE = (8.6, 7.2)
TIME_DECIMALS = 0       # 0: 整数 kyr。例: 10 kyr
SHOW_IMAGES = False     # 全画像を Notebook に表示するなら True
GLOBAL_SCALE_SAMPLE_COUNT = 10

# 実セルグリッド線の表示設定
FULL_GRID_LINEWIDTH = 0.18
FULL_GRID_ALPHA = 0.28
ZOOM_GRID_LINEWIDTH = 0.35
ZOOM_GRID_ALPHA = 0.45
FULL_GRID_COLOR = "white"
ZOOM_GRID_COLOR = "black"

# シンク半径を手動指定 [AU]
SINK_RADIUS_AU = 1000.0


# ============================================================
# 小関数
# ============================================================
def read_vtk_time_code(filename):
    with open(filename, "rb") as handle:
        header = handle.read(2048).decode("ascii", errors="ignore")
    match = re.search(
        r"time\s*=\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        header,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTK header time not found: {filename}")
    return float(match.group(1))


def get_step_number(path):
    match = re.search(r"(?:prim\.)?out2\.(\d+)", Path(path).name)
    return int(match.group(1)) if match else None


def find_density_name(grid):
    for name in ("dens", "density", "rho", "prim_dens", "prim_density"):
        if name in grid.cell_data:
            return name
    raise KeyError(
        "Density field not found. "
        f"cell_data={list(grid.cell_data.keys())}"
    )


def log_limits(values, percentiles, fallback=None):
    values = np.asarray(values).ravel()
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        if fallback is not None:
            return fallback
        raise RuntimeError("No positive finite density values for color scale.")
    vmin, vmax = np.percentile(values, percentiles)
    if not np.isfinite(vmin) or vmin <= 0.0:
        vmin = np.min(values)
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = np.max(values)
    if vmax <= vmin:
        vmin = max(vmin * 0.5, np.finfo(float).tiny)
        vmax = vmax * 2.0
    return float(vmin), float(vmax)


def extract_xy_midplane(grid, density_name):
    """z=0 と交差する block から、z=0 に最も近いセル中心層を抽出。"""
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    tolerance = 1.0e-12 * max(abs(zmin), abs(zmax), 1.0)
    if not (zmin - tolerance <= 0.0 <= zmax + tolerance):
        return None, None

    centers = grid.cell_centers().points
    if centers.size == 0:
        return None, None
    rho_code = np.asarray(grid.cell_data[density_name]).reshape(-1)

    z_values = centers[:, 2]
    unique_z = np.unique(z_values)
    z_nearest = unique_z[np.argmin(np.abs(unique_z))]
    if unique_z.size >= 2:
        dz = np.min(np.diff(np.sort(unique_z)))
        layer_tolerance = max(0.25 * abs(dz), 1.0e-12)
    else:
        layer_tolerance = 1.0e-10
    mask = np.abs(z_values - z_nearest) <= layer_tolerance

    points_au = centers[mask] * LENGTH_UNIT_AU
    rho = rho_code[mask] * DENSITY_UNIT
    return points_au, rho


def domain_from_step(files):
    """実際の VTK block bounds から x-y 計算領域を得る。"""
    xmins, xmaxs, ymins, ymaxs = [], [], [], []
    for filename in files:
        grid = pv.read(filename)
        xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
        tol = 1.0e-12 * max(abs(zmin), abs(zmax), 1.0)
        if zmin - tol <= 0.0 <= zmax + tol:
            xmins.append(xmin * LENGTH_UNIT_AU)
            xmaxs.append(xmax * LENGTH_UNIT_AU)
            ymins.append(ymin * LENGTH_UNIT_AU)
            ymaxs.append(ymax * LENGTH_UNIT_AU)
        del grid
    if not xmins:
        raise RuntimeError("No blocks intersect the z=0 plane.")
    return min(xmins), max(xmaxs), min(ymins), max(ymaxs)


def fractional_bounds(full_bounds, fraction):
    """原点を中心に、全域の各辺を fraction 倍した表示範囲を返す。"""
    xmin, xmax, ymin, ymax = full_bounds
    x_half = 0.5 * (xmax - xmin) * fraction
    y_half = 0.5 * (ymax - ymin) * fraction
    return -x_half, x_half, -y_half, y_half


def extract_xy_cell_grid_segments(grid):
    """z=0断面と交差するblockの実セル境界線をx-y面で返す。"""
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    tolerance = 1.0e-12 * max(abs(zmin), abs(zmax), 1.0)
    if not (zmin - tolerance <= 0.0 <= zmax + tolerance):
        return []

    # Athena++ legacy VTKは通常RectilinearGrid。座標配列が得られない
    # 場合にもpointsから境界座標を復元できるようにする。
    if hasattr(grid, "x") and hasattr(grid, "y"):
        x_edges = np.asarray(grid.x, dtype=float)
        y_edges = np.asarray(grid.y, dtype=float)
    else:
        x_edges = np.unique(np.asarray(grid.points[:, 0], dtype=float))
        y_edges = np.unique(np.asarray(grid.points[:, 1], dtype=float))

    x_edges = x_edges[np.isfinite(x_edges)] * LENGTH_UNIT_AU
    y_edges = y_edges[np.isfinite(y_edges)] * LENGTH_UNIT_AU
    xmin_au, xmax_au = xmin * LENGTH_UNIT_AU, xmax * LENGTH_UNIT_AU
    ymin_au, ymax_au = ymin * LENGTH_UNIT_AU, ymax * LENGTH_UNIT_AU

    segments = []
    for x_edge in x_edges:
        segments.append(((x_edge, ymin_au), (x_edge, ymax_au)))
    for y_edge in y_edges:
        segments.append(((xmin_au, y_edge), (xmax_au, y_edge)))
    return segments


def deduplicate_segments(segments):
    """隣接block間で重複する同一線分を除く。"""
    unique = {}
    for segment in segments:
        array = np.asarray(segment, dtype=float)
        key = tuple(np.round(array.ravel(), decimals=8))
        unique[key] = array
    if not unique:
        return np.empty((0, 2, 2), dtype=float)
    return np.stack(list(unique.values()))


def load_step_slice(files, density_name, collect_grid=False):
    point_arrays, density_arrays = [], []
    grid_segments = []
    for filename in files:
        try:
            grid = pv.read(filename)
            points, density = extract_xy_midplane(grid, density_name)
            if collect_grid:
                grid_segments.extend(extract_xy_cell_grid_segments(grid))
            del grid
            if points is not None and points.size:
                point_arrays.append(points)
                density_arrays.append(density)
        except Exception as error:
            print(f"[WARNING] {filename}: {error}")
    if not point_arrays:
        raise RuntimeError("No x-y midplane cells were extracted.")
    points = np.vstack(point_arrays)
    density = np.concatenate(density_arrays)
    valid = (
        np.all(np.isfinite(points), axis=1)
        & np.isfinite(density)
        & (density > 0.0)
    )
    segments = (
        deduplicate_segments(grid_segments)
        if collect_grid
        else np.empty((0, 2, 2), dtype=float)
    )
    return points[valid], density[valid], segments


def interpolate_density(points, density, bounds, resolution):
    xmin, xmax, ymin, ymax = bounds
    x_axis = np.linspace(xmin, xmax, resolution)
    y_axis = np.linspace(ymin, ymax, resolution)
    X, Y = np.meshgrid(x_axis, y_axis)

    # 表示域近傍だけを補間に使用し、ズーム時のメモリを節約する。
    xpad = 0.05 * (xmax - xmin)
    ypad = 0.05 * (ymax - ymin)
    use = (
        (points[:, 0] >= xmin - xpad)
        & (points[:, 0] <= xmax + xpad)
        & (points[:, 1] >= ymin - ypad)
        & (points[:, 1] <= ymax + ypad)
    )
    source_points = points[use][:, (0, 1)]
    source_density = density[use]
    if source_points.shape[0] < 4:
        source_points = points[:, (0, 1)]
        source_density = density

    linear = griddata(source_points, source_density, (X, Y), method="linear")
    if np.any(~np.isfinite(linear)):
        nearest = griddata(source_points, source_density, (X, Y), method="nearest")
        density_map = np.where(np.isfinite(linear), linear, nearest)
        del nearest
    else:
        density_map = linear
    density_map[(~np.isfinite(density_map)) | (density_map <= 0.0)] = np.nan
    return X, Y, density_map, int(source_points.shape[0])


def scientific_axis_scale(bounds):
    largest = max(abs(value) for value in bounds)
    if not np.isfinite(largest) or largest <= 0.0:
        return 1.0, 0
    exponent = int(np.floor(np.log10(largest)))
    scale = 10.0**exponent
    return scale, exponent


def plot_and_save(
    X, Y, density_map, bounds, norm, cmap_name, interpolation,
    step, time_kyr, sink_radius_au, source_cells, grid_segments,
    label, output_path,
):
    # シンクは補間後に幾何学的にマスクするため、常に同じ真円になる。
    if sink_radius_au > 0.0:
        sink_mask = np.hypot(X, Y) <= sink_radius_au
        density_plot = np.ma.array(density_map, mask=sink_mask)
    else:
        density_plot = np.ma.masked_invalid(density_map)

    cmap = copy.copy(plt.get_cmap(cmap_name))
    cmap.set_bad("white")

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if label == "full":
        image = ax.pcolormesh(
            X, Y, density_plot, cmap=cmap, norm=norm,
            shading="auto", rasterized=True,
        )
    else:
        image = ax.imshow(
            density_plot,
            origin="lower",
            extent=[bounds[0], bounds[1], bounds[2], bounds[3]],
            cmap=cmap,
            norm=norm,
            interpolation=interpolation,
            aspect="equal",
        )


    # 表示範囲と交差する実セル境界だけを描画する。
    if grid_segments.size:
        xmin, xmax, ymin, ymax = bounds
        sxmin = np.min(grid_segments[:, :, 0], axis=1)
        sxmax = np.max(grid_segments[:, :, 0], axis=1)
        symin = np.min(grid_segments[:, :, 1], axis=1)
        symax = np.max(grid_segments[:, :, 1], axis=1)
        visible = (
            (sxmax >= xmin) & (sxmin <= xmax)
            & (symax >= ymin) & (symin <= ymax)
        )
        if label == "full":
            grid_linewidth = FULL_GRID_LINEWIDTH
            grid_alpha = FULL_GRID_ALPHA
            grid_color = FULL_GRID_COLOR
        else:
            grid_linewidth = ZOOM_GRID_LINEWIDTH
            grid_alpha = ZOOM_GRID_ALPHA
            grid_color = ZOOM_GRID_COLOR
        collection = LineCollection(
            grid_segments[visible],
            colors=grid_color,
            linewidths=grid_linewidth,
            alpha=grid_alpha,
            zorder=3,
            rasterized=True,
        )
        ax.add_collection(collection)

    # シンク内部では密度だけでなくグリッド線も白抜きにする。
    if sink_radius_au > 0.0:
        ax.add_patch(Circle(
            (0.0, 0.0), sink_radius_au,
            facecolor="white", edgecolor="none", zorder=4,
        ))

    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect("equal")
    ax.axhline(0.0, color="gray", lw=0.6, alpha=0.55)
    ax.axvline(0.0, color="gray", lw=0.6, alpha=0.55)

    scale, exponent = scientific_axis_scale(bounds)
    formatter = FuncFormatter(lambda value, position: f"{value / scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\,{{\rm AU}}]$")
    ax.set_ylabel(rf"$y\ [10^{{{exponent}}}\,{{\rm AU}}]$")

    time_text = f"{time_kyr:.{TIME_DECIMALS}f}"
    title_suffix = {
        "full": "full domain",
        "zoom_1over10": "1/10 zoom",
    }[label]
    ax.set_title(
        "Density map with actual cell grid on the x-y midplane\n"
        f"step={step:05d}, t={time_text} kyr, {title_suffix}",
        fontsize=14,
    )

    colorbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    colorbar.ax.yaxis.set_major_formatter(LogFormatterSciNotation())
    colorbar.set_label(r"$\rho\ [M_\odot\, {\rm AU}^{-3}]$")

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    if SHOW_IMAGES:
        plt.show()
    plt.close(fig)
    del fig, density_plot, image


# ============================================================
# VTK ファイルの整理
# ============================================================
files_by_step = defaultdict(list)
for path in vtk_dir.glob("*.vtk"):
    if "Toyouchi.block" not in path.name:
        continue
    step = get_step_number(path)
    if step is not None:
        files_by_step[step].append(str(path))

steps = sorted(files_by_step)
for step in steps:
    files_by_step[step].sort()
if not steps:
    raise FileNotFoundError(f"No Athena++ block VTK files found in {vtk_dir}")

test_grid = pv.read(files_by_step[steps[0]][0])
density_name = find_density_name(test_grid)
del test_grid

full_bounds = domain_from_step(files_by_step[steps[0]])
plot_bounds = {
    "full": full_bounds,
    "zoom_1over10": fractional_bounds(full_bounds, 1.0 / 10.0),
}

sink_radius_au = float(SINK_RADIUS_AU)
if not np.isfinite(sink_radius_au) or sink_radius_au < 0.0:
    raise ValueError("SINK_RADIUS_AU must be a finite value >= 0 [AU].")

print(f"[INFO] Input          : {vtk_dir}")
print(f"[INFO] Steps          : {len(steps)}")
print(f"[INFO] Density field  : {density_name}")
print(f"[INFO] Full bounds AU : {full_bounds}")
for key, directory in output_dirs.items():
    print(f"[INFO] Output {key:13s}: {directory}")


# ============================================================
# 全体図用の共通カラースケールをサンプル時刻から決定
# シンク内部は統計から除外する。
# ============================================================
sample_indices = np.unique(
    np.linspace(0, len(steps) - 1, min(GLOBAL_SCALE_SAMPLE_COUNT, len(steps))).astype(int)
)
sample_density = []
for sample_index in sample_indices:
    step = steps[sample_index]
    points, density, _ = load_step_slice(
        files_by_step[step], density_name, collect_grid=False
    )
    outside_sink = np.hypot(points[:, 0], points[:, 1]) > sink_radius_au
    if np.any(outside_sink):
        sample_density.append(density[outside_sink])
    del points, density, outside_sink
    gc.collect()

full_vmin, full_vmax = log_limits(
    np.concatenate(sample_density), FULL_PERCENTILES
)
full_norm = LogNorm(vmin=full_vmin, vmax=full_vmax)
del sample_density
gc.collect()

print(f"[INFO] Full color range: {full_vmin:.6e} -- {full_vmax:.6e}")


# ============================================================
# 各 step を1つずつ読み込み、描画後すぐ解放
# ============================================================
for index, step in enumerate(steps, start=1):
    files = files_by_step[step]
    points, density, grid_segments = load_step_slice(
        files, density_name, collect_grid=True
    )
    time_code = read_vtk_time_code(files[0])
    time_kyr = time_code * TIME_UNIT_KYR

    radius_xy = np.hypot(points[:, 0], points[:, 1])
    outside_sink = radius_xy > sink_radius_au

    print(
        f"[{index:4d}/{len(steps):4d}] step={step:05d}, "
        f"t={time_kyr:.{TIME_DECIMALS}f} kyr, cells={len(density)}"
    )

    for label in ("full", "zoom_1over10"):
        bounds = plot_bounds[label]
        resolution = FULL_RESOLUTION if label == "full" else ZOOM_RESOLUTION
        X, Y, density_map, source_cells = interpolate_density(
            points, density, bounds, resolution
        )

        if label == "full":
            norm = full_norm
            cmap_name = FULL_CMAP
            interpolation = "none"
        else:
            xmin, xmax, ymin, ymax = bounds
            in_view = (
                (points[:, 0] >= xmin) & (points[:, 0] <= xmax)
                & (points[:, 1] >= ymin) & (points[:, 1] <= ymax)
                & outside_sink
            )
            zoom_vmin, zoom_vmax = log_limits(
                density[in_view],
                ZOOM_PERCENTILES,
                fallback=(full_vmin, full_vmax),
            )
            norm = LogNorm(vmin=zoom_vmin, vmax=zoom_vmax)
            cmap_name = ZOOM_CMAP
            interpolation = ZOOM_INTERPOLATION

        output_path = output_dirs[label] / (
            f"density_xy_step_{step:05d}_{label}.png"
        )
        plot_and_save(
            X=X,
            Y=Y,
            density_map=density_map,
            bounds=bounds,
            norm=norm,
            cmap_name=cmap_name,
            interpolation=interpolation,
            step=step,
            time_kyr=time_kyr,
            sink_radius_au=sink_radius_au,
            source_cells=source_cells,
            grid_segments=grid_segments,
            label=label,
            output_path=output_path,
        )
        print(f"    saved: {output_path}")
        del X, Y, density_map, norm
        gc.collect()

    del points, density, grid_segments, radius_xy, outside_sink
    gc.collect()

print("[DONE] All density maps were saved.")




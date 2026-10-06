# -*- coding: utf-8 -*-
# Converted from プラズマβ.ipynb

# %% cell 1
# ============================================================
# 全タイムステップ：x-z平面 プラズマβマップ
#   ・y=0ミッドプレーン
#   ・全計算領域 + 各軸幅1/10 zoom
#   ・距離：AU（例 x [10^4 AU]）
#   ・時間：Myr
#   ・等温音速 cs=0.707（コード単位、L0 = rb / 2）
#   ・Athena++規格：Pmag=B^2/2
#
# 必要パッケージ：numpy, scipy, matplotlib, pyvista
# ============================================================

import gc
import os
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, LogLocator


# ============================================================
# 1. ユーザー設定
# ============================================================

# VTKファイルが置かれているディレクトリに変更する
vtk_dir = Path(
    os.path.expanduser("~/athena-project/results/〇〇")
).resolve()

output_dir_full = vtk_dir / "xz_plasma_beta_full"
output_dir_zoom10 = vtk_dir / "xz_plasma_beta_zoom_1over10"
output_dir_full.mkdir(parents=True, exist_ok=True)
output_dir_zoom10.mkdir(parents=True, exist_ok=True)

# 等温音速（コード単位）
CS_CODE = 0.707

# βの共通表示範囲。添付図と同じく 10^-1 ～ 10^2
# 自動範囲にしたい場合は両方をNoneにする
BETA_VMIN = None
BETA_VMAX = None

# カラーマップ：低β=濃紫、高β=黄
BETA_CMAP = "plasma"

# 補間後の画像解像度
MAP_N = 700

# 1/10 zoom：x幅、z幅を全域の1/10にする
ZOOM_FRACTION = 1.0 / 10.0

# y=0を中心とする断面
Y_SLICE_CODE = 0.0

# 図と保存の設定
FIGSIZE = (8.5, 8.0)
SAVE_DPI = 200

# 表示だけ行い保存しない場合はFalse
SAVE_IMAGES = True

# Notebookにも画像を表示する場合はTrue
# 全stepでTrueにすると大量に表示されるので通常はFalseを推奨
SHOW_IMAGES = False


# ============================================================
# 2. コード単位
# ============================================================

M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15  # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
YEAR_CGS = 365.25 * 24.0 * 3600.0

LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
TIME_UNIT_MYR = T_UNIT_CGS / YEAR_CGS / 1.0e6

print("=== Plasma-beta x-z maps ===")
print(f"Input          : {vtk_dir}")
print(f"Full output    : {output_dir_full}")
print(f"1/10 output    : {output_dir_zoom10}")
print(f"1 code length  : {LENGTH_UNIT_AU:.6e} AU")
print(f"1 code time    : {TIME_UNIT_MYR:.6e} Myr")
print(f"isothermal cs  : {CS_CODE:.6e} [code velocity]")

if not vtk_dir.exists():
    raise FileNotFoundError(f"VTKディレクトリが見つかりません: {vtk_dir}")


# ============================================================
# 3. VTKファイルを出力stepごとに整理
# ============================================================

def get_step_number(path):
    """ファイル名末尾の .out2.00000.vtk からstepを取得する。"""
    match = re.search(r"\.out2\.(\d+)\.vtk$", path.name)
    return int(match.group(1)) if match else None


def read_vtk_time_code(path):
    """Athena++ legacy VTKのヘッダからcode timeを読む。"""
    with open(path, "rb") as handle:
        header = handle.read(2048).decode("ascii", errors="ignore")

    match = re.search(
        r"time\s*=\s*([+\-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+\-]?\d+)?)",
        header,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTKヘッダからtimeを取得できません: {path}")
    return float(match.group(1))


vtk_files = sorted(vtk_dir.glob("Toyouchi.block*.out2.*.vtk"))

if not vtk_files:
    raise FileNotFoundError(
        "VTKファイルが見つかりません。\n"
        f"検索条件: {vtk_dir / 'Toyouchi.block*.out2.*.vtk'}"
    )

files_by_step = defaultdict(list)
for path in vtk_files:
    step = get_step_number(path)
    if step is not None:
        files_by_step[step].append(path)

steps = sorted(files_by_step)
if not steps:
    raise RuntimeError("有効なout2 step番号を取得できませんでした。")

for step in steps:
    files_by_step[step].sort()

print(f"Detected VTK files : {len(vtk_files)}")
print(f"Detected steps     : {len(steps)} ({steps[0]:05d} -- {steps[-1]:05d})")


# ============================================================
# 4. VTK物理量名の検索
# ============================================================

def find_cell_array_name(grid, candidates):
    lower_to_original = {
        name.lower(): name for name in grid.cell_data.keys()
    }
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    return None


RHO_CANDIDATES = ("rho", "dens", "density", "prim_dens", "prim_density")
B_CANDIDATES = ("Bcc", "bcc", "B", "magnetic_field")


# ============================================================
# 5. y=0近傍のx-z断面を抽出
# ============================================================

def load_xz_slice(step):
    """
    y=Y_SLICE_CODEと交差するleaf MeshBlockのみを使用し、各blockで
    その面に最も近いセル中心層を抽出する。

    y=0を挟む等距離の2層がある場合は両方を採用し、後段で同一(x,z)
    座標ごとに平均する。
    """
    x_parts, z_parts = [], []
    rho_parts, B_parts = [], []

    domain_xmin = np.inf
    domain_xmax = -np.inf
    domain_zmin = np.inf
    domain_zmax = -np.inf

    files = files_by_step[step]
    time_code = read_vtk_time_code(files[0])

    intersecting_blocks = 0

    for path in files:
        grid = pv.read(path)
        bounds = grid.bounds  # xmin,xmax,ymin,ymax,zmin,zmax

        ymin, ymax = bounds[2], bounds[3]
        tolerance = 1.0e-12 * max(
            abs(ymin), abs(ymax), abs(Y_SLICE_CODE), 1.0
        )

        if not (ymin - tolerance <= Y_SLICE_CODE <= ymax + tolerance):
            del grid
            continue

        intersecting_blocks += 1
        domain_xmin = min(domain_xmin, bounds[0])
        domain_xmax = max(domain_xmax, bounds[1])
        domain_zmin = min(domain_zmin, bounds[4])
        domain_zmax = max(domain_zmax, bounds[5])

        rho_name = find_cell_array_name(grid, RHO_CANDIDATES)
        B_name = find_cell_array_name(grid, B_CANDIDATES)

        if rho_name is None or B_name is None:
            available = list(grid.cell_data.keys())
            raise KeyError(
                f"rhoまたはBccがありません: {path.name}\n"
                f"cell_data={available}"
            )

        centers = np.asarray(grid.cell_centers().points)
        rho = np.asarray(grid.cell_data[rho_name]).reshape(-1)
        B = np.asarray(grid.cell_data[B_name])

        if B.ndim != 2 or B.shape[1] < 3:
            raise ValueError(f"Bccのshapeが想定外です: {B.shape}, file={path.name}")
        B = B[:, :3]

        y = centers[:, 1]
        distance = np.abs(y - Y_SLICE_CODE)
        nearest_distance = np.min(distance)
        layer_tolerance = max(1.0e-12, nearest_distance * 1.0e-6)
        plane_mask = np.abs(distance - nearest_distance) <= layer_tolerance

        valid = (
            plane_mask
            & np.isfinite(rho)
            & (rho > 0.0)
            & np.all(np.isfinite(B), axis=1)
        )

        if np.any(valid):
            x_parts.append(centers[valid, 0])
            z_parts.append(centers[valid, 2])
            rho_parts.append(rho[valid])
            B_parts.append(B[valid])

        del centers, rho, B, grid

    if not x_parts:
        raise RuntimeError(f"step={step:05d}からy=0近傍セルを抽出できません。")

    return {
        "step": step,
        "time_code": time_code,
        "x": np.concatenate(x_parts),
        "z": np.concatenate(z_parts),
        "rho": np.concatenate(rho_parts),
        "B": np.vstack(B_parts),
        "domain_xmin": domain_xmin,
        "domain_xmax": domain_xmax,
        "domain_zmin": domain_zmin,
        "domain_zmax": domain_zmax,
        "intersecting_blocks": intersecting_blocks,
    }


# ============================================================
# 6. y=0上下層など、同一(x,z)座標の値を平均
# ============================================================

def average_duplicate_xz(data):
    coords = np.column_stack((data["x"], data["z"]))
    unique_coords, inverse = np.unique(coords, axis=0, return_inverse=True)
    counts = np.bincount(inverse).astype(float)

    data["x"] = unique_coords[:, 0]
    data["z"] = unique_coords[:, 1]
    data["rho"] = np.bincount(inverse, weights=data["rho"]) / counts

    B_average = np.empty((len(unique_coords), 3), dtype=float)
    for component in range(3):
        B_average[:, component] = (
            np.bincount(inverse, weights=data["B"][:, component]) / counts
        )
    data["B"] = B_average
    return data


# ============================================================
# 7. プラズマβを計算
# ============================================================

def calculate_plasma_beta(rho_code, B_code):
    """
    Athena++で磁場がsqrt(4π)吸収型の場合：
        Pgas = rho * cs^2              （等温）
        Pmag = (Bx^2+By^2+Bz^2) / 2
        beta = Pgas/Pmag = 2*rho*cs^2/B^2

    βは無次元なのでCGS変換は不要。
    """
    B_squared = np.sum(B_code**2, axis=1)
    beta = np.full_like(rho_code, np.nan, dtype=float)

    valid = (
        np.isfinite(rho_code)
        & (rho_code > 0.0)
        & np.isfinite(B_squared)
        & (B_squared > 0.0)
    )
    beta[valid] = 2.0 * rho_code[valid] * CS_CODE**2 / B_squared[valid]
    return beta


# ============================================================
# 8. 自動β表示範囲（BETA_VMIN/VMaxがNoneの場合のみ）
# ============================================================

def determine_beta_limits():
    if BETA_VMIN is not None and BETA_VMAX is not None:
        return float(BETA_VMIN), float(BETA_VMAX)

    # 全セルを保持せず、各stepの対数βヒストグラムだけを加算する
    log_edges = np.linspace(-12.0, 12.0, 1201)
    histogram = np.zeros(len(log_edges) - 1, dtype=np.int64)

    print("[INFO] Determining common beta range...")
    for index, step in enumerate(steps, start=1):
        data = average_duplicate_xz(load_xz_slice(step))
        beta = calculate_plasma_beta(data["rho"], data["B"])
        log_beta = np.log10(beta[np.isfinite(beta) & (beta > 0.0)])
        if log_beta.size:
            histogram += np.histogram(log_beta, bins=log_edges)[0]
        del data, beta, log_beta
        gc.collect()
        print(f"  beta scan {index:3d}/{len(steps):3d}: step={step:05d}")

    cumulative = np.cumsum(histogram)
    total = cumulative[-1]
    if total == 0:
        raise RuntimeError("有限かつ正のプラズマβがありません。")

    def histogram_percentile(percent):
        target = percent / 100.0 * total
        bin_index = np.searchsorted(cumulative, target)
        bin_index = np.clip(bin_index, 0, len(log_edges) - 2)
        return 10.0 ** (0.5 * (log_edges[bin_index] + log_edges[bin_index + 1]))

    automatic_min = histogram_percentile(1.0)
    automatic_max = histogram_percentile(99.0)
    vmin = automatic_min if BETA_VMIN is None else float(BETA_VMIN)
    vmax = automatic_max if BETA_VMAX is None else float(BETA_VMAX)
    return vmin, vmax


beta_vmin, beta_vmax = determine_beta_limits()
if not (np.isfinite(beta_vmin) and np.isfinite(beta_vmax) and 0 < beta_vmin < beta_vmax):
    raise ValueError(f"不正なβ表示範囲です: {beta_vmin}, {beta_vmax}")

print(f"Common beta range : {beta_vmin:.6e} -- {beta_vmax:.6e}")


# ============================================================
# 9. 軸を x [10^n AU] の形式にする
# ============================================================

def apply_scaled_au_axes(ax, xmin, xmax, zmin, zmax):
    maximum = max(abs(xmin), abs(xmax), abs(zmin), abs(zmax))
    exponent = int(np.floor(np.log10(maximum))) if maximum > 0.0 else 0
    scale = 10.0**exponent
    formatter = FuncFormatter(lambda value, position: f"{value / scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\ {{\rm AU}}]$")
    ax.set_ylabel(rf"$z\ [10^{{{exponent}}}\ {{\rm AU}}]$")


# ============================================================
# 10. 1枚のβマップを作る
# ============================================================

def plot_beta_map(data, zoom_fraction, output_dir, suffix):
    step = data["step"]
    time_myr = data["time_code"] * TIME_UNIT_MYR

    x_au = data["x"] * LENGTH_UNIT_AU
    z_au = data["z"] * LENGTH_UNIT_AU
    beta = calculate_plasma_beta(data["rho"], data["B"])

    full_xmin = data["domain_xmin"] * LENGTH_UNIT_AU
    full_xmax = data["domain_xmax"] * LENGTH_UNIT_AU
    full_zmin = data["domain_zmin"] * LENGTH_UNIT_AU
    full_zmax = data["domain_zmax"] * LENGTH_UNIT_AU

    if zoom_fraction is None:
        xmin, xmax = full_xmin, full_xmax
        zmin, zmax = full_zmin, full_zmax
        range_name = "full domain"
    else:
        x_center = 0.0 if full_xmin <= 0.0 <= full_xmax else 0.5 * (full_xmin + full_xmax)
        z_center = 0.0 if full_zmin <= 0.0 <= full_zmax else 0.5 * (full_zmin + full_zmax)
        x_half = 0.5 * (full_xmax - full_xmin) * zoom_fraction
        z_half = 0.5 * (full_zmax - full_zmin) * zoom_fraction
        xmin, xmax = x_center - x_half, x_center + x_half
        zmin, zmax = z_center - z_half, z_center + z_half
        range_name = f"1/{int(round(1.0 / zoom_fraction))} zoom"

    margin_x = 0.03 * (xmax - xmin)
    margin_z = 0.03 * (zmax - zmin)
    source_mask = (
        (x_au >= xmin - margin_x) & (x_au <= xmax + margin_x)
        & (z_au >= zmin - margin_z) & (z_au <= zmax + margin_z)
        & np.isfinite(beta) & (beta > 0.0)
    )

    source_count = np.count_nonzero(source_mask)
    if source_count < 4:
        raise RuntimeError(
            f"step={step:05d}, {range_name}: 補間元セルが不足しています "
            f"({source_count} cells)。"
        )

    x_grid = np.linspace(xmin, xmax, MAP_N)
    z_grid = np.linspace(zmin, zmax, MAP_N)
    X, Z = np.meshgrid(x_grid, z_grid)
    points = np.column_stack((x_au[source_mask], z_au[source_mask]))

    # βを直接線形補間せず、log10(beta)を補間する
    log_beta_source = np.log10(beta[source_mask])
    try:
        log_beta_grid = griddata(points, log_beta_source, (X, Z), method="linear")
    except Exception as error:
        print(f"[WARN] linear interpolation failed ({error}); nearest only")
        log_beta_grid = np.full_like(X, np.nan, dtype=float)

    if np.any(~np.isfinite(log_beta_grid)):
        nearest = griddata(points, log_beta_source, (X, Z), method="nearest")
        log_beta_grid = np.where(np.isfinite(log_beta_grid), log_beta_grid, nearest)

    beta_grid = 10.0**log_beta_grid

    fig, ax = plt.subplots(figsize=FIGSIZE)
    image = ax.pcolormesh(
        X,
        Z,
        beta_grid,
        shading="auto",
        cmap=BETA_CMAP,
        norm=LogNorm(vmin=beta_vmin, vmax=beta_vmax),
        rasterized=True,
    )

    colorbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    colorbar.set_label(r"Plasma $\beta=P_{\rm gas}/P_{\rm mag}$")
    colorbar.locator = LogLocator(base=10.0)
    colorbar.update_ticks()

    # β=1（ガス圧と磁気圧が等しい）の境界を細い黒線で表示
    if beta_vmin < 1.0 < beta_vmax:
        finite_grid = beta_grid[np.isfinite(beta_grid)]
        if finite_grid.size and np.min(finite_grid) <= 1.0 <= np.max(finite_grid):
            ax.contour(X, Z, beta_grid, levels=[1.0], colors="black", linewidths=0.7)

    ax.axhline(0.0, color="white", lw=0.5, alpha=0.35)
    ax.axvline(0.0, color="white", lw=0.5, alpha=0.35)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(zmin, zmax)
    ax.set_aspect("equal", adjustable="box")
    apply_scaled_au_axes(ax, xmin, xmax, zmin, zmax)

    ax.set_title(
        r"Plasma $\beta$ on the x-z plane" + "\n"
        f"step={step:05d}, t={time_myr:.6e} Myr, {range_name}"
    )

    ax.text(
        0.02,
        0.98,
        f"source cells = {source_count}\n"
        f"β range = {beta_vmin:.1e} -- {beta_vmax:.1e}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox=dict(boxstyle="round", facecolor="black", edgecolor="white", alpha=0.55),
    )

    fig.tight_layout()
    output_path = output_dir / f"xz_plasma_beta_{step:05d}_{suffix}.png"

    if SAVE_IMAGES:
        fig.savefig(output_path, dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
        print(f"[SAVED] {output_path}")

    if SHOW_IMAGES:
        plt.show()
    plt.close(fig)

    del X, Z, beta_grid, log_beta_grid, beta
    return output_path


# ============================================================
# 11. 全stepを逐次処理（stepごとにメモリ解放）
# ============================================================

saved_full = []
saved_zoom10 = []

for index, step in enumerate(steps, start=1):
    print(f"\n[STEP {index:3d}/{len(steps):3d}] step={step:05d}")

    slice_data = average_duplicate_xz(load_xz_slice(step))
    print(
        f"  t={slice_data['time_code'] * TIME_UNIT_MYR:.6e} Myr, "
        f"blocks={slice_data['intersecting_blocks']}, "
        f"unique cells={len(slice_data['x'])}"
    )

    saved_full.append(
        plot_beta_map(
            slice_data,
            zoom_fraction=None,
            output_dir=output_dir_full,
            suffix="full",
        )
    )

    saved_zoom10.append(
        plot_beta_map(
            slice_data,
            zoom_fraction=ZOOM_FRACTION,
            output_dir=output_dir_zoom10,
            suffix="zoom10",
        )
    )

    del slice_data
    gc.collect()

print("\n=== Completed ===")
print(f"Full-domain images : {len(saved_full)} -> {output_dir_full}")
print(f"1/10 zoom images   : {len(saved_zoom10)} -> {output_dir_zoom10}")


# %% cell 2
# ============================================================
# 全タイムステップ：x-z平面 プラズマβ + 投影磁力線マップ
#   ・y=0ミッドプレーン
#   ・全計算領域 + 各軸幅1/10 zoom
#   ・距離：AU（例 x [10^4 AU]）
#   ・時間：Myr
#   ・等温音速 cs=0.707（コード単位、L0 = rb / 2）
#   ・Athena++規格：Pmag=B^2/2
#   ・x-z面への投影磁場(Bx,Bz)を細い黒線で重ねる
#
# 必要パッケージ：numpy, scipy, matplotlib, pyvista
# ============================================================

import gc
import os
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, LogLocator


# ============================================================
# 1. ユーザー設定
# ============================================================

# VTKファイルが置かれているディレクトリに変更する
vtk_dir = Path(
    os.path.expanduser("~/athena-project/results/〇〇")
).resolve()

output_dir_full = vtk_dir / "xz_plasma_beta_with_Blines_full"
output_dir_zoom10 = vtk_dir / "xz_plasma_beta_with_Blines_zoom_1over10"
output_dir_full.mkdir(parents=True, exist_ok=True)
output_dir_zoom10.mkdir(parents=True, exist_ok=True)

# 等温音速（コード単位）
CS_CODE = 0.707

# βの共通表示範囲：全stepの1--99 percentileから自動決定
BETA_VMIN = None
BETA_VMAX = None

# カラーマップ：低β=濃紫、高β=黄
BETA_CMAP = "plasma"

# 投影磁力線の表示設定
STREAM_DENSITY = 1.1
STREAM_LINEWIDTH = 0.45
STREAM_ARROWSIZE = 0.55

# 補間後の画像解像度
MAP_N = 700

# 1/10 zoom：x幅、z幅を全域の1/10にする
ZOOM_FRACTION = 1.0 / 10.0

# y=0を中心とする断面
Y_SLICE_CODE = 0.0

# 図と保存の設定
FIGSIZE = (8.5, 8.0)
SAVE_DPI = 200

# 表示だけ行い保存しない場合はFalse
SAVE_IMAGES = True

# Notebookにも画像を表示する場合はTrue
# 全stepでTrueにすると大量に表示されるので通常はFalseを推奨
SHOW_IMAGES = False


# ============================================================
# 2. コード単位
# ============================================================

M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15  # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
YEAR_CGS = 365.25 * 24.0 * 3600.0

LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
TIME_UNIT_MYR = T_UNIT_CGS / YEAR_CGS / 1.0e6

print("=== Plasma-beta x-z maps ===")
print(f"Input          : {vtk_dir}")
print(f"Full output    : {output_dir_full}")
print(f"1/10 output    : {output_dir_zoom10}")
print(f"1 code length  : {LENGTH_UNIT_AU:.6e} AU")
print(f"1 code time    : {TIME_UNIT_MYR:.6e} Myr")
print(f"isothermal cs  : {CS_CODE:.6e} [code velocity]")

if not vtk_dir.exists():
    raise FileNotFoundError(f"VTKディレクトリが見つかりません: {vtk_dir}")


# ============================================================
# 3. VTKファイルを出力stepごとに整理
# ============================================================

def get_step_number(path):
    """ファイル名末尾の .out2.00000.vtk からstepを取得する。"""
    match = re.search(r"\.out2\.(\d+)\.vtk$", path.name)
    return int(match.group(1)) if match else None


def read_vtk_time_code(path):
    """Athena++ legacy VTKのヘッダからcode timeを読む。"""
    with open(path, "rb") as handle:
        header = handle.read(2048).decode("ascii", errors="ignore")

    match = re.search(
        r"time\s*=\s*([+\-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+\-]?\d+)?)",
        header,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTKヘッダからtimeを取得できません: {path}")
    return float(match.group(1))


vtk_files = sorted(vtk_dir.glob("Toyouchi.block*.out2.*.vtk"))

if not vtk_files:
    raise FileNotFoundError(
        "VTKファイルが見つかりません。\n"
        f"検索条件: {vtk_dir / 'Toyouchi.block*.out2.*.vtk'}"
    )

files_by_step = defaultdict(list)
for path in vtk_files:
    step = get_step_number(path)
    if step is not None:
        files_by_step[step].append(path)

steps = sorted(files_by_step)
if not steps:
    raise RuntimeError("有効なout2 step番号を取得できませんでした。")

for step in steps:
    files_by_step[step].sort()

print(f"Detected VTK files : {len(vtk_files)}")
print(f"Detected steps     : {len(steps)} ({steps[0]:05d} -- {steps[-1]:05d})")


# ============================================================
# 4. VTK物理量名の検索
# ============================================================

def find_cell_array_name(grid, candidates):
    lower_to_original = {
        name.lower(): name for name in grid.cell_data.keys()
    }
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    return None


RHO_CANDIDATES = ("rho", "dens", "density", "prim_dens", "prim_density")
B_CANDIDATES = ("Bcc", "bcc", "B", "magnetic_field")


# ============================================================
# 5. y=0近傍のx-z断面を抽出
# ============================================================

def load_xz_slice(step):
    """
    y=Y_SLICE_CODEと交差するleaf MeshBlockのみを使用し、各blockで
    その面に最も近いセル中心層を抽出する。

    y=0を挟む等距離の2層がある場合は両方を採用し、後段で同一(x,z)
    座標ごとに平均する。
    """
    x_parts, z_parts = [], []
    rho_parts, B_parts = [], []

    domain_xmin = np.inf
    domain_xmax = -np.inf
    domain_zmin = np.inf
    domain_zmax = -np.inf

    files = files_by_step[step]
    time_code = read_vtk_time_code(files[0])

    intersecting_blocks = 0

    for path in files:
        grid = pv.read(path)
        bounds = grid.bounds  # xmin,xmax,ymin,ymax,zmin,zmax

        ymin, ymax = bounds[2], bounds[3]
        tolerance = 1.0e-12 * max(
            abs(ymin), abs(ymax), abs(Y_SLICE_CODE), 1.0
        )

        if not (ymin - tolerance <= Y_SLICE_CODE <= ymax + tolerance):
            del grid
            continue

        intersecting_blocks += 1
        domain_xmin = min(domain_xmin, bounds[0])
        domain_xmax = max(domain_xmax, bounds[1])
        domain_zmin = min(domain_zmin, bounds[4])
        domain_zmax = max(domain_zmax, bounds[5])

        rho_name = find_cell_array_name(grid, RHO_CANDIDATES)
        B_name = find_cell_array_name(grid, B_CANDIDATES)

        if rho_name is None or B_name is None:
            available = list(grid.cell_data.keys())
            raise KeyError(
                f"rhoまたはBccがありません: {path.name}\n"
                f"cell_data={available}"
            )

        centers = np.asarray(grid.cell_centers().points)
        rho = np.asarray(grid.cell_data[rho_name]).reshape(-1)
        B = np.asarray(grid.cell_data[B_name])

        if B.ndim != 2 or B.shape[1] < 3:
            raise ValueError(f"Bccのshapeが想定外です: {B.shape}, file={path.name}")
        B = B[:, :3]

        y = centers[:, 1]
        distance = np.abs(y - Y_SLICE_CODE)
        nearest_distance = np.min(distance)
        layer_tolerance = max(1.0e-12, nearest_distance * 1.0e-6)
        plane_mask = np.abs(distance - nearest_distance) <= layer_tolerance

        valid = (
            plane_mask
            & np.isfinite(rho)
            & (rho > 0.0)
            & np.all(np.isfinite(B), axis=1)
        )

        if np.any(valid):
            x_parts.append(centers[valid, 0])
            z_parts.append(centers[valid, 2])
            rho_parts.append(rho[valid])
            B_parts.append(B[valid])

        del centers, rho, B, grid

    if not x_parts:
        raise RuntimeError(f"step={step:05d}からy=0近傍セルを抽出できません。")

    return {
        "step": step,
        "time_code": time_code,
        "x": np.concatenate(x_parts),
        "z": np.concatenate(z_parts),
        "rho": np.concatenate(rho_parts),
        "B": np.vstack(B_parts),
        "domain_xmin": domain_xmin,
        "domain_xmax": domain_xmax,
        "domain_zmin": domain_zmin,
        "domain_zmax": domain_zmax,
        "intersecting_blocks": intersecting_blocks,
    }


# ============================================================
# 6. y=0上下層など、同一(x,z)座標の値を平均
# ============================================================

def average_duplicate_xz(data):
    coords = np.column_stack((data["x"], data["z"]))
    unique_coords, inverse = np.unique(coords, axis=0, return_inverse=True)
    counts = np.bincount(inverse).astype(float)

    data["x"] = unique_coords[:, 0]
    data["z"] = unique_coords[:, 1]
    data["rho"] = np.bincount(inverse, weights=data["rho"]) / counts

    B_average = np.empty((len(unique_coords), 3), dtype=float)
    for component in range(3):
        B_average[:, component] = (
            np.bincount(inverse, weights=data["B"][:, component]) / counts
        )
    data["B"] = B_average
    return data


# ============================================================
# 7. プラズマβを計算
# ============================================================

def calculate_plasma_beta(rho_code, B_code):
    """
    Athena++で磁場がsqrt(4π)吸収型の場合：
        Pgas = rho * cs^2              （等温）
        Pmag = (Bx^2+By^2+Bz^2) / 2
        beta = Pgas/Pmag = 2*rho*cs^2/B^2

    βは無次元なのでCGS変換は不要。
    """
    B_squared = np.sum(B_code**2, axis=1)
    beta = np.full_like(rho_code, np.nan, dtype=float)

    valid = (
        np.isfinite(rho_code)
        & (rho_code > 0.0)
        & np.isfinite(B_squared)
        & (B_squared > 0.0)
    )
    beta[valid] = 2.0 * rho_code[valid] * CS_CODE**2 / B_squared[valid]
    return beta


# ============================================================
# 8. 自動β表示範囲（BETA_VMIN/VMaxがNoneの場合のみ）
# ============================================================

def determine_beta_limits():
    if BETA_VMIN is not None and BETA_VMAX is not None:
        return float(BETA_VMIN), float(BETA_VMAX)

    # 全セルを保持せず、各stepの対数βヒストグラムだけを加算する
    log_edges = np.linspace(-12.0, 12.0, 1201)
    histogram = np.zeros(len(log_edges) - 1, dtype=np.int64)

    print("[INFO] Determining common beta range...")
    for index, step in enumerate(steps, start=1):
        data = average_duplicate_xz(load_xz_slice(step))
        beta = calculate_plasma_beta(data["rho"], data["B"])
        log_beta = np.log10(beta[np.isfinite(beta) & (beta > 0.0)])
        if log_beta.size:
            histogram += np.histogram(log_beta, bins=log_edges)[0]
        del data, beta, log_beta
        gc.collect()
        print(f"  beta scan {index:3d}/{len(steps):3d}: step={step:05d}")

    cumulative = np.cumsum(histogram)
    total = cumulative[-1]
    if total == 0:
        raise RuntimeError("有限かつ正のプラズマβがありません。")

    def histogram_percentile(percent):
        target = percent / 100.0 * total
        bin_index = np.searchsorted(cumulative, target)
        bin_index = np.clip(bin_index, 0, len(log_edges) - 2)
        return 10.0 ** (0.5 * (log_edges[bin_index] + log_edges[bin_index + 1]))

    automatic_min = histogram_percentile(1.0)
    automatic_max = histogram_percentile(99.0)
    vmin = automatic_min if BETA_VMIN is None else float(BETA_VMIN)
    vmax = automatic_max if BETA_VMAX is None else float(BETA_VMAX)
    return vmin, vmax


beta_vmin, beta_vmax = determine_beta_limits()
if not (np.isfinite(beta_vmin) and np.isfinite(beta_vmax) and 0 < beta_vmin < beta_vmax):
    raise ValueError(f"不正なβ表示範囲です: {beta_vmin}, {beta_vmax}")

print(f"Common beta range : {beta_vmin:.6e} -- {beta_vmax:.6e}")


# ============================================================
# 9. 軸を x [10^n AU] の形式にする
# ============================================================

def apply_scaled_au_axes(ax, xmin, xmax, zmin, zmax):
    maximum = max(abs(xmin), abs(xmax), abs(zmin), abs(zmax))
    exponent = int(np.floor(np.log10(maximum))) if maximum > 0.0 else 0
    scale = 10.0**exponent
    formatter = FuncFormatter(lambda value, position: f"{value / scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\ {{\rm AU}}]$")
    ax.set_ylabel(rf"$z\ [10^{{{exponent}}}\ {{\rm AU}}]$")


# ============================================================
# 10. 1枚のβマップを作る
# ============================================================

def plot_beta_map(data, zoom_fraction, output_dir, suffix):
    step = data["step"]
    time_myr = data["time_code"] * TIME_UNIT_MYR

    x_au = data["x"] * LENGTH_UNIT_AU
    z_au = data["z"] * LENGTH_UNIT_AU
    beta = calculate_plasma_beta(data["rho"], data["B"])

    full_xmin = data["domain_xmin"] * LENGTH_UNIT_AU
    full_xmax = data["domain_xmax"] * LENGTH_UNIT_AU
    full_zmin = data["domain_zmin"] * LENGTH_UNIT_AU
    full_zmax = data["domain_zmax"] * LENGTH_UNIT_AU

    if zoom_fraction is None:
        xmin, xmax = full_xmin, full_xmax
        zmin, zmax = full_zmin, full_zmax
        range_name = "full domain"
    else:
        x_center = 0.0 if full_xmin <= 0.0 <= full_xmax else 0.5 * (full_xmin + full_xmax)
        z_center = 0.0 if full_zmin <= 0.0 <= full_zmax else 0.5 * (full_zmin + full_zmax)
        x_half = 0.5 * (full_xmax - full_xmin) * zoom_fraction
        z_half = 0.5 * (full_zmax - full_zmin) * zoom_fraction
        xmin, xmax = x_center - x_half, x_center + x_half
        zmin, zmax = z_center - z_half, z_center + z_half
        range_name = f"1/{int(round(1.0 / zoom_fraction))} zoom"

    margin_x = 0.03 * (xmax - xmin)
    margin_z = 0.03 * (zmax - zmin)
    source_mask = (
        (x_au >= xmin - margin_x) & (x_au <= xmax + margin_x)
        & (z_au >= zmin - margin_z) & (z_au <= zmax + margin_z)
        & np.isfinite(beta) & (beta > 0.0)
    )

    source_count = np.count_nonzero(source_mask)
    if source_count < 4:
        raise RuntimeError(
            f"step={step:05d}, {range_name}: 補間元セルが不足しています "
            f"({source_count} cells)。"
        )

    x_grid = np.linspace(xmin, xmax, MAP_N)
    z_grid = np.linspace(zmin, zmax, MAP_N)
    X, Z = np.meshgrid(x_grid, z_grid)
    points = np.column_stack((x_au[source_mask], z_au[source_mask]))

    # βを直接線形補間せず、log10(beta)を補間する
    log_beta_source = np.log10(beta[source_mask])
    try:
        log_beta_grid = griddata(points, log_beta_source, (X, Z), method="linear")
    except Exception as error:
        print(f"[WARN] linear interpolation failed ({error}); nearest only")
        log_beta_grid = np.full_like(X, np.nan, dtype=float)

    if np.any(~np.isfinite(log_beta_grid)):
        nearest = griddata(points, log_beta_source, (X, Z), method="nearest")
        log_beta_grid = np.where(np.isfinite(log_beta_grid), log_beta_grid, nearest)

    beta_grid = 10.0**log_beta_grid

    # --------------------------------------------------------
    # x-z面への投影磁場 (Bx, Bz) を同じ規則格子へ補間
    # 磁力線は方向を示すため、単位変換は不要
    # --------------------------------------------------------
    Bx_source = data["B"][source_mask, 0]
    Bz_source = data["B"][source_mask, 2]

    def interpolate_field_component(values):
        try:
            linear = griddata(points, values, (X, Z), method="linear")
        except Exception:
            linear = np.full_like(X, np.nan, dtype=float)

        if np.any(~np.isfinite(linear)):
            nearest = griddata(points, values, (X, Z), method="nearest")
            linear = np.where(np.isfinite(linear), linear, nearest)

        return np.nan_to_num(linear, nan=0.0, posinf=0.0, neginf=0.0)

    Bx_grid = interpolate_field_component(Bx_source)
    Bz_grid = interpolate_field_component(Bz_source)
    Bproj_grid = np.hypot(Bx_grid, Bz_grid)

    finite_Bproj = Bproj_grid[np.isfinite(Bproj_grid)]
    has_projected_field = (
        finite_Bproj.size > 0
        and np.max(finite_Bproj) > 0.0
    )

    fig, ax = plt.subplots(figsize=FIGSIZE)
    image = ax.pcolormesh(
        X,
        Z,
        beta_grid,
        shading="auto",
        cmap=BETA_CMAP,
        norm=LogNorm(vmin=beta_vmin, vmax=beta_vmax),
        rasterized=True,
    )

    # 細い黒線：x-z面に投影した磁力線
    if has_projected_field:
        ax.streamplot(
            x_grid,
            z_grid,
            Bx_grid,
            Bz_grid,
            color="black",
            density=STREAM_DENSITY,
            linewidth=STREAM_LINEWIDTH,
            arrowsize=STREAM_ARROWSIZE,
            arrowstyle="->",
            minlength=0.03,
            maxlength=20.0,
            integration_direction="both",
            broken_streamlines=False,
            zorder=3,
        )

    colorbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    colorbar.set_label(r"Plasma $\beta=P_{\rm gas}/P_{\rm mag}$")
    colorbar.locator = LogLocator(base=10.0)
    colorbar.update_ticks()

    # β=1（ガス圧と磁気圧が等しい）の境界。
    # 黒い磁力線と区別するため白い破線にする。
    if beta_vmin < 1.0 < beta_vmax:
        finite_grid = beta_grid[np.isfinite(beta_grid)]
        if finite_grid.size and np.min(finite_grid) <= 1.0 <= np.max(finite_grid):
            ax.contour(
                X,
                Z,
                beta_grid,
                levels=[1.0],
                colors="white",
                linewidths=0.75,
                linestyles="--",
                zorder=4,
            )

    ax.axhline(0.0, color="white", lw=0.5, alpha=0.35)
    ax.axvline(0.0, color="white", lw=0.5, alpha=0.35)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(zmin, zmax)
    ax.set_aspect("equal", adjustable="box")
    apply_scaled_au_axes(ax, xmin, xmax, zmin, zmax)

    ax.set_title(
        r"Plasma $\beta$ and projected magnetic field lines" + "\n"
        f"step={step:05d}, t={time_myr:.6e} Myr, {range_name}"
    )

    ax.text(
        0.02,
        0.98,
        f"source cells = {source_count}\n"
        f"β range = {beta_vmin:.1e} -- {beta_vmax:.1e}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox=dict(boxstyle="round", facecolor="black", edgecolor="white", alpha=0.55),
    )

    fig.tight_layout()
    output_path = output_dir / f"xz_plasma_beta_Blines_{step:05d}_{suffix}.png"

    if SAVE_IMAGES:
        fig.savefig(output_path, dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
        print(f"[SAVED] {output_path}")

    if SHOW_IMAGES:
        plt.show()
    plt.close(fig)

    del X, Z, beta_grid, log_beta_grid, beta
    del Bx_grid, Bz_grid, Bproj_grid, Bx_source, Bz_source
    return output_path


# ============================================================
# 11. 全stepを逐次処理（stepごとにメモリ解放）
# ============================================================

saved_full = []
saved_zoom10 = []

for index, step in enumerate(steps, start=1):
    print(f"\n[STEP {index:3d}/{len(steps):3d}] step={step:05d}")

    slice_data = average_duplicate_xz(load_xz_slice(step))
    print(
        f"  t={slice_data['time_code'] * TIME_UNIT_MYR:.6e} Myr, "
        f"blocks={slice_data['intersecting_blocks']}, "
        f"unique cells={len(slice_data['x'])}"
    )

    saved_full.append(
        plot_beta_map(
            slice_data,
            zoom_fraction=None,
            output_dir=output_dir_full,
            suffix="full",
        )
    )

    saved_zoom10.append(
        plot_beta_map(
            slice_data,
            zoom_fraction=ZOOM_FRACTION,
            output_dir=output_dir_zoom10,
            suffix="zoom10",
        )
    )

    del slice_data
    gc.collect()

print("\n=== Completed ===")
print(f"Full-domain images : {len(saved_full)} -> {output_dir_full}")
print(f"1/10 zoom images   : {len(saved_zoom10)} -> {output_dir_zoom10}")


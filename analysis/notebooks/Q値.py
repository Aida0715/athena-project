# -*- coding: utf-8 -*-
"""Toyouchi SMR VTKからface-on Toomre-Qマップを全時刻分作成する。

定義:
  Sigma(x,y) = integral rho dz  (|z| <= Z_HALF_WIDTH_AU)
  Omega(R)   = sum(m_i vphi_i/R_i) / sum(m_i)
  kappa^2    = 2 Omega^2 [2 + d ln(Omega)/d ln(R)]
  Q(x,y)     = cs kappa(R) / [pi G Sigma(x,y)]

SMRセルはセル体積で質量加重し、面密度マップは各セルの投影矩形を
2次元difference arrayへ加えるため、点補間より質量保存性が高い
中心付近はΩ=vφ/Rが不安定になるため、ソフトニング長または最小セル幅を基準に除外しているため白抜きになる
また所々リング状に白抜きになるのはΩが半径方向に急落し、κ^2が負になりQ値を定義できない場所を除外しているため
（リング状になるのはある半径でΩを円環平均しているため）
本解析コードでは局所的な面密度と、円環平均した面密度がありそれぞれ、Q値map、半径プロファイルに使用されている
"""

from pathlib import Path
from collections import defaultdict
import csv
import gc
import json
import re

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter1d


# ============================================================
# 1. ユーザー設定
# ============================================================
DATADIR = Path("~/athena-project/results/〇〇").expanduser().resolve()
VTK_GLOB = "Toyouchi.block*.out2.*.vtk"
STEPS = None                    # None: 全時刻、例: [0, 25, 50, 100]

# True: metadataからraccを取得し、中心のr < raccを白抜きにする
MASK_RACC = True

# 共通親ディレクトリ。この下に出力種類別ディレクトリを作る。
OUTPUT_ROOT = DATADIR / "toomre_Q_analysis"

# 面密度の積分範囲 |z| <= 2500 AU
Z_HALF_WIDTH_AU = 2500.0

# face-on表示幅。0.1なら全計算領域の各軸幅の1/10。
MAP_FRACTION = 0.1
MAP_RESOLUTION = 700
MAP_AXIS_SCALE_AU = 1.0e4

# 外側でOmega/kappaが定義できない場合、値を外挿せず、有効半径に
# 内接する正方形へ表示を自動クロップする。
CROP_MAP_TO_VALID_Q = True

# 円環プロファイル設定
N_RADIAL_BINS = 60
MIN_CELLS_PER_ANNULUS = 8
INNER_EXCLUSION_CELL_WIDTHS = 2.0

# kappaはOmegaの微分を含み、SMR境界・円環サンプリングの微小な段差を
# 増幅しやすい。等間隔log(R)ビン上でOmegaとその傾きを別々に平滑化する。
OMEGA_SMOOTH_SIGMA_BINS = 3.0
SLOPE_SMOOTH_SIGMA_BINS = 2.0

# Q表示と円盤外マスク
Q_DISPLAY_MIN = 0.1
Q_DISPLAY_MAX = 10.0
SIGMA_MASK_FRACTION = 1.0e-4  # 各時刻のmax(Sigma)に対する比

# Q=1線はSMRセル境界を多数拾って格子状に見える場合があるため既定OFF。
# 必要な診断時だけTrueにする。NPZには未平滑Qを常に保存する。
DRAW_Q1_CONTOUR = False
MAP_INTERPOLATION = "bilinear"

FIGSIZE_SINGLE = (8.2, 7.0)
FIGSIZE_SUMMARY = (13.0, 10.0)
DPI = 180
OVERWRITE = False
SHOW_IMAGES = False

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0
KYR_CGS = 1.0e3 * YEAR_CGS


# ============================================================
# 2. 出力先とmetadata
# ============================================================
OUTPUT_DIRS = {
    "surface_density": OUTPUT_ROOT / "surface_density",
    "q_map": OUTPUT_ROOT / "Q_maps",
    "rotation": OUTPUT_ROOT / "rotation_profiles",
    "summary": OUTPUT_ROOT / "summary_4panel",
    "tables": OUTPUT_ROOT / "tables",
}
for directory in OUTPUT_DIRS.values():
    directory.mkdir(parents=True, exist_ok=True)


def load_metadata(directory):
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(
        directory.glob("*.metadata.json")
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one metadata JSON in {directory}; "
            f"found {len(candidates)}"
        )
    metadata = json.loads(candidates[0].read_text(encoding="utf-8"))
    if metadata.get("format") != "Toyouchi runtime metadata":
        raise ValueError(f"Unsupported metadata: {candidates[0]}")
    metadata["_path"] = str(candidates[0])
    return metadata


METADATA = load_metadata(DATADIR)
UNITS = METADATA["units"]
LENGTH_CGS = float(UNITS["length_cgs"])
MASS_CGS = float(UNITS["mass_cgs"])
TIME_CGS = float(UNITS["time_cgs"])
LENGTH_AU = LENGTH_CGS / AU_CGS
SIGMA_UNIT_MSUN_AU2 = (
    (MASS_CGS / MSUN_CGS) / (LENGTH_AU**2)
)
OMEGA_UNIT_KYR_INV = KYR_CGS / TIME_CGS
CS_CODE = float(METADATA["hydro"]["sound_speed_code"])
G_CODE = float(UNITS["gravity_constant_code"])
EPSILON_CODE = float(METADATA["gravity"]["softening_code"])
MSTAR_CODE = float(METADATA["gravity"]["mstar_code"])
CENTRAL_GRAVITY = bool(METADATA["gravity"]["central_enabled"])
Z_HALF_WIDTH_CODE = Z_HALF_WIDTH_AU / LENGTH_AU
RACC_AU = (
    float(METADATA["virtual_sphere"]["radius_au"])
    if MASK_RACC else 0.0
)
BZ0_MICROGAUSS = float(METADATA["magnetic_field"]["initial_bz_microgauss"])
BZ0_LABEL = rf"$B_{{z0}}={BZ0_MICROGAUSS:g}\,\mu{{\rm G}}$"

mesh_meta = METADATA["mesh"]
FULL_BOUNDS = np.array([
    [mesh_meta["x1min_code"], mesh_meta["x1max_code"]],
    [mesh_meta["x2min_code"], mesh_meta["x2max_code"]],
    [mesh_meta["x3min_code"], mesh_meta["x3max_code"]],
], dtype=float)
CENTER = np.zeros(3, dtype=float)

if not DATADIR.is_dir():
    raise FileNotFoundError(DATADIR)
if not (0.0 < MAP_FRACTION <= 1.0):
    raise ValueError("MAP_FRACTION must be in (0, 1]")
if Z_HALF_WIDTH_AU <= 0.0:
    raise ValueError("Z_HALF_WIDTH_AU must be positive")


# ============================================================
# 3. VTK整理・基本情報
# ============================================================
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
    match = re.search(
        r"\btime\s*=\s*([-+0-9.eEdD]+)", header, flags=re.IGNORECASE
    )
    if not match:
        raise ValueError(f"VTK time was not found in {path}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def block_geometry(grid, path):
    bounds = np.asarray(grid.bounds, dtype=float).reshape(3, 2)
    try:
        cell_dims = np.asarray(grid.dimensions, dtype=int) - 1
    except Exception as exc:
        raise ValueError(f"Structured-grid dimensions unavailable: {path}") from exc
    if np.any(cell_dims <= 0) or int(np.prod(cell_dims)) != grid.n_cells:
        raise ValueError(
            f"Unexpected dimensions at {path}: dimensions={grid.dimensions}, "
            f"n_cells={grid.n_cells}"
        )
    spacing = (bounds[:, 1] - bounds[:, 0]) / cell_dims
    if np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError(f"Invalid cell spacing at {path}: {spacing}")
    return bounds, spacing


def map_bounds():
    full_width = FULL_BOUNDS[:2, 1] - FULL_BOUNDS[:2, 0]
    half_width = 0.5 * MAP_FRACTION * full_width
    return np.array([
        [CENTER[0] - half_width[0], CENTER[0] + half_width[0]],
        [CENTER[1] - half_width[1], CENTER[1] + half_width[1]],
    ])


# ============================================================
# 4. 1時刻のSMR集計
# ============================================================
def inspect_resolution(paths):
    minimum = np.full(3, np.inf)
    for path in paths:
        grid = pv.read(path)
        try:
            _, spacing = block_geometry(grid, path)
            minimum = np.minimum(minimum, spacing)
        finally:
            del grid
    if np.any(~np.isfinite(minimum)):
        raise RuntimeError("Failed to determine the finest cell size")
    return minimum


def choose_radial_edges(bounds_xy, finest_spacing):
    # 表示正方形の四隅までkappaを与える。計算全領域は表示領域より
    # 広いため、この半径までの円環セルをVTK全体から集計できる。
    half_x = max(
        bounds_xy[0, 1] - CENTER[0], CENTER[0] - bounds_xy[0, 0]
    )
    half_y = max(
        bounds_xy[1, 1] - CENTER[1], CENTER[1] - bounds_xy[1, 0]
    )
    map_rmax = np.hypot(half_x, half_y)
    inner = max(
        EPSILON_CODE,
        INNER_EXCLUSION_CELL_WIDTHS * max(finest_spacing[0], finest_spacing[1]),
    )
    if not 0.0 < inner < map_rmax:
        raise ValueError(
            f"Invalid radial interval: inner={inner}, outer={map_rmax} [code]"
        )
    return np.geomspace(inner, map_rmax, N_RADIAL_BINS + 1)


def add_rectangles(diff, ix0, ix1, iy0, iy1, values):
    np.add.at(diff, (iy0, ix0), values)
    np.add.at(diff, (iy1, ix0), -values)
    np.add.at(diff, (iy0, ix1), -values)
    np.add.at(diff, (iy1, ix1), values)


def aggregate_step(paths):
    bounds_xy = map_bounds()
    finest = inspect_resolution(paths)
    radial_edges = choose_radial_edges(bounds_xy, finest)
    radial_centers = np.sqrt(radial_edges[:-1] * radial_edges[1:])

    nx = ny = int(MAP_RESOLUTION)
    x_edges = np.linspace(bounds_xy[0, 0], bounds_xy[0, 1], nx + 1)
    y_edges = np.linspace(bounds_xy[1, 0], bounds_xy[1, 1], ny + 1)
    pixel_dx = x_edges[1] - x_edges[0]
    pixel_dy = y_edges[1] - y_edges[0]
    sigma_diff = np.zeros((ny + 1, nx + 1), dtype=np.float64)

    ann_mass = np.zeros(N_RADIAL_BINS, dtype=np.float64)
    ann_omega_mass = np.zeros(N_RADIAL_BINS, dtype=np.float64)
    ann_count = np.zeros(N_RADIAL_BINS, dtype=np.int64)

    for block_number, path in enumerate(paths, start=1):
        grid = pv.read(path)
        try:
            missing = [name for name in ("rho", "vel") if name not in grid.cell_data]
            if missing:
                raise KeyError(f"{path}: missing cell arrays {missing}")
            rho = np.asarray(grid.cell_data["rho"], dtype=np.float64).reshape(-1)
            velocity = np.asarray(grid.cell_data["vel"], dtype=np.float64)
            xyz = np.asarray(grid.cell_centers().points, dtype=np.float64)
            bounds, spacing = block_geometry(grid, path)
            dx, dy, dz = spacing

            if velocity.shape != (grid.n_cells, 3) or rho.size != grid.n_cells:
                raise ValueError(f"Unexpected rho/vel shape: {path}")
            finite = np.isfinite(rho) & (rho > 0.0) & np.all(np.isfinite(velocity), axis=1)

            zlo = np.maximum(xyz[:, 2] - 0.5 * dz, -Z_HALF_WIDTH_CODE)
            zhi = np.minimum(xyz[:, 2] + 0.5 * dz, Z_HALF_WIDTH_CODE)
            dz_overlap = np.maximum(zhi - zlo, 0.0)
            selected = finite & (dz_overlap > 0.0)
            if not np.any(selected):
                continue

            x = xyz[selected, 0] - CENTER[0]
            y = xyz[selected, 1] - CENTER[1]
            rho_s = rho[selected]
            vel_s = velocity[selected]
            dz_s = dz_overlap[selected]

            # 円環量：SMRセルの実際の部分体積で質量加重する。
            radius = np.hypot(x, y)
            mass = rho_s * dx * dy * dz_s
            vphi = np.divide(
                -y * vel_s[:, 0] + x * vel_s[:, 1], radius,
                out=np.zeros_like(radius), where=radius > 0.0,
            )
            omega = np.divide(
                vphi, radius, out=np.zeros_like(radius), where=radius > 0.0
            )
            ibin = np.searchsorted(radial_edges, radius, side="right") - 1
            # 逆回転成分も除外せず、符号付きOmegaを質量平均する。
            valid_bin = (
                (ibin >= 0) & (ibin < N_RADIAL_BINS)
                & np.isfinite(omega)
            )
            np.add.at(ann_mass, ibin[valid_bin], mass[valid_bin])
            np.add.at(
                ann_omega_mass, ibin[valid_bin],
                mass[valid_bin] * omega[valid_bin],
            )
            np.add.at(ann_count, ibin[valid_bin], 1)

            # 面密度：各セルの投影矩形へ rho*dz_overlap を加える。
            cell_xlo = x - 0.5 * dx
            cell_xhi = x + 0.5 * dx
            cell_ylo = y - 0.5 * dy
            cell_yhi = y + 0.5 * dy
            in_map = (
                (cell_xhi > bounds_xy[0, 0]) & (cell_xlo < bounds_xy[0, 1])
                & (cell_yhi > bounds_xy[1, 0]) & (cell_ylo < bounds_xy[1, 1])
            )
            if np.any(in_map):
                ix0 = np.floor((cell_xlo[in_map] - bounds_xy[0, 0]) / pixel_dx).astype(int)
                ix1 = np.ceil((cell_xhi[in_map] - bounds_xy[0, 0]) / pixel_dx).astype(int)
                iy0 = np.floor((cell_ylo[in_map] - bounds_xy[1, 0]) / pixel_dy).astype(int)
                iy1 = np.ceil((cell_yhi[in_map] - bounds_xy[1, 0]) / pixel_dy).astype(int)
                ix0 = np.clip(ix0, 0, nx); ix1 = np.clip(ix1, 0, nx)
                iy0 = np.clip(iy0, 0, ny); iy1 = np.clip(iy1, 0, ny)
                nonempty = (ix1 > ix0) & (iy1 > iy0)
                add_rectangles(
                    sigma_diff,
                    ix0[nonempty], ix1[nonempty],
                    iy0[nonempty], iy1[nonempty],
                    (rho_s[in_map] * dz_s[in_map])[nonempty],
                )
        finally:
            del grid
        if block_number % 128 == 0:
            gc.collect()

    sigma_map = np.cumsum(np.cumsum(sigma_diff[:-1, :-1], axis=0), axis=1)
    ann_area = np.pi * (radial_edges[1:]**2 - radial_edges[:-1]**2)
    sigma_annulus = np.divide(
        ann_mass, ann_area, out=np.full_like(ann_mass, np.nan), where=ann_area > 0
    )
    omega_raw = np.divide(
        ann_omega_mass, ann_mass, out=np.full_like(ann_mass, np.nan),
        where=(ann_mass > 0.0) & (ann_count >= MIN_CELLS_PER_ANNULUS),
    )
    return {
        "bounds_xy": bounds_xy, "x_edges": x_edges, "y_edges": y_edges,
        "sigma_map": sigma_map, "radial_edges": radial_edges,
        "radius": radial_centers, "sigma_annulus": sigma_annulus,
        "omega_raw": omega_raw, "ann_count": ann_count, "finest": finest,
    }


# ============================================================
# 5. Omega平滑化、kappa、Q
# ============================================================
def smooth_rotation_and_kappa(radius, omega_raw):
    valid = np.isfinite(radius) & np.isfinite(omega_raw) & (omega_raw > 0.0)
    if np.count_nonzero(valid) < 5:
        raise RuntimeError("Too few valid annuli to calculate kappa")

    log_r_all = np.log(radius)
    log_omega_interp = np.interp(
        log_r_all, log_r_all[valid], np.log(omega_raw[valid]),
        left=np.nan, right=np.nan,
    )
    usable = np.isfinite(log_omega_interp)
    log_r_usable = log_r_all[usable]
    raw_log = log_omega_interp[usable]
    smooth_log = gaussian_filter1d(
        raw_log, sigma=OMEGA_SMOOTH_SIGMA_BINS,
        mode="nearest", truncate=3.0,
    )

    omega_smooth = np.full_like(radius, np.nan)
    omega_smooth[usable] = np.exp(smooth_log)
    slope = np.full_like(radius, np.nan)
    slope_raw = np.gradient(smooth_log, log_r_usable)
    slope[usable] = gaussian_filter1d(
        slope_raw, sigma=SLOPE_SMOOTH_SIGMA_BINS,
        mode="nearest", truncate=3.0,
    )
    kappa2 = 2.0 * omega_smooth**2 * (2.0 + slope)
    kappa = np.full_like(kappa2, np.nan)
    stable = np.isfinite(kappa2) & (kappa2 > 0.0)
    kappa[stable] = np.sqrt(kappa2[stable])
    return omega_smooth, kappa, slope


def make_q_products(result):
    radius = result["radius"]
    omega_smooth, kappa, slope = smooth_rotation_and_kappa(
        radius, result["omega_raw"]
    )
    q_annulus = CS_CODE * kappa / (np.pi * G_CODE * result["sigma_annulus"])

    x_centers = 0.5 * (result["x_edges"][:-1] + result["x_edges"][1:])
    y_centers = 0.5 * (result["y_edges"][:-1] + result["y_edges"][1:])
    xx, yy = np.meshgrid(x_centers - CENTER[0], y_centers - CENTER[1])
    rr = np.hypot(xx, yy)
    valid_kappa = np.isfinite(kappa) & (kappa > 0.0)
    kappa_map = np.full_like(rr, np.nan)
    # kappa^2 <= 0 または未定義の円環を跨いで補間すると、本来Qを
    # 定義できない半径に人工的なkappaを作り、同心円状の偽構造になる。
    # したがって連続した有効区間ごとにだけ補間し、欠損円環はNaNのまま残す。
    valid_indices = np.flatnonzero(valid_kappa)
    if valid_indices.size >= 2:
        split_at = np.flatnonzero(np.diff(valid_indices) > 1) + 1
        valid_runs = np.split(valid_indices, split_at)
        for run in valid_runs:
            if run.size < 2:
                continue
            rlo = radius[run[0]]
            rhi = radius[run[-1]]
            inside = (rr >= rlo) & (rr <= rhi)
            kappa_map[inside] = np.interp(
                np.log(rr[inside]), np.log(radius[run]), kappa[run]
            )
    sigma = result["sigma_map"]
    sigma_threshold = SIGMA_MASK_FRACTION * np.nanmax(sigma)
    disk = np.isfinite(sigma) & (sigma > sigma_threshold) & np.isfinite(kappa_map)
    q_map = np.full_like(sigma, np.nan)
    q_map[disk] = CS_CODE * kappa_map[disk] / (np.pi * G_CODE * sigma[disk])

    result.update({
        "omega_smooth": omega_smooth, "kappa": kappa, "slope": slope,
        "q_annulus": q_annulus, "q_map": q_map,
        "sigma_threshold": sigma_threshold,
        "q_valid_rmax": (
            float(radius[valid_kappa].max())
            if np.any(valid_kappa) else np.nan
        ),
    })
    return result


# ============================================================
# 6. 保存と描画
# ============================================================
def axis_extent_au(result):
    bounds = result["bounds_xy"] * LENGTH_AU
    return [
        bounds[0, 0] / MAP_AXIS_SCALE_AU,
        bounds[0, 1] / MAP_AXIS_SCALE_AU,
        bounds[1, 0] / MAP_AXIS_SCALE_AU,
        bounds[1, 1] / MAP_AXIS_SCALE_AU,
    ]


def save_table(step, time_kyr, result):
    path = OUTPUT_DIRS["tables"] / f"Q_profile_{step:05d}.csv"
    omega_k = np.full_like(result["radius"], np.nan)
    if CENTRAL_GRAVITY and MSTAR_CODE > 0.0:
        omega_k = np.sqrt(
            G_CODE * MSTAR_CODE
            / (result["radius"]**2 + EPSILON_CODE**2)**1.5
        )
    columns = np.column_stack([
        result["radius"] * LENGTH_AU,
        result["sigma_annulus"] * SIGMA_UNIT_MSUN_AU2,
        result["omega_raw"] * OMEGA_UNIT_KYR_INV,
        result["omega_smooth"] * OMEGA_UNIT_KYR_INV,
        omega_k * OMEGA_UNIT_KYR_INV,
        result["kappa"] * OMEGA_UNIT_KYR_INV,
        result["q_annulus"], result["ann_count"],
    ])
    header = (
        "R_AU,Sigma_Msun_AU-2,Omega_raw_kyr-1,Omega_smooth_kyr-1,"
        "OmegaK_softened_kyr-1,kappa_kyr-1,Q_annulus,cell_count\n"
        f"step={step}, time_kyr={time_kyr:.8e}, "
        f"z_half_width_AU={Z_HALF_WIDTH_AU:.8e}"
    )
    np.savetxt(path, columns, delimiter=",", header=header)
    np.savez_compressed(
        OUTPUT_DIRS["tables"] / f"Q_maps_{step:05d}.npz",
        x_edges_au=result["x_edges"] * LENGTH_AU,
        y_edges_au=result["y_edges"] * LENGTH_AU,
        sigma_Msun_AU2=result["sigma_map"] * SIGMA_UNIT_MSUN_AU2,
        q_map=result["q_map"], time_kyr=time_kyr,
    )


def draw_sigma(ax, result):
    sigma_msun_au2 = result["sigma_map"] * SIGMA_UNIT_MSUN_AU2
    if MASK_RACC:
        x = 0.5 * (result["x_edges"][:-1] + result["x_edges"][1:]) * LENGTH_AU
        y = 0.5 * (result["y_edges"][:-1] + result["y_edges"][1:]) * LENGTH_AU
        xx, yy = np.meshgrid(x, y)
        sigma_msun_au2 = np.where(xx*xx + yy*yy < RACC_AU**2, np.nan, sigma_msun_au2)
    positive = sigma_msun_au2[
        np.isfinite(sigma_msun_au2) & (sigma_msun_au2 > 0.0)
    ]
    if positive.size == 0:
        raise RuntimeError("Surface-density map has no positive values")
    vmin, vmax = np.nanpercentile(positive, [1.0, 99.8])
    if not vmax > vmin:
        vmin, vmax = positive.min(), positive.max() * 1.001
    image = ax.imshow(
        sigma_msun_au2, origin="lower", extent=axis_extent_au(result),
        cmap="magma", norm=LogNorm(vmin=max(vmin, np.finfo(float).tiny), vmax=vmax),
        interpolation=MAP_INTERPOLATION, aspect="equal",
    )
    return image


def draw_q(ax, result):
    q_map = result["q_map"]
    if MASK_RACC:
        x = 0.5 * (result["x_edges"][:-1] + result["x_edges"][1:]) * LENGTH_AU
        y = 0.5 * (result["y_edges"][:-1] + result["y_edges"][1:]) * LENGTH_AU
        xx, yy = np.meshgrid(x, y)
        q_map = np.where(xx*xx + yy*yy < RACC_AU**2, np.nan, q_map)
    log_q = np.log10(q_map)
    image = ax.imshow(
        log_q, origin="lower", extent=axis_extent_au(result), cmap="RdBu",
        # log10(Q)=0、すなわちQ=1を発散色の白い中心に固定する。
        norm=TwoSlopeNorm(
            vmin=np.log10(Q_DISPLAY_MIN), vcenter=0.0,
            vmax=np.log10(Q_DISPLAY_MAX),
        ),
        interpolation=MAP_INTERPOLATION, aspect="equal",
    )
    finite = q_map[np.isfinite(q_map)]
    if (
        DRAW_Q1_CONTOUR and finite.size
        and finite.min() <= 1.0 <= finite.max()
    ):
        x = (
            0.5 * (result["x_edges"][:-1] + result["x_edges"][1:])
            * LENGTH_AU / MAP_AXIS_SCALE_AU
        )
        y = (
            0.5 * (result["y_edges"][:-1] + result["y_edges"][1:])
            * LENGTH_AU / MAP_AXIS_SCALE_AU
        )
        ax.contour(x, y, q_map, levels=[1.0], colors="black", linewidths=1.0)
    return image


def format_map_axis(ax, result):
    exponent = int(np.log10(MAP_AXIS_SCALE_AU))
    ax.set_xlabel(rf"$x$ [$10^{{{exponent}}}$ AU]")
    ax.set_ylabel(rf"$y$ [$10^{{{exponent}}}$ AU]")
    if CROP_MAP_TO_VALID_Q and np.isfinite(result["q_valid_rmax"]):
        # 円形のkappa有効領域に完全に入る最大の軸平行正方形。
        half = (
            0.995 * result["q_valid_rmax"] / np.sqrt(2.0)
            * LENGTH_AU / MAP_AXIS_SCALE_AU
        )
        full_extent = axis_extent_au(result)
        half = min(
            half,
            abs(full_extent[0]), abs(full_extent[1]),
            abs(full_extent[2]), abs(full_extent[3]),
        )
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)


def plot_outputs(step, time_kyr, result):
    suffix = f"{step:05d}"

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    image = draw_sigma(ax, result); format_map_axis(ax, result)
    ax.set_title(rf"Surface density, $t={time_kyr:.1f}$ kyr, $|z|\leq{Z_HALF_WIDTH_AU:g}$ AU, {BZ0_LABEL}")
    fig.colorbar(
        image, ax=ax,
        label=r"$\Sigma$ [$M_\odot\,{\rm AU}^{-2}$]",
    )
    fig.tight_layout(); fig.savefig(OUTPUT_DIRS["surface_density"] / f"Sigma_{suffix}.png", dpi=DPI)
    if SHOW_IMAGES: plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    image = draw_q(ax, result); format_map_axis(ax, result)
    contour_note = "; black: $Q=1$" if DRAW_Q1_CONTOUR else ""
    ax.set_title(rf"Toomre $Q$, $t={time_kyr:.1f}$ kyr, {BZ0_LABEL}{contour_note}")
    cbar = fig.colorbar(image, ax=ax, label=r"Toomre $Q$")
    ticks = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
    cbar.set_ticks(np.log10(ticks)); cbar.set_ticklabels([f"{v:g}" for v in ticks])
    fig.tight_layout(); fig.savefig(OUTPUT_DIRS["q_map"] / f"Q_{suffix}.png", dpi=DPI)
    if SHOW_IMAGES: plt.show()
    plt.close(fig)

    radius_au = result["radius"] * LENGTH_AU
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 8.0), sharex=True)
    axes[0].loglog(radius_au, result["omega_raw"] * OMEGA_UNIT_KYR_INV, ".", ms=3, alpha=0.55, label="mass-weighted raw")
    axes[0].loglog(radius_au, result["omega_smooth"] * OMEGA_UNIT_KYR_INV, lw=1.7, label="smoothed")
    if CENTRAL_GRAVITY and MSTAR_CODE > 0.0:
        omega_k = np.sqrt(G_CODE * MSTAR_CODE / (result["radius"]**2 + EPSILON_CODE**2)**1.5)
        axes[0].loglog(radius_au, omega_k * OMEGA_UNIT_KYR_INV, "--", lw=1.1, label="softened central-gravity reference")
    axes[0].loglog(radius_au, result["kappa"] * OMEGA_UNIT_KYR_INV, lw=1.4, label=r"$\kappa$")
    axes[0].set_ylabel(r"Frequency [kyr$^{-1}$]"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.25, which="both")
    axes[1].semilogx(radius_au, result["q_annulus"], color="tab:purple", lw=1.5)
    axes[1].axhline(1.0, color="black", ls="--", lw=1.0)
    axes[1].set(xlabel="R [AU]", ylabel="Annular Toomre Q")
    axes[1].grid(alpha=0.25, which="both")
    fig.suptitle(f"Rotation diagnostics: step={step:05d}, t={time_kyr:.1f} kyr, {BZ0_LABEL}")
    fig.tight_layout(); fig.savefig(OUTPUT_DIRS["rotation"] / f"rotation_Q_{suffix}.png", dpi=DPI)
    if SHOW_IMAGES: plt.show()
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE_SUMMARY)
    im0 = draw_sigma(axes[0, 0], result); format_map_axis(axes[0, 0], result); axes[0, 0].set_title(r"$\Sigma(x,y)$")
    fig.colorbar(
        im0, ax=axes[0, 0],
        label=r"$M_\odot\,{\rm AU}^{-2}$", fraction=0.046,
    )
    im1 = draw_q(axes[0, 1], result); format_map_axis(axes[0, 1], result)
    q_title = r"$Q(x,y)$; black: $Q=1$" if DRAW_Q1_CONTOUR else r"$Q(x,y)$"
    axes[0, 1].set_title(q_title)
    cb1 = fig.colorbar(im1, ax=axes[0, 1], fraction=0.046); cb1.set_ticks(np.log10(ticks)); cb1.set_ticklabels([f"{v:g}" for v in ticks])
    axes[1, 0].loglog(radius_au, result["omega_raw"] * OMEGA_UNIT_KYR_INV, ".", ms=3, alpha=0.5, label="raw")
    axes[1, 0].loglog(radius_au, result["omega_smooth"] * OMEGA_UNIT_KYR_INV, label=r"$\Omega$")
    axes[1, 0].loglog(radius_au, result["kappa"] * OMEGA_UNIT_KYR_INV, label=r"$\kappa$")
    axes[1, 0].set(xlabel="R [AU]", ylabel=r"Frequency [kyr$^{-1}$]"); axes[1, 0].legend(); axes[1, 0].grid(alpha=0.25, which="both")
    axes[1, 1].semilogx(radius_au, result["q_annulus"], lw=1.5)
    axes[1, 1].axhline(1.0, color="black", ls="--", lw=1.0)
    axes[1, 1].set(xlabel="R [AU]", ylabel="Annular Q"); axes[1, 1].grid(alpha=0.25, which="both")
    fig.suptitle(rf"Toomre analysis: step={step:05d}, $t={time_kyr:.1f}$ kyr, $|z|\leq{Z_HALF_WIDTH_AU:g}$ AU, {BZ0_LABEL}")
    fig.tight_layout(); fig.savefig(OUTPUT_DIRS["summary"] / f"Q_summary_{suffix}.png", dpi=DPI)
    if SHOW_IMAGES: plt.show()
    plt.close(fig)


# ============================================================
# 7. 全時刻実行
# ============================================================
def main():
    groups = collect_steps()
    print(f"[INFO] Metadata : {METADATA['_path']}")
    print(f"[INFO] Input    : {DATADIR}")
    print(f"[INFO] Output   : {OUTPUT_ROOT}")
    print(f"[INFO] Steps    : {len(groups)}")
    print(f"[INFO] z range  : +/- {Z_HALF_WIDTH_AU:g} AU")
    print(f"[INFO] cs(code) : {CS_CODE:.8e}; G(code)={G_CODE:.8e}")

    for number, (step, paths) in enumerate(groups.items(), start=1):
        sentinel = OUTPUT_DIRS["summary"] / f"Q_summary_{step:05d}.png"
        if sentinel.exists() and not OVERWRITE:
            print(f"[SKIP {number:3d}/{len(groups):3d}] step={step:05d}")
            continue
        print(f"[RUN  {number:3d}/{len(groups):3d}] step={step:05d}, blocks={len(paths)}")
        time_code = read_vtk_time(paths[0])
        result = make_q_products(aggregate_step(paths))
        time_kyr = time_code * TIME_CGS / KYR_CGS
        save_table(step, time_kyr, result)
        plot_outputs(step, time_kyr, result)
        finite_q = result["q_map"][np.isfinite(result["q_map"])]
        if finite_q.size:
            print(
                f"       t={time_kyr:.1f} kyr, "
                f"Q[p16,median,p84]={np.percentile(finite_q, [16,50,84])}"
            )
        del result
        gc.collect()

    print(f"[DONE] {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()

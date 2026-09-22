# -*- coding: utf-8 -*-
# Converted from 磁力線.ipynb

# %% cell 1
#x-z磁場のストリームラインon密度map

import os
import re
import glob
import gc
from collections import defaultdict

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter


# ============================================================
# 1. ユーザー設定
# ============================================================

# VTKファイルがあるディレクトリ
vtk_dir = os.path.expanduser(
    "~/athena-project/results/〇〇"
)

# 親出力ディレクトリ
output_root = os.path.join(
    vtk_dir,
    "xz_density_with_magnetic_fieldlines"
)

# 全体図の保存先
output_dir_full = os.path.join(
    output_root,
    "full"
)

# 1/10ズーム図の保存先
output_dir_zoom10 = os.path.join(
    output_root,
    "zoom_1over10"
)

os.makedirs(output_dir_full, exist_ok=True)
os.makedirs(output_dir_zoom10, exist_ok=True)

# 1/10ズームの各軸幅
zoom_fraction_10 = 1.0 / 10.0

# 補間後の画像解像度
plot_resolution = 700

# 磁力線の密度
stream_density = 1.5

# x-z面の法線（y方向）に平均する半厚み。
# 2.0なら中心面の上下それぞれ約2セル、合計約4セル厚を平均する。
slice_half_thickness_cells = 2.0

# 白抜きにする球対称シンク領域の半径 [AU]
sink_radius_au = 1000.0

# 図のサイズ
figsize = (9, 8)

# 保存画像のdpi
dpi = 200

print("=== Output directories ===")
print(f"Output root : {output_root}")
print(f"Full images : {output_dir_full}")
print(f"1/10 zoom   : {output_dir_zoom10}")


# ============================================================
# 2. Toyouchi.cppのコード単位
# L0 = rb / 2
# ============================================================

Munit = 4.0e33       # g
Lunit = 7.03e15      # cm
Tunit = 3.61e10      # s

Vunit = Lunit / Tunit
Rhounit = Munit / Lunit**3
Punit = Rhounit * Vunit**2

# Athena++の磁場単位 [Gauss]
Bunit = np.sqrt(4.0 * np.pi * Punit)

AU = 1.495978707e13
YEAR = 365.25 * 24.0 * 3600.0

print("\n=== Code units ===")
print(f"Lunit   = {Lunit:.6e} cm")
print(f"Rhounit = {Rhounit:.6e} g cm^-3")
print(f"Bunit   = {Bunit:.6e} G")
print(f"Tunit   = {Tunit / (1.0e6 * YEAR):.6e} Myr")


# ============================================================
# 3. VTKファイルを出力時刻ごとに整理
# ============================================================

timestep_dict = defaultdict(list)

vtk_pattern = os.path.join(
    vtk_dir,
    "Toyouchi.block*.out2.*.vtk"
)

vtk_files = glob.glob(vtk_pattern)

if not vtk_files:
    raise FileNotFoundError(
        "VTK files were not found.\n"
        f"Checked pattern:\n{vtk_pattern}"
    )


def get_block_number(filename):
    match = re.search(
        r"block(\d+)",
        os.path.basename(filename)
    )

    if match is None:
        return -1

    return int(match.group(1))


for filename in vtk_files:
    match = re.search(
        r"\.out2\.(\d+)\.vtk$",
        filename
    )

    if match is not None:
        step = int(match.group(1))
        timestep_dict[step].append(filename)

timesteps = sorted(timestep_dict)

if not timesteps:
    raise RuntimeError(
        "Could not extract timestep numbers "
        "from VTK filenames."
    )

for step in timesteps:
    timestep_dict[step] = sorted(
        timestep_dict[step],
        key=get_block_number
    )

print("\n=== VTK information ===")
print(f"Number of VTK files    : {len(vtk_files)}")
print(f"Number of output steps : {len(timesteps)}")
print(f"First step             : {timesteps[0]:05d}")
print(f"Last step              : {timesteps[-1]:05d}")
print(
    f"Blocks at first step   : "
    f"{len(timestep_dict[timesteps[0]])}"
)
print(
    f"Blocks at last step    : "
    f"{len(timestep_dict[timesteps[-1]])}"
)


# ============================================================
# 4. Athena++ VTKヘッダーから時刻を取得
# ============================================================

def read_athena_vtk_time(filename):
    """Athena++ legacy VTKの2行目からコード時刻を取得する。"""

    with open(filename, "rb") as file:
        file.readline()

        title = file.readline().decode(
            "ascii",
            errors="ignore"
        )

    match = re.search(
        r"time\s*=\s*"
        r"([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        title,
        flags=re.IGNORECASE
    )

    if match is None:
        return np.nan

    return float(match.group(1))


# ============================================================
# 5. VTK配列名の検索
# ============================================================

def find_array_name(grid, candidates):
    """候補の中からVTKに存在する配列名を返す。"""

    for name in candidates:
        if name in grid.array_names:
            return name

    raise KeyError(
        "Required array was not found.\n"
        f"Candidates: {candidates}\n"
        f"Available arrays: {grid.array_names}"
    )


# ============================================================
# 6. x-z断面の抽出
# ============================================================

def load_xz_slice(step):
    """
    y=0面と交差するMeshBlockから、
    各MeshBlockのy方向セル幅を基準に有限厚みのスラブを抽出する。
    同一(x,z)座標に含まれる複数のy層は後段で平均する。
    """

    x_list = []
    z_list = []
    rho_list = []
    B_list = []

    domain_xmin = np.inf
    domain_xmax = -np.inf
    domain_zmin = np.inf
    domain_zmax = -np.inf

    files = timestep_dict[step]
    time_code = read_athena_vtk_time(files[0])

    for filename in files:
        grid = pv.read(filename)

        # xmin, xmax, ymin, ymax, zmin, zmax
        bounds = grid.bounds

        ymin = bounds[2]
        ymax = bounds[3]

        tolerance_bounds = 1.0e-12 * max(
            abs(ymin),
            abs(ymax),
            1.0
        )

        # y=0と交差しないブロックを除外
        if not (
            ymin - tolerance_bounds
            <= 0.0
            <= ymax + tolerance_bounds
        ):
            del grid
            continue

        domain_xmin = min(domain_xmin, bounds[0])
        domain_xmax = max(domain_xmax, bounds[1])
        domain_zmin = min(domain_zmin, bounds[4])
        domain_zmax = max(domain_zmax, bounds[5])

        points = grid.cell_centers().points

        rho_name = find_array_name(
            grid,
            [
                "rho",
                "dens",
                "density",
                "prim_dens",
                "prim_density",
            ]
        )

        B_name = find_array_name(
            grid,
            [
                "Bcc",
                "bcc",
                "B",
                "magnetic_field",
            ]
        )

        rho = np.asarray(
            grid[rho_name]
        ).reshape(-1)

        B = np.asarray(
            grid[B_name]
        ).reshape(-1, 3)

        y = points[:, 1]

        unique_y = np.unique(y)
        dy_values = np.diff(unique_y)
        dy_values = dy_values[dy_values > 0.0]
        if dy_values.size:
            dy_local = np.min(dy_values)
        else:
            dy_local = max(ymax - ymin, 1.0e-12)

        slab_half_thickness = (
            slice_half_thickness_cells * dy_local
        )
        tolerance = max(1.0e-12, 1.0e-10 * dy_local)
        plane_mask = np.abs(y) <= slab_half_thickness + tolerance

        valid_mask = (
            plane_mask
            & np.isfinite(points[:, 0])
            & np.isfinite(points[:, 2])
            & np.isfinite(rho)
            & np.all(np.isfinite(B), axis=1)
        )

        if np.any(valid_mask):
            x_list.append(points[valid_mask, 0])
            z_list.append(points[valid_mask, 2])
            rho_list.append(rho[valid_mask])
            B_list.append(B[valid_mask])

        del points, rho, B, grid

    if not x_list:
        raise RuntimeError(
            "No cells in the finite-thickness slab around y=0 were found.\n"
            f"step={step:05d}"
        )

    return {
        "step": step,
        "time_code": time_code,
        "x_code": np.concatenate(x_list),
        "z_code": np.concatenate(z_list),
        "rho_code": np.concatenate(rho_list),
        "B_code": np.vstack(B_list),
        "domain_xmin_code": domain_xmin,
        "domain_xmax_code": domain_xmax,
        "domain_zmin_code": domain_zmin,
        "domain_zmax_code": domain_zmax,
    }


# ============================================================
# 7. 同一x-z座標にある値を平均
# ============================================================

def average_duplicate_xz(x, z, rho, B):
    """
    y=0まわりの複数層やMeshBlock境界によって同一(x,z)座標に
    複数の値がある場合、それらを平均する。
    """

    coordinates = np.column_stack((x, z))

    unique_coordinates, inverse = np.unique(
        coordinates,
        axis=0,
        return_inverse=True
    )

    counts = np.bincount(
        inverse
    ).astype(float)

    rho_average = (
        np.bincount(inverse, weights=rho)
        / counts
    )

    B_average = np.column_stack([
        np.bincount(inverse, weights=B[:, 0]) / counts,
        np.bincount(inverse, weights=B[:, 1]) / counts,
        np.bincount(inverse, weights=B[:, 2]) / counts,
    ])

    return (
        unique_coordinates[:, 0],
        unique_coordinates[:, 1],
        rho_average,
        B_average
    )


# ============================================================
# 8. 全出力共通の密度カラースケール
# ============================================================

print(
    "\n[INFO] Determining common density color scale..."
)

rho_vmin_candidates = []
rho_vmax_candidates = []

for index, step in enumerate(timesteps):
    slice_data = load_xz_slice(step)

    rho_cgs = (
        slice_data["rho_code"]
        * Rhounit
    )

    valid_density = rho_cgs[
        np.isfinite(rho_cgs)
        & (rho_cgs > 0.0)
    ]

    if valid_density.size > 0:
        rho_vmin_candidates.append(
            np.percentile(valid_density, 0.5)
        )

        rho_vmax_candidates.append(
            np.max(valid_density)
        )

    print(
        f"  [{index + 1:3d}/{len(timesteps):3d}] "
        f"step={step:05d}"
    )

    del slice_data, rho_cgs, valid_density
    gc.collect()

if not rho_vmin_candidates:
    raise RuntimeError(
        "No positive finite density values were found."
    )

rho_vmin = float(
    np.min(rho_vmin_candidates)
)

rho_vmax = float(
    np.max(rho_vmax_candidates)
)

if (
    not np.isfinite(rho_vmin)
    or not np.isfinite(rho_vmax)
    or rho_vmin <= 0.0
    or rho_vmax <= rho_vmin
):
    raise RuntimeError(
        "Invalid density color range:\n"
        f"rho_vmin={rho_vmin}\n"
        f"rho_vmax={rho_vmax}"
    )

print(f"\nrho_vmin = {rho_vmin:.6e} g cm^-3")
print(f"rho_vmax = {rho_vmax:.6e} g cm^-3")


# ============================================================
# 9. 1時刻分の描画・保存
# ============================================================

def plot_xz_density_with_fieldlines(
    step,
    output_dir,
    zoom_fraction=None,
    filename_suffix=""
):
    """
    x-z密度マップへx-z面に射影した磁力線を重ねる。

    zoom_fraction=None:
        実データ全体を表示する。

    zoom_fraction=0.1:
        x方向・z方向とも全体幅の1/10を表示する。
    """

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    data = load_xz_slice(step)

    (
        x_code,
        z_code,
        rho_code,
        B_code
    ) = average_duplicate_xz(
        data["x_code"],
        data["z_code"],
        data["rho_code"],
        data["B_code"]
    )

    # --------------------------------------------------------
    # コード単位から表示単位へ変換
    # --------------------------------------------------------

    x_au = x_code * Lunit / AU
    z_au = z_code * Lunit / AU

    rho_cgs = rho_code * Rhounit

    Bx_microgauss = (
        B_code[:, 0]
        * Bunit
        * 1.0e6
    )

    By_microgauss = (
        B_code[:, 1]
        * Bunit
        * 1.0e6
    )

    Bz_microgauss = (
        B_code[:, 2]
        * Bunit
        * 1.0e6
    )

    Bmag_microgauss = np.sqrt(
        Bx_microgauss**2
        + By_microgauss**2
        + Bz_microgauss**2
    )

    time_kyr = (
        data["time_code"]
        * Tunit
        / (1.0e3 * YEAR)
    )

    # --------------------------------------------------------
    # 実データから全計算領域を取得
    # --------------------------------------------------------

    data_xmin = (
        data["domain_xmin_code"]
        * Lunit
        / AU
    )

    data_xmax = (
        data["domain_xmax_code"]
        * Lunit
        / AU
    )

    data_zmin = (
        data["domain_zmin_code"]
        * Lunit
        / AU
    )

    data_zmax = (
        data["domain_zmax_code"]
        * Lunit
        / AU
    )

    if zoom_fraction is None:
        xmin = data_xmin
        xmax = data_xmax
        zmin = data_zmin
        zmax = data_zmax

        range_label = "full domain"

    else:
        if not (0.0 < zoom_fraction <= 1.0):
            raise ValueError(
                "zoom_fraction must satisfy "
                "0 < zoom_fraction <= 1."
            )

        x_center = (
            0.0
            if data_xmin <= 0.0 <= data_xmax
            else 0.5 * (data_xmin + data_xmax)
        )

        z_center = (
            0.0
            if data_zmin <= 0.0 <= data_zmax
            else 0.5 * (data_zmin + data_zmax)
        )

        x_half_width = (
            0.5
            * (data_xmax - data_xmin)
            * zoom_fraction
        )

        z_half_width = (
            0.5
            * (data_zmax - data_zmin)
            * zoom_fraction
        )

        xmin = max(
            data_xmin,
            x_center - x_half_width
        )

        xmax = min(
            data_xmax,
            x_center + x_half_width
        )

        zmin = max(
            data_zmin,
            z_center - z_half_width
        )

        zmax = min(
            data_zmax,
            z_center + z_half_width
        )

        range_label = (
            f"zoom={zoom_fraction:g} "
            "of full width"
        )

    print(
        f"[RANGE] step={step:05d}, {range_label}\n"
        f"        full x={data_xmin:.6e} -- "
        f"{data_xmax:.6e} AU\n"
        f"        full z={data_zmin:.6e} -- "
        f"{data_zmax:.6e} AU\n"
        f"        plot x={xmin:.6e} -- "
        f"{xmax:.6e} AU\n"
        f"        plot z={zmin:.6e} -- "
        f"{zmax:.6e} AU"
    )

    if xmax <= xmin or zmax <= zmin:
        raise ValueError(
            "Invalid plotting range:\n"
            f"x={xmin} -- {xmax}\n"
            f"z={zmin} -- {zmax}"
        )

    # --------------------------------------------------------
    # 表示領域内のデータを抽出
    # --------------------------------------------------------

    margin_x = 0.05 * (xmax - xmin)
    margin_z = 0.05 * (zmax - zmin)

    source_mask = (
        (x_au >= xmin - margin_x)
        & (x_au <= xmax + margin_x)
        & (z_au >= zmin - margin_z)
        & (z_au <= zmax + margin_z)
        & np.isfinite(rho_cgs)
        & (rho_cgs > 0.0)
        & np.isfinite(Bx_microgauss)
        & np.isfinite(Bz_microgauss)
    )

    number_of_source_points = int(
        np.count_nonzero(source_mask)
    )

    if number_of_source_points < 10:
        raise RuntimeError(
            "Too few source cells in the requested range.\n"
            f"step={step:05d}\n"
            f"source cells={number_of_source_points}\n"
            f"x range={xmin:.3e} -- {xmax:.3e} AU\n"
            f"z range={zmin:.3e} -- {zmax:.3e} AU"
        )

    source_points = np.column_stack([
        x_au[source_mask],
        z_au[source_mask]
    ])

    # --------------------------------------------------------
    # 規則格子を作成
    # --------------------------------------------------------

    x_grid = np.linspace(
        xmin,
        xmax,
        plot_resolution
    )

    z_grid = np.linspace(
        zmin,
        zmax,
        plot_resolution
    )

    X, Z = np.meshgrid(
        x_grid,
        z_grid
    )

    # --------------------------------------------------------
    # 密度をlog空間で補間
    # --------------------------------------------------------

    log_rho_source = np.log10(
        rho_cgs[source_mask]
    )

    log_rho_linear = griddata(
        source_points,
        log_rho_source,
        (X, Z),
        method="linear"
    )

    log_rho_nearest = griddata(
        source_points,
        log_rho_source,
        (X, Z),
        method="nearest"
    )

    log_rho_grid = np.where(
        np.isfinite(log_rho_linear),
        log_rho_linear,
        log_rho_nearest
    )

    rho_grid = 10.0**log_rho_grid

    # --------------------------------------------------------
    # Bx・Bzを補間
    # --------------------------------------------------------

    def interpolate_field(values):
        linear = griddata(
            source_points,
            values[source_mask],
            (X, Z),
            method="linear"
        )

        nearest = griddata(
            source_points,
            values[source_mask],
            (X, Z),
            method="nearest"
        )

        result = np.where(
            np.isfinite(linear),
            linear,
            nearest
        )

        return np.nan_to_num(
            result,
            nan=0.0,
            posinf=0.0,
            neginf=0.0
        )

    Bx_grid = interpolate_field(
        Bx_microgauss
    )

    Bz_grid = interpolate_field(
        Bz_microgauss
    )

    Bproj_grid = np.hypot(
        Bx_grid,
        Bz_grid
    )

    # シンク球と断面の交差領域は、密度だけを白抜きにする。
    # 磁力線はシンク内部も含めて連続的に描画する。
    sink_mask_grid = np.hypot(X, Z) < sink_radius_au
    rho_grid_plot = np.ma.masked_where(sink_mask_grid, rho_grid)

    has_projected_field = np.any(
        Bproj_grid > 0.0
    )

    density_cmap = plt.get_cmap("YlGnBu").copy()
    density_cmap.set_bad("white")

    # --------------------------------------------------------
    # 描画
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=figsize
    )

    density_map = ax.pcolormesh(
        X,
        Z,
        rho_grid_plot,
        cmap=density_cmap,
        norm=LogNorm(
            vmin=rho_vmin,
            vmax=rho_vmax
        ),
        shading="auto",
        rasterized=True
    )

    if has_projected_field:
        ax.streamplot(
            x_grid,
            z_grid,
            Bx_grid,
            Bz_grid,
            color="black",
            density=stream_density,
            linewidth=0.65,
            arrowsize=0.8,
            arrowstyle="->",
            minlength=0.02,
            maxlength=20.0,
            integration_direction="both",
            broken_streamlines=False
        )
    else:
        ax.text(
            0.5,
            0.03,
            r"$B_x=B_z=0$: no in-plane field lines",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
            color="black",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.75,
                "edgecolor": "gray",
            }
        )

    colorbar = fig.colorbar(
        density_map,
        ax=ax,
        pad=0.02
    )

    colorbar.set_label(
        r"Density [g cm$^{-3}$]"
    )

    # --------------------------------------------------------
    # 座標軸を10^n AU形式で表示
    # --------------------------------------------------------

    maximum_absolute_coordinate = max(
        abs(xmin),
        abs(xmax),
        abs(zmin),
        abs(zmax)
    )

    if maximum_absolute_coordinate > 0.0:
        scale_exponent = int(
            np.floor(
                np.log10(
                    maximum_absolute_coordinate
                )
            )
        )
    else:
        scale_exponent = 0

    axis_scale_au = 10.0**scale_exponent

    axis_formatter = FuncFormatter(
        lambda value, position: (
            f"{value / axis_scale_au:g}"
        )
    )

    ax.xaxis.set_major_formatter(
        axis_formatter
    )

    ax.yaxis.set_major_formatter(
        axis_formatter
    )

    ax.set_xlabel(
        rf"$x\ [10^{{{scale_exponent}}}\ "
        rf"{{\rm AU}}]$"
    )

    ax.set_ylabel(
        rf"$z\ [10^{{{scale_exponent}}}\ "
        rf"{{\rm AU}}]$"
    )

    zoom_title = (
        ""
        if zoom_fraction is None
        else ", 1/10 zoom"
    )

    ax.set_title(
        "Density and projected magnetic field lines: x-z\n"
        f"step={step:05d}, "
        f"t={time_kyr:.1f} kyr"
        f"{zoom_title}"
    )

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(zmin, zmax)

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.grid(
        True,
        alpha=0.15,
        linestyle="--"
    )

    # --------------------------------------------------------
    # 表示範囲内の磁場強度
    # --------------------------------------------------------

    Bmax_visible = np.max(
        Bmag_microgauss[source_mask]
    )

    Bxmax_visible = np.max(
        np.abs(
            Bx_microgauss[source_mask]
        )
    )

    Bzmax_visible = np.max(
        np.abs(
            Bz_microgauss[source_mask]
        )
    )

    information = (
        f"|B|max  = {Bmax_visible:.3e} μG\n"
        f"|Bx|max = {Bxmax_visible:.3e} μG\n"
        f"|Bz|max = {Bzmax_visible:.3e} μG"
    )

    ax.text(
        0.02,
        0.98,
        information,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox={
            "boxstyle": "round",
            "facecolor": "black",
            "alpha": 0.55,
            "edgecolor": "white",
        }
    )

    plt.tight_layout()

    suffix = (
        f"_{filename_suffix}"
        if filename_suffix
        else ""
    )

    output_file = os.path.join(
        output_dir,
        (
            f"xz_density_Blines_"
            f"{step:05d}"
            f"{suffix}.png"
        )
    )

    fig.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(fig)

    # メモリ解放
    del (
        data,
        x_code,
        z_code,
        rho_code,
        B_code,
        x_au,
        z_au,
        rho_cgs,
        Bx_microgauss,
        By_microgauss,
        Bz_microgauss,
        Bmag_microgauss,
        source_points,
        X,
        Z,
        rho_grid,
        Bx_grid,
        Bz_grid
    )

    gc.collect()

    print(
        f"[SAVED] step={step:05d}, "
        f"time={time_kyr:.1f} kyr\n"
        f"        {output_file}"
    )

    return output_file


# ============================================================
# 10. 全体図を保存
# ============================================================

print("\n=== Generating full-domain images ===")

full_output_files = []

for index, step in enumerate(timesteps):
    print(
        f"[FULL {index + 1:3d}/{len(timesteps):3d}] "
        f"Processing step={step:05d}"
    )

    output_file = plot_xz_density_with_fieldlines(
        step=step,
        output_dir=output_dir_full,
        zoom_fraction=None,
        filename_suffix="full"
    )

    full_output_files.append(
        output_file
    )


# ============================================================
# 11. 1/10ズーム図を保存
# ============================================================

print("\n=== Generating 1/10 zoom images ===")

zoom10_output_files = []

for index, step in enumerate(timesteps):
    print(
        f"[ZOOM 1/10 "
        f"{index + 1:3d}/{len(timesteps):3d}] "
        f"Processing step={step:05d}"
    )

    output_file = plot_xz_density_with_fieldlines(
        step=step,
        output_dir=output_dir_zoom10,
        zoom_fraction=zoom_fraction_10,
        filename_suffix="zoom_1over10"
    )

    zoom10_output_files.append(
        output_file
    )


# ============================================================
# 12. 保存結果
# ============================================================

print("\n=== Finished ===")
print(f"Full images : {len(full_output_files)}")
print(f"1/10 zoom   : {len(zoom10_output_files)}")

print(
    "\nFull output directory:\n"
    f"{output_dir_full}"
)

print(
    "\n1/10 zoom output directory:\n"
    f"{output_dir_zoom10}"
)

print("\nFull images:")

for filename in full_output_files:
    print(
        "  ",
        os.path.basename(filename)
    )

print("\n1/10 zoom images:")

for filename in zoom10_output_files:
    print(
        "  ",
        os.path.basename(filename)
    )

# %% cell 2
#x-y磁場のストリームラインon密度map

import os
import re
import glob
from collections import defaultdict

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter


# ============================================================
# 1. ユーザー設定
# ============================================================

# VTKファイルがあるディレクトリ
vtk_dir = os.path.expanduser(
    "~/athena-project/results/〇〇"
)

# 親出力ディレクトリ
output_root_xy = os.path.join(
    vtk_dir,
    "xy_density_with_magnetic_fieldlines"
)

# x-z版と同じく、全体図とズーム図を同じ親ディレクトリに格納
output_dir1 = os.path.join(output_root_xy, "full")
output_dir2 = os.path.join(output_root_xy, "zoom_1over10")

os.makedirs(output_dir1, exist_ok=True)
os.makedirs(output_dir2, exist_ok=True)

# zoomの各軸幅を全体領域の1/10にする
zoom_fraction_10 = 1.0 / 10.0

# 補間後の画像解像度
plot_resolution = 700

# 磁力線の密度
stream_density = 1.5

# x-y面の法線（z方向）に平均する共通の半厚み [AU]。
# AMRレベルによらず、全MeshBlockで同じ物理厚みを使う。
slice_half_thickness_au = 1000.0

# 白抜きにする球対称シンク領域の半径 [AU]
sink_radius_au = 1000.0

# 図のサイズ
figsize = (9, 8)

# 保存画像のdpi
dpi = 200

print("=== Output directories ===")
print(f"Output root : {output_root_xy}")
print(f"Full images : {output_dir1}")
print(f"1/10 zoom   : {output_dir2}")


# ============================================================
# 2. Toyouchi.cppのコード単位
# ============================================================

Munit = 4.0e33       # g
Lunit = 7.03e15      # cm (L0 = rb / 2)
Tunit = 3.61e10      # s

Vunit = Lunit / Tunit
Rhounit = Munit / Lunit**3
Punit = Rhounit * Vunit**2

# Athena++の磁場単位 [Gauss]
Bunit = np.sqrt(4.0 * np.pi * Punit)

AU = 1.496e13
YEAR = 3.15576e7

print("\n=== Code units ===")
print(f"Lunit   = {Lunit:.6e} cm")
print(f"Rhounit = {Rhounit:.6e} g cm^-3")
print(f"Bunit   = {Bunit:.6e} Gauss")
print(f"Tunit   = {Tunit/(1.0e6*YEAR):.6e} Myr")


# ============================================================
# 3. VTKファイルを出力時刻ごとに整理
# ============================================================

timestep_dict = defaultdict(list)

vtk_pattern = os.path.join(
    vtk_dir,
    "Toyouchi.block*.out2.*.vtk"
)

vtk_files = glob.glob(vtk_pattern)

if not vtk_files:
    raise FileNotFoundError(
        "VTK files were not found.\n"
        f"Checked pattern:\n{vtk_pattern}"
    )

for filename in vtk_files:
    match = re.search(
        r"\.out2\.(\d+)\.vtk$",
        filename
    )

    if match:
        step = int(match.group(1))
        timestep_dict[step].append(filename)

timesteps = sorted(timestep_dict)

if not timesteps:
    raise RuntimeError(
        "Could not extract timestep numbers "
        "from VTK filenames."
    )

# block番号順に並べる
for step in timesteps:
    timestep_dict[step] = sorted(
        timestep_dict[step],
        key=lambda filename: int(
            re.search(
                r"block(\d+)",
                filename
            ).group(1)
        )
    )

print("\n=== VTK information ===")
print(f"Number of VTK files     : {len(vtk_files)}")
print(f"Number of output steps  : {len(timesteps)}")
print(f"First step              : {timesteps[0]:05d}")
print(f"Last step               : {timesteps[-1]:05d}")
print(
    f"Blocks at first step    : "
    f"{len(timestep_dict[timesteps[0]])}"
)
print(
    f"Blocks at last step     : "
    f"{len(timestep_dict[timesteps[-1]])}"
)


# ============================================================
# 4. Athena++ VTKヘッダーから時刻を読む
# ============================================================

def read_athena_vtk_time(filename):
    """
    Athena++ legacy VTKの2行目からコード時刻を取得する。
    """

    with open(filename, "rb") as file:
        file.readline()
        title = file.readline().decode(
            "ascii",
            errors="ignore"
        )

    match = re.search(
        r"time=([0-9eE+\-.]+)",
        title
    )

    if match is None:
        return np.nan

    return float(match.group(1))


# ============================================================
# 5. 配列名を探す
# ============================================================

def find_array_name(grid, candidates):
    """
    candidatesの中から、VTKに存在する配列名を返す。
    """

    for name in candidates:
        if name in grid.array_names:
            return name

    raise KeyError(
        "Required array was not found.\n"
        f"Candidates: {candidates}\n"
        f"Available arrays: {grid.array_names}"
    )


# ============================================================
# 6. x-y断面を読み込む
# ============================================================

def load_xy_slice(step):
    """
    z=0まわりの固定物理厚みのスラブと交差するセルを抽出する。
    後段では各セルとスラブが重なるz方向の長さを重みとして、
    同一(x,y)座標の密度と磁場を平均する。
    """

    x_list = []
    y_list = []
    rho_list = []
    B_list = []
    weight_list = []

    # z=0と交差する全MeshBlockの実境界
    domain_xmin = np.inf
    domain_xmax = -np.inf
    domain_ymin = np.inf
    domain_ymax = -np.inf

    files = timestep_dict[step]

    time_code = read_athena_vtk_time(
        files[0]
    )

    slab_half_thickness_code = (
        slice_half_thickness_au * AU / Lunit
    )

    for filename in files:
        grid = pv.read(filename)

        # PyVista bounds:
        # xmin, xmax, ymin, ymax, zmin, zmax
        bounds = grid.bounds

        zmin = bounds[4]
        zmax = bounds[5]

        # 固定厚みのスラブと交差しないMeshBlockを除外
        if (
            zmax < -slab_half_thickness_code
            or zmin > slab_half_thickness_code
        ):
            del grid
            continue

        domain_xmin = min(domain_xmin, bounds[0])
        domain_xmax = max(domain_xmax, bounds[1])
        domain_ymin = min(domain_ymin, bounds[2])
        domain_ymax = max(domain_ymax, bounds[3])

        points = grid.cell_centers().points

        rho_name = find_array_name(
            grid,
            [
                "rho",
                "dens",
                "density",
                "prim_dens",
                "prim_density",
            ]
        )

        B_name = find_array_name(
            grid,
            [
                "Bcc",
                "bcc",
                "B",
                "magnetic_field",
            ]
        )

        rho = np.asarray(
            grid[rho_name]
        ).reshape(-1)

        B = np.asarray(
            grid[B_name]
        ).reshape(-1, 3)

        z = points[:, 2]

        unique_z = np.unique(z)
        dz_values = np.diff(unique_z)
        dz_values = dz_values[dz_values > 0.0]
        if dz_values.size:
            dz_local = np.min(dz_values)
        else:
            dz_local = max(zmax - zmin, 1.0e-12)

        # 各セル区間と共通スラブとの重なり長さ。
        # 境界セルは重なっている割合だけ平均へ寄与する。
        cell_lower = z - 0.5 * dz_local
        cell_upper = z + 0.5 * dz_local
        overlap = np.minimum(
            cell_upper,
            slab_half_thickness_code,
        ) - np.maximum(
            cell_lower,
            -slab_half_thickness_code,
        )
        overlap = np.clip(overlap, 0.0, None)
        mask = overlap > 0.0

        if not np.any(mask):
            continue

        x_list.append(
            points[mask, 0]
        )

        y_list.append(
            points[mask, 1]
        )

        rho_list.append(
            rho[mask]
        )

        B_list.append(
            B[mask]
        )

        weight_list.append(
            overlap[mask]
        )

        del points, rho, B, grid

    if not x_list:
        raise RuntimeError(
            "No cells in the finite-thickness slab around z=0 were found.\n"
            f"step = {step:05d}"
        )

    return {
        "step": step,
        "time_code": time_code,
        "x_code": np.concatenate(x_list),
        "y_code": np.concatenate(y_list),
        "rho_code": np.concatenate(rho_list),
        "B_code": np.vstack(B_list),
        "slab_weight_code": np.concatenate(weight_list),
        "domain_xmin_code": domain_xmin,
        "domain_xmax_code": domain_xmax,
        "domain_ymin_code": domain_ymin,
        "domain_ymax_code": domain_ymax,
    }


# ============================================================
# 7. 同じx-y座標の値を平均する
# ============================================================

def average_duplicate_xy(
    x,
    y,
    rho,
    B,
    weights,
):
    """
    同じ(x,y)座標に複数のセル値がある場合に平均する。

    固定厚みスラブとの重なり長さで重み付けし、
    z=0まわりの複数層やMeshBlock境界の重複を処理する。
    """

    coords = np.column_stack([
        x,
        y
    ])

    unique_coords, inverse = np.unique(
        coords,
        axis=0,
        return_inverse=True
    )

    weight_sum = np.bincount(inverse, weights=weights)

    rho_sum = np.bincount(
        inverse,
        weights=rho * weights
    )

    Bx_sum = np.bincount(
        inverse,
        weights=B[:, 0] * weights
    )

    By_sum = np.bincount(
        inverse,
        weights=B[:, 1] * weights
    )

    Bnormal_sum = np.bincount(
        inverse,
        weights=B[:, 2] * weights
    )

    rho_avg = rho_sum / weight_sum

    B_avg = np.column_stack([
        Bx_sum / weight_sum,
        By_sum / weight_sum,
        Bnormal_sum / weight_sum
    ])

    return (
        unique_coords[:, 0],
        unique_coords[:, 1],
        rho_avg,
        B_avg
    )


# ============================================================
# 8. 全出力で共通の密度カラースケールを決める
# ============================================================

print(
    "\n[INFO] Determining common "
    "density color scale..."
)

density_samples = []

for n, step in enumerate(timesteps):
    slice_data = load_xy_slice(step)

    _, _, rho_average_code, _ = average_duplicate_xy(
        slice_data["x_code"],
        slice_data["y_code"],
        slice_data["rho_code"],
        slice_data["B_code"],
        slice_data["slab_weight_code"],
    )

    rho_cgs = (
        rho_average_code * Rhounit
    )

    valid = (
        np.isfinite(rho_cgs)
        & (rho_cgs > 0.0)
    )

    if np.any(valid):
        density_samples.append(
            rho_cgs[valid]
        )

    print(
        f"  [{n+1:3d}/{len(timesteps):3d}] "
        f"step={step:05d}"
    )

if not density_samples:
    raise RuntimeError(
        "No positive finite density values "
        "were found."
    )

all_slice_density = np.concatenate(
    density_samples
)

rho_vmin = np.percentile(
    all_slice_density,
    0.5
)

# 全時刻の最大密度を含める
rho_vmax = np.max(
    all_slice_density
)

if not np.isfinite(rho_vmin):
    raise RuntimeError(
        "rho_vmin is not finite."
    )

if not np.isfinite(rho_vmax):
    raise RuntimeError(
        "rho_vmax is not finite."
    )

if rho_vmin <= 0.0:
    raise RuntimeError(
        "rho_vmin must be positive "
        "for LogNorm."
    )

if rho_vmax <= rho_vmin:
    raise RuntimeError(
        "Invalid density color range:\n"
        f"rho_vmin={rho_vmin}\n"
        f"rho_vmax={rho_vmax}"
    )

print(
    f"\nrho_vmin = {rho_vmin:.6e} g cm^-3"
)

print(
    f"rho_vmax = {rho_vmax:.6e} g cm^-3"
)


# ============================================================
# 9. 1時刻分を描画・保存する
# ============================================================

def plot_xy_density_with_fieldlines(
    step,
    output_dir,
    zoom_fraction=None,
    filename_suffix=""
):
    """
    x-y密度マップへ、x-y面に射影した磁力線を重ねる。

    Parameters
    ----------
    step : int
        out2番号

    output_dir : str
        保存先

    zoom_fraction : float or None
        Noneなら実データ全体。0.1なら各軸幅を全体の1/10にする。

    filename_suffix : str
        出力ファイル名に付ける文字列
    """

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    data = load_xy_slice(step)

    x_code = data["x_code"]
    y_code = data["y_code"]
    rho_code = data["rho_code"]
    B_code = data["B_code"]
    slab_weight_code = data["slab_weight_code"]

    # 重複するx-y座標を平均
    (
        x_code,
        y_code,
        rho_code,
        B_code
    ) = average_duplicate_xy(
        x_code,
        y_code,
        rho_code,
        B_code,
        slab_weight_code,
    )

    # --------------------------------------------------------
    # コード単位からCGS・表示単位へ変換
    # --------------------------------------------------------

    x_AU = (
        x_code
        * Lunit
        / AU
    )

    y_AU = (
        y_code
        * Lunit
        / AU
    )

    rho_cgs = (
        rho_code
        * Rhounit
    )

    Bx_microG = (
        B_code[:, 0]
        * Bunit
        * 1.0e6
    )

    By_microG = (
        B_code[:, 1]
        * Bunit
        * 1.0e6
    )

    Bnormal_microG = (
        B_code[:, 2]
        * Bunit
        * 1.0e6
    )

    Bmag_microG = np.sqrt(
        Bx_microG**2
        + Bnormal_microG**2
        + By_microG**2
    )

    time_code = data["time_code"]

    time_kyr = (
        time_code
        * Tunit
        / (1.0e3 * YEAR)
    )

    # --------------------------------------------------------
    # 描画範囲
    # --------------------------------------------------------

    # セル中心の最小・最大ではなく、VTKブロック境界から実領域を取得
    data_xmin = data["domain_xmin_code"] * Lunit / AU
    data_xmax = data["domain_xmax_code"] * Lunit / AU
    data_ymin = data["domain_ymin_code"] * Lunit / AU
    data_ymax = data["domain_ymax_code"] * Lunit / AU

    if zoom_fraction is None:
        xmin, xmax = data_xmin, data_xmax
        ymin, ymax = data_ymin, data_ymax
        range_label = "full domain"
    else:
        if not (0.0 < zoom_fraction <= 1.0):
            raise ValueError(
                "zoom_fraction must satisfy 0 < zoom_fraction <= 1"
            )

        x_center = 0.0 if data_xmin <= 0.0 <= data_xmax else 0.5 * (data_xmin + data_xmax)
        y_center = 0.0 if data_ymin <= 0.0 <= data_ymax else 0.5 * (data_ymin + data_ymax)

        x_half = 0.5 * (data_xmax - data_xmin) * zoom_fraction
        y_half = 0.5 * (data_ymax - data_ymin) * zoom_fraction

        xmin = max(data_xmin, x_center - x_half)
        xmax = min(data_xmax, x_center + x_half)
        ymin = max(data_ymin, y_center - y_half)
        ymax = min(data_ymax, y_center + y_half)
        range_label = f"zoom = {zoom_fraction:.3g} of full width"

    print(
        f"[RANGE] step={step:05d}, {range_label}\n"
        f"        full x={data_xmin:.6e} -- {data_xmax:.6e} AU\n"
        f"        full y={data_ymin:.6e} -- {data_ymax:.6e} AU\n"
        f"        plot x={xmin:.6e} -- {xmax:.6e} AU\n"
        f"        plot y={ymin:.6e} -- {ymax:.6e} AU"
    )

    if xmax <= xmin:
        raise ValueError(
            f"Invalid x range: {xmin}, {xmax}"
        )

    if ymax <= ymin:
        raise ValueError(
            f"Invalid y range: {ymin}, {ymax}"
        )

    # 表示範囲の外側を少し含めて補間する
    margin_x = 0.05 * (
        xmax - xmin
    )

    margin_z = 0.05 * (
        ymax - ymin
    )

    mask = (
        (x_AU >= xmin - margin_x)
        & (x_AU <= xmax + margin_x)
        & (y_AU >= ymin - margin_z)
        & (y_AU <= ymax + margin_z)
        & np.isfinite(rho_cgs)
        & (rho_cgs > 0.0)
        & np.isfinite(Bx_microG)
        & np.isfinite(By_microG)
    )

    number_of_source_points = (
        np.count_nonzero(mask)
    )

    if number_of_source_points < 10:
        raise RuntimeError(
            "Too few source cells in the "
            "requested plotting range.\n"
            f"step={step:05d}\n"
            f"source cells={number_of_source_points}\n"
            f"x range={xmin:.3e} to {xmax:.3e} AU\n"
            f"y range={ymin:.3e} to {ymax:.3e} AU"
        )

    # --------------------------------------------------------
    # 規則格子
    # --------------------------------------------------------

    x_grid = np.linspace(
        xmin,
        xmax,
        plot_resolution
    )

    y_grid = np.linspace(
        ymin,
        ymax,
        plot_resolution
    )

    X, Y = np.meshgrid(
        x_grid,
        y_grid
    )

    source_points = np.column_stack([
        x_AU[mask],
        y_AU[mask]
    ])

    # --------------------------------------------------------
    # 密度をlog空間で補間
    # --------------------------------------------------------

    log_rho_source = np.log10(
        rho_cgs[mask]
    )

    log_rho_linear = griddata(
        source_points,
        log_rho_source,
        (X, Y),
        method="linear"
    )

    log_rho_nearest = griddata(
        source_points,
        log_rho_source,
        (X, Y),
        method="nearest"
    )

    log_rho_grid = np.where(
        np.isfinite(log_rho_linear),
        log_rho_linear,
        log_rho_nearest
    )

    rho_grid = (
        10.0**log_rho_grid
    )

    # --------------------------------------------------------
    # Bx, Byを補間
    # --------------------------------------------------------

    Bx_linear = griddata(
        source_points,
        Bx_microG[mask],
        (X, Y),
        method="linear"
    )

    By_linear = griddata(
        source_points,
        By_microG[mask],
        (X, Y),
        method="linear"
    )

    Bx_nearest = griddata(
        source_points,
        Bx_microG[mask],
        (X, Y),
        method="nearest"
    )

    By_nearest = griddata(
        source_points,
        By_microG[mask],
        (X, Y),
        method="nearest"
    )

    Bx_grid = np.where(
        np.isfinite(Bx_linear),
        Bx_linear,
        Bx_nearest
    )

    By_grid = np.where(
        np.isfinite(By_linear),
        By_linear,
        By_nearest
    )

    Bx_grid = np.nan_to_num(
        Bx_grid,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    By_grid = np.nan_to_num(
        By_grid,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    # x-y面内の磁場強度。弱磁場を除外する閾値は設けない。
    # 方向だけをstreamplotへ渡すことで、絶対強度が非常に小さくても
    # 同じ磁力線形状を追跡できるようにする。
    Bproj_grid = np.hypot(Bx_grid, By_grid)
    nonzero_projected_field = Bproj_grid > 0.0

    Bx_direction = np.divide(
        Bx_grid,
        Bproj_grid,
        out=np.zeros_like(Bx_grid),
        where=nonzero_projected_field,
    )
    By_direction = np.divide(
        By_grid,
        Bproj_grid,
        out=np.zeros_like(By_grid),
        where=nonzero_projected_field,
    )

    # シンク球と断面の交差領域は、密度だけを白抜きにする。
    # 磁力線はシンク内部も含めて連続的に描画する。
    sink_mask_grid = np.hypot(X, Y) < sink_radius_au
    rho_grid_plot = np.ma.masked_where(sink_mask_grid, rho_grid)

    # 全領域で完全にBx=By=0の場合だけ磁力線を描けない。
    has_projected_field = np.any(nonzero_projected_field)

    density_cmap = plt.get_cmap("YlGnBu").copy()
    density_cmap.set_bad("white")

    # --------------------------------------------------------
    # 描画
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=figsize
    )

    density_map = ax.pcolormesh(
        X,
        Y,
        rho_grid_plot,
        # 論文図風：低密度=淡い黄緑、高密度=濃い青
        cmap=density_cmap,
        norm=LogNorm(
            vmin=rho_vmin,
            vmax=rho_vmax
        ),
        shading="auto",
        rasterized=True
    )

    if has_projected_field:
        ax.streamplot(
            x_grid,
            y_grid,
            Bx_direction,
            By_direction,
            color="black",
            density=stream_density,
            linewidth=0.65,
            arrowsize=0.8,
            arrowstyle="->",
            minlength=0.02,
            maxlength=20.0,
            integration_direction="both",
            # 他の磁力線に近づいても積分を打ち切らない
            broken_streamlines=False
        )
    else:
        ax.text(
            0.50,
            0.03,
            r"$B_x=B_y=0$: no in-plane field lines",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
            color="black",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.75,
                "edgecolor": "gray",
            },
        )

    colorbar = fig.colorbar(
        density_map,
        ax=ax,
        pad=0.02
    )

    colorbar.set_label(
        r"Density [g cm$^{-3}$]"
    )

    # --------------------------------------------------------
    # 表示範囲に応じて、目盛を10^n AU単位で表示
    # 例：表示範囲が10^5 AU級なら x [10^5 AU]
    # --------------------------------------------------------
    maximum_absolute_coordinate = max(
        abs(xmin),
        abs(xmax),
        abs(ymin),
        abs(ymax),
    )

    if maximum_absolute_coordinate > 0.0:
        scale_exponent = int(
            np.floor(
                np.log10(maximum_absolute_coordinate)
            )
        )
    else:
        scale_exponent = 0

    axis_scale_AU = 10.0**scale_exponent

    axis_formatter = FuncFormatter(
        lambda value, position: (
            f"{value / axis_scale_AU:g}"
        )
    )

    ax.xaxis.set_major_formatter(
        axis_formatter
    )
    ax.yaxis.set_major_formatter(
        axis_formatter
    )

    ax.set_xlabel(
        rf"$x\ [10^{{{scale_exponent}}}\ {{\rm AU}}]$"
    )
    ax.set_ylabel(
        rf"$y\ [10^{{{scale_exponent}}}\ {{\rm AU}}]$"
    )

    ax.set_title(
        "Density and projected magnetic field lines: x-y midplane\n"
        f"step={step:05d}, "
        f"t={time_kyr:.1f} kyr"
    )

    ax.set_xlim(
        xmin,
        xmax
    )

    ax.set_ylim(
        ymin,
        ymax
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.grid(
        True,
        alpha=0.15,
        linestyle="--"
    )

    # 表示範囲内の磁場最大値
    Bmax_visible = np.max(
        Bmag_microG[mask]
    )

    Bxmax_visible = np.max(
        np.abs(Bx_microG[mask])
    )

    Bymax_visible = np.max(
        np.abs(By_microG[mask])
    )

    # 情報ボックスは通常文字列にして、環境ごとのmathtext差を避ける。
    info = (
        f"Bmax = {Bmax_visible:.3e} μG\n"
        f"|Bx|max = {Bxmax_visible:.3e} μG\n"
        f"|By|max = {Bymax_visible:.3e} μG"
    )

    ax.text(
        0.02,
        0.98,
        info,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox={
            "boxstyle": "round",
            "facecolor": "black",
            "alpha": 0.55,
            "edgecolor": "white",
        }
    )

    plt.tight_layout()

    if filename_suffix:
        suffix = (
            "_" + filename_suffix
        )
    else:
        suffix = ""

    output_file = os.path.join(
        output_dir,
        (
            f"xy_density_Blines_"
            f"{step:05d}"
            f"{suffix}.png"
        )
    )

    fig.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white"
    )

    # メモリ解放
    plt.close(fig)

    print(
        f"[SAVED] step={step:05d}, "
        f"time={time_kyr:.1f} kyr\n"
        f"        {output_file}"
    )

    return output_file


# ============================================================
# 10. 全体図をoutput_dir1へ保存
# ============================================================

print("\n=== Generating full-domain images ===")

full_output_files = []

for n, step in enumerate(timesteps):
    print(
        f"[FULL {n+1:3d}/{len(timesteps):3d}] "
        f"Processing step={step:05d}"
    )

    output_file = (
        plot_xy_density_with_fieldlines(
            step=step,
            output_dir=output_dir1,
            zoom_fraction=None,
            filename_suffix="full"
        )
    )

    full_output_files.append(
        output_file
    )


# ============================================================
# 11. 1/10ズーム図をoutput_dir2へ保存
# ============================================================

print("\n=== Generating 1/10 zoom images ===")

zoom10_output_files = []

for n, step in enumerate(timesteps):
    print(
        f"[ZOOM 1/10 {n+1:3d}/{len(timesteps):3d}] "
        f"Processing step={step:05d}"
    )

    output_file = (
        plot_xy_density_with_fieldlines(
            step=step,
            output_dir=output_dir2,
            zoom_fraction=zoom_fraction_10,
            filename_suffix="zoom_1over10"
        )
    )

    zoom10_output_files.append(
        output_file
    )


# ============================================================
# 12. 保存結果を確認
# ============================================================

print("\n=== Finished ===")
print(
    f"Full images : "
    f"{len(full_output_files)}"
)

print(
    f"1/10 zoom   : "
    f"{len(zoom10_output_files)}"
)

print(f"\nFull output directory:\n{output_dir1}")
print(f"\n1/10 zoom output directory:\n{output_dir2}")

print("\nFull images:")
for filename in full_output_files:
    print(
        "  ",
        os.path.basename(filename)
    )

print("\n1/10 zoom images:")
for filename in zoom10_output_files:
    print(
        "  ",
        os.path.basename(filename)
    )


     

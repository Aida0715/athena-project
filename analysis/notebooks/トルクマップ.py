# -*- coding: utf-8 -*-
# Converted from トルクマップ.ipynb

# %% cell 1
# ============================================================
# 全タイムステップ：
# x-y平面 密度・速度方向＋角運動量輸送診断
#
# パネル
#   1. ガス自己重力トルク密度
#   2. 圧力勾配トルク密度
#   3. 移流角運動量流束
#
# 表示
#   背景：密度（白黒）
#   赤：正
#   青：負
#   黒矢印：
#     自己重力・圧力パネルでは規格化した(vx, vy)方向
#     移流フラックスパネルでは規格化した動径速度方向
#       sign(v_R) e_R（外向き流は外向き、降着流は中心向き）
#
# 単位
#   距離：10^4 AU
#   時間：Myr
#   密度：M_sun AU^-3
#   トルク診断量：code units
# ============================================================

import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import griddata


# ============================================================
# 入出力
# ============================================================
vtk_dir = Path(
    os.path.expanduser(
        "~/athena-project/results/〇〇"
    )
).resolve()

output_dir = (
    vtk_dir
    / "xy_torque_maps"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

STREAM_TOKEN = ".out2."

print(f"[INFO] Input : {vtk_dir}")
print(f"[INFO] Output: {output_dir}")

if not vtk_dir.is_dir():
    raise FileNotFoundError(
        f"VTK directory does not exist: {vtk_dir}"
    )


# ============================================================
# コード単位
# ============================================================
M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15  # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33

MYR_CGS = (
    1.0e6
    * 365.25
    * 24.0
    * 3600.0
)

LENGTH_UNIT_AU = (
    L_UNIT_CGS / AU_CGS
)

TIME_UNIT_MYR = (
    T_UNIT_CGS / MYR_CGS
)

MASS_UNIT_MSUN = (
    M_UNIT_CGS / MSUN_CGS
)

DENSITY_UNIT_MSUN_AU3 = (
    MASS_UNIT_MSUN
    / LENGTH_UNIT_AU**3
)

print(
    f"[INFO] 1 code length  = "
    f"{LENGTH_UNIT_AU:.6e} AU"
)
print(
    f"[INFO] 1 code time    = "
    f"{TIME_UNIT_MYR:.6e} Myr"
)
print(
    f"[INFO] 1 code density = "
    f"{DENSITY_UNIT_MSUN_AU3:.6e} "
    "M_sun AU^-3"
)


# ============================================================
# 物理・描画設定
# ============================================================

# 表示範囲：-10^5 ～ +10^5 AU
PLOT_RADIUS_AU = 1.0e5

# 等温音速（code velocity）
ISO_SOUND_SPEED_CODE = 0.707  # L0 = rb / 2

# 本描画解像度
MAP_RESOLUTION = 700

# 全時刻のトルク範囲を調べる低解像度
SCALE_SCAN_RESOLUTION = 250

# 矢印数
QUIVER_N = 31

# 密度の全時刻共通範囲
DENSITY_PERCENTILES = (
    5.0,
    99.0,
)

# トルク・流束の全時刻共通範囲
TORQUE_ABS_PERCENTILE = 99.0

# シンク内部を統計から除外
SINK_RADIUS_AU = 1000.0

# 極端に遅い領域の矢印を非表示
MIN_SPEED_FRACTION = 0.01

# トルクの重ね描き透明度
TORQUE_ALPHA = 0.68

DENSITY_CMAP = "Greys"
TORQUE_CMAP = "RdBu_r"

FIGSIZE = (18, 6)
DPI = 200


# ============================================================
# VTK変数名候補
# ============================================================
DENSITY_CANDIDATES = (
    "rho",
    "dens",
    "density",
)

VELOCITY_CANDIDATES = (
    "vel",
    "velocity",
    "v",
)

MOMENTUM_CANDIDATES = (
    "mom",
    "momentum",
)

POTENTIAL_CANDIDATES = (
    "phi",
    "potential",
    "grav_potential",
)

# 将来の磁場出力に対応する候補
MAGNETIC_VECTOR_CANDIDATES = (
    "Bcc",
    "Bcc_xyz",
    "bcc",
    "magnetic_field",
    "B",
)

MAGNETIC_COMPONENT_CANDIDATES = (
    ("Bcc1", "Bcc2", "Bcc3"),
    ("bcc1", "bcc2", "bcc3"),
    ("bx", "by", "bz"),
    ("Bx", "By", "Bz"),
)


# ============================================================
# VTKヘッダから時刻取得
# ============================================================
def read_vtk_time_code(filename):
    with open(filename, "rb") as vtk_file:
        header = vtk_file.read(2048).decode(
            "ascii",
            errors="ignore",
        )

    match = re.search(
        r"time\s*=\s*"
        r"([+-]?(?:\d+\.?\d*|\.\d+)"
        r"(?:[eE][+-]?\d+)?)",
        header,
        flags=re.IGNORECASE,
    )

    if match is None:
        raise ValueError(
            f"VTK time not found: {filename}"
        )

    return float(match.group(1))


def get_step_number(filename):
    match = re.search(
        r"(?:prim\.)?out2\.(\d+)",
        Path(filename).name,
    )

    if match is None:
        return None

    return int(match.group(1))


# ============================================================
# フィールド検索
# ============================================================
def find_cell_field(grid, candidates):
    lower_to_original = {
        name.lower(): name
        for name in grid.cell_data.keys()
    }

    for candidate in candidates:
        key = candidate.lower()

        if key in lower_to_original:
            return lower_to_original[key]

    return None


def extract_base_fields(grid, filename):
    density_name = find_cell_field(
        grid,
        DENSITY_CANDIDATES,
    )

    potential_name = find_cell_field(
        grid,
        POTENTIAL_CANDIDATES,
    )

    if density_name is None:
        raise KeyError(
            f"Density field not found: {filename}"
        )

    if potential_name is None:
        raise KeyError(
            f"Self-gravity potential phi "
            f"not found: {filename}"
        )

    rho = np.asarray(
        grid.cell_data[density_name]
    ).reshape(-1)

    phi = np.asarray(
        grid.cell_data[potential_name]
    ).reshape(-1)

    velocity_name = find_cell_field(
        grid,
        VELOCITY_CANDIDATES,
    )

    if velocity_name is not None:
        velocity = np.asarray(
            grid.cell_data[velocity_name]
        )

    else:
        momentum_name = find_cell_field(
            grid,
            MOMENTUM_CANDIDATES,
        )

        if momentum_name is None:
            raise KeyError(
                f"Velocity/momentum field "
                f"not found: {filename}"
            )

        momentum = np.asarray(
            grid.cell_data[momentum_name]
        )

        velocity = (
            momentum
            / np.maximum(
                rho[:, None],
                1.0e-300,
            )
        )

    if (
        velocity.ndim != 2
        or velocity.shape[1] < 3
    ):
        raise ValueError(
            f"Unexpected velocity shape: "
            f"{velocity.shape}"
        )

    return rho, phi, velocity[:, :3]


# ============================================================
# 将来用：磁場の検索関数
# 現時点ではトルク計算には使用しない
# ============================================================
def extract_optional_magnetic_field(grid):
    vector_name = find_cell_field(
        grid,
        MAGNETIC_VECTOR_CANDIDATES,
    )

    if vector_name is not None:
        magnetic = np.asarray(
            grid.cell_data[vector_name]
        )

        if (
            magnetic.ndim == 2
            and magnetic.shape[1] >= 3
        ):
            return magnetic[:, :3]

    lower_to_original = {
        name.lower(): name
        for name in grid.cell_data.keys()
    }

    for component_names in (
        MAGNETIC_COMPONENT_CANDIDATES
    ):
        if all(
            name.lower() in lower_to_original
            for name in component_names
        ):
            actual_names = [
                lower_to_original[name.lower()]
                for name in component_names
            ]

            return np.column_stack([
                np.asarray(
                    grid.cell_data[name]
                ).reshape(-1)
                for name in actual_names
            ])

    return None


# ============================================================
# 正確なz=0面を抽出
# ============================================================
def extract_xy_midplane(filename):
    grid = pv.read(filename)

    xmin, xmax, ymin, ymax, zmin, zmax = (
        grid.bounds
    )

    tolerance = (
        1.0e-12
        * max(
            abs(zmin),
            abs(zmax),
            1.0,
        )
    )

    if not (
        zmin - tolerance
        <= 0.0
        <= zmax + tolerance
    ):
        return None

    rho_code, phi_code, velocity = (
        extract_base_fields(
            grid,
            filename,
        )
    )

    working_grid = grid.copy()

    working_grid.cell_data[
        "__rho__"
    ] = rho_code

    working_grid.cell_data[
        "__phi__"
    ] = phi_code

    working_grid.cell_data[
        "__velocity__"
    ] = velocity

    # 将来磁場が追加されたときに認識
    magnetic = extract_optional_magnetic_field(
        grid
    )

    if magnetic is not None:
        working_grid.cell_data[
            "__magnetic__"
        ] = magnetic

    point_grid = (
        working_grid.cell_data_to_point_data(
            pass_cell_data=False,
        )
    )

    sliced = point_grid.slice(
        normal=(0.0, 0.0, 1.0),
        origin=(0.0, 0.0, 0.0),
    )

    if sliced.n_points == 0:
        return None

    points_code = np.asarray(
        sliced.points
    )

    rho_slice = np.asarray(
        sliced.point_data["__rho__"]
    ).reshape(-1)

    phi_slice = np.asarray(
        sliced.point_data["__phi__"]
    ).reshape(-1)

    velocity_slice = np.asarray(
        sliced.point_data["__velocity__"]
    )

    x_code = points_code[:, 0]
    y_code = points_code[:, 1]

    x_au = (
        x_code * LENGTH_UNIT_AU
    )

    y_au = (
        y_code * LENGTH_UNIT_AU
    )

    valid = (
        np.isfinite(x_code)
        & np.isfinite(y_code)
        & np.isfinite(rho_slice)
        & (rho_slice > 0.0)
        & np.isfinite(phi_slice)
        & np.all(
            np.isfinite(velocity_slice),
            axis=1,
        )
        & (np.abs(x_au) <= PLOT_RADIUS_AU)
        & (np.abs(y_au) <= PLOT_RADIUS_AU)
    )

    if not np.any(valid):
        return None

    result = {
        "x_code": x_code[valid],
        "y_code": y_code[valid],
        "rho_code": rho_slice[valid],
        "phi_code": phi_slice[valid],
        "vx": velocity_slice[valid, 0],
        "vy": velocity_slice[valid, 1],
    }

    if "__magnetic__" in sliced.point_data:
        magnetic_slice = np.asarray(
            sliced.point_data["__magnetic__"]
        )

        result["bx"] = magnetic_slice[
            valid, 0
        ]
        result["by"] = magnetic_slice[
            valid, 1
        ]
        result["bz"] = magnetic_slice[
            valid, 2
        ]

    return result


# ============================================================
# AMRブロック境界の重複点を統合
# ============================================================
def merge_duplicate_points(data):
    x = np.asarray(data["x_code"])
    y = np.asarray(data["y_code"])

    coordinate_tolerance = max(
        (
            PLOT_RADIUS_AU
            / LENGTH_UNIT_AU
        )
        * 1.0e-10,
        1.0e-10,
    )

    ix = np.rint(
        x / coordinate_tolerance
    ).astype(np.int64)

    iy = np.rint(
        y / coordinate_tolerance
    ).astype(np.int64)

    keys = np.column_stack((ix, iy))

    _, inverse = np.unique(
        keys,
        axis=0,
        return_inverse=True,
    )

    count = np.bincount(
        inverse
    ).astype(float)

    def average(values):
        return (
            np.bincount(
                inverse,
                weights=np.asarray(values),
            )
            / count
        )

    merged = {
        "x_code": average(data["x_code"]),
        "y_code": average(data["y_code"]),
        "phi_code": average(data["phi_code"]),
        "vx": average(data["vx"]),
        "vy": average(data["vy"]),
    }

    # 密度は対数平均
    merged_log_rho = average(
        np.log10(
            np.maximum(
                data["rho_code"],
                1.0e-300,
            )
        )
    )

    merged["rho_code"] = (
        10.0**merged_log_rho
    )

    # 将来用磁場
    if all(
        key in data
        for key in ("bx", "by", "bz")
    ):
        merged["bx"] = average(data["bx"])
        merged["by"] = average(data["by"])
        merged["bz"] = average(data["bz"])

    return merged


# ============================================================
# 1タイムステップ分を抽出
# ============================================================
def extract_snapshot(step_info):
    required_keys = (
        "x_code",
        "y_code",
        "rho_code",
        "phi_code",
        "vx",
        "vy",
    )

    collected = {
        key: []
        for key in required_keys
    }

    magnetic_collected = {
        "bx": [],
        "by": [],
        "bz": [],
    }

    all_blocks_have_magnetic = True
    used_blocks = 0

    for filename in step_info["files"]:
        try:
            block = extract_xy_midplane(
                filename
            )

            if block is None:
                continue

            used_blocks += 1

            for key in required_keys:
                collected[key].append(
                    block[key]
                )

            if all(
                key in block
                for key in ("bx", "by", "bz")
            ):
                for key in magnetic_collected:
                    magnetic_collected[key].append(
                        block[key]
                    )
            else:
                all_blocks_have_magnetic = False

        except Exception as error:
            print(
                f"[WARNING] "
                f"{Path(filename).name}: "
                f"{error}"
            )

    if not collected["x_code"]:
        raise RuntimeError(
            f"No x-y data at "
            f"step={step_info['step']:05d}"
        )

    raw = {
        key: np.concatenate(values)
        for key, values in collected.items()
    }

    if (
        all_blocks_have_magnetic
        and magnetic_collected["bx"]
    ):
        for key, values in (
            magnetic_collected.items()
        ):
            raw[key] = np.concatenate(values)

    merged = merge_duplicate_points(raw)

    merged["step"] = step_info["step"]
    merged["time_code"] = step_info[
        "time_code"
    ]
    merged["time_myr"] = step_info[
        "time_myr"
    ]

    print(
        f"[INFO] step={step_info['step']:05d}: "
        f"blocks={used_blocks}/"
        f"{len(step_info['files'])}, "
        f"points={len(merged['x_code']):,}, "
        f"B={'yes' if 'bx' in merged else 'no'}"
    )

    return merged


# ============================================================
# 補間
# ============================================================
def interpolate_field(
    points,
    values,
    X,
    Y,
):
    linear = griddata(
        points,
        values,
        (X, Y),
        method="linear",
    )

    if np.any(~np.isfinite(linear)):
        nearest = griddata(
            points,
            values,
            (X, Y),
            method="nearest",
        )

        return np.where(
            np.isfinite(linear),
            linear,
            nearest,
        )

    return linear


# ============================================================
# 各診断量を計算
# ============================================================
def calculate_diagnostics(
    snapshot,
    resolution,
):
    radius_code = (
        PLOT_RADIUS_AU
        / LENGTH_UNIT_AU
    )

    axis_code = np.linspace(
        -radius_code,
        radius_code,
        resolution,
    )

    X_code, Y_code = np.meshgrid(
        axis_code,
        axis_code,
    )

    spacing_code = (
        axis_code[1]
        - axis_code[0]
    )

    points = np.column_stack(
        (
            snapshot["x_code"],
            snapshot["y_code"],
        )
    )

    # 密度は対数で補間
    log_rho_map = interpolate_field(
        points,
        np.log10(
            np.maximum(
                snapshot["rho_code"],
                1.0e-300,
            )
        ),
        X_code,
        Y_code,
    )

    rho_code = (
        10.0**log_rho_map
    )

    phi_code = interpolate_field(
        points,
        snapshot["phi_code"],
        X_code,
        Y_code,
    )

    vx = interpolate_field(
        points,
        snapshot["vx"],
        X_code,
        Y_code,
    )

    vy = interpolate_field(
        points,
        snapshot["vy"],
        X_code,
        Y_code,
    )

    # --------------------------------------------------------
    # 自己重力加速度：g = -grad(phi)
    # np.gradientの返り値は(y方向, x方向)
    # --------------------------------------------------------
    dphi_dy, dphi_dx = np.gradient(
        phi_code,
        spacing_code,
        spacing_code,
        edge_order=2,
    )

    gx_self = -dphi_dx
    gy_self = -dphi_dy

    torque_self = (
        rho_code
        * (
            X_code * gy_self
            - Y_code * gx_self
        )
    )

    # --------------------------------------------------------
    # 等温圧力：P = cs^2 rho
    # --------------------------------------------------------
    pressure_code = (
        ISO_SOUND_SPEED_CODE**2
        * rho_code
    )

    dP_dy, dP_dx = np.gradient(
        pressure_code,
        spacing_code,
        spacing_code,
        edge_order=2,
    )

    # f_P = -grad(P)
    torque_pressure = (
        -X_code * dP_dy
        + Y_code * dP_dx
    )

    # --------------------------------------------------------
    # 移流角運動量流束
    # F_J,adv = rho v_R j_z
    # --------------------------------------------------------
    R_code = np.hypot(
        X_code,
        Y_code,
    )

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        v_R = (
            X_code * vx
            + Y_code * vy
        ) / R_code

    v_R[R_code == 0.0] = 0.0

    j_z = (
        X_code * vy
        - Y_code * vx
    )

    angular_momentum_flux = (
        rho_code
        * v_R
        * j_z
    )

    return {
        "X_code": X_code,
        "Y_code": Y_code,
        "rho_code": rho_code,
        "vx": vx,
        "vy": vy,
        "R_code": R_code,
        "v_R": v_R,
        "j_z": j_z,
        "torque_self": torque_self,
        "torque_pressure": torque_pressure,
        "angular_momentum_flux": (
            angular_momentum_flux
        ),
    }


# ============================================================
# VTKファイル探索
# ============================================================
vtk_files = sorted(
    path
    for path in vtk_dir.glob("*.vtk")
    if (
        "Toyouchi.block" in path.name
        and STREAM_TOKEN in path.name
    )
)

if not vtk_files:
    raise FileNotFoundError(
        f"No VTK files found in {vtk_dir}"
    )

files_by_step = defaultdict(list)

for filename in vtk_files:
    step = get_step_number(filename)

    if step is not None:
        files_by_step[step].append(filename)

if not files_by_step:
    raise RuntimeError(
        "No out2 timestep numbers found."
    )


# ============================================================
# スナップショット一覧
# ============================================================
step_table = []

for step in sorted(files_by_step):
    files = sorted(files_by_step[step])

    block_times = np.array([
        read_vtk_time_code(filename)
        for filename in files
    ])

    time_code = float(
        np.nanmedian(block_times)
    )

    step_table.append({
        "step": step,
        "time_code": time_code,
        "time_myr": (
            time_code
            * TIME_UNIT_MYR
        ),
        "files": files,
    })

step_table.sort(
    key=lambda item: (
        item["time_code"],
        item["step"],
    )
)

print(
    f"[INFO] Snapshot count: "
    f"{len(step_table)}"
)


# ============================================================
# Pass 1：x-y面を抽出
# ============================================================
snapshots = []
density_samples = []

print(
    "\n[INFO] Pass 1/3: "
    "extracting x-y slices"
)

for index, step_info in enumerate(
    step_table
):
    print(
        f"[INFO] Extract "
        f"{index + 1}/{len(step_table)}"
    )

    try:
        snapshot = extract_snapshot(
            step_info
        )
    except Exception as error:
        print(
            f"[WARNING] step="
            f"{step_info['step']:05d}: "
            f"{error}"
        )
        continue

    snapshots.append(snapshot)

    radius_au = (
        np.hypot(
            snapshot["x_code"],
            snapshot["y_code"],
        )
        * LENGTH_UNIT_AU
    )

    density_mask = (
        np.isfinite(snapshot["rho_code"])
        & (snapshot["rho_code"] > 0.0)
        & (radius_au >= SINK_RADIUS_AU)
    )

    density_samples.append(
        snapshot["rho_code"][
            density_mask
        ]
        * DENSITY_UNIT_MSUN_AU3
    )

if not snapshots:
    raise RuntimeError(
        "No valid snapshots extracted."
    )

global_density = np.concatenate(
    density_samples
)

rho_vmin, rho_vmax = np.percentile(
    global_density,
    DENSITY_PERCENTILES,
)

density_norm = LogNorm(
    vmin=rho_vmin,
    vmax=rho_vmax,
    clip=True,
)

print(
    f"[INFO] Density range: "
    f"{rho_vmin:.6e} -- "
    f"{rho_vmax:.6e} M_sun AU^-3"
)


# ============================================================
# Pass 2：全時刻共通のトルク表示範囲
# ============================================================
diagnostic_samples = {
    "torque_self": [],
    "torque_pressure": [],
    "angular_momentum_flux": [],
}

print(
    "\n[INFO] Pass 2/3: "
    "scanning diagnostic ranges"
)

for index, snapshot in enumerate(
    snapshots
):
    print(
        f"[INFO] Scale scan "
        f"{index + 1}/{len(snapshots)}"
    )

    diagnostic = calculate_diagnostics(
        snapshot,
        SCALE_SCAN_RESOLUTION,
    )

    radius_au = (
        np.hypot(
            diagnostic["X_code"],
            diagnostic["Y_code"],
        )
        * LENGTH_UNIT_AU
    )

    outside_sink = (
        radius_au >= SINK_RADIUS_AU
    )

    for key in diagnostic_samples:
        values = diagnostic[key]

        mask = (
            outside_sink
            & np.isfinite(values)
        )

        diagnostic_samples[key].append(
            np.abs(values[mask])
        )

diagnostic_limits = {}

for key, samples in (
    diagnostic_samples.items()
):
    values = np.concatenate(samples)

    limit = np.percentile(
        values,
        TORQUE_ABS_PERCENTILE,
    )

    if (
        not np.isfinite(limit)
        or limit <= 0.0
    ):
        limit = 1.0

    diagnostic_limits[key] = limit

    print(
        f"[INFO] {key}: "
        f"±{limit:.6e}"
    )


# ============================================================
# 描画準備
# ============================================================
axis_formatter = FuncFormatter(
    lambda value, position: (
        f"{value / 1.0e4:g}"
    )
)

panel_definitions = [
    (
        "torque_self",
        "Self-gravity torque density",
        r"$\mathcal{T}_{z,\rm self}$"
        " [code units]",
    ),
    (
        "torque_pressure",
        "Pressure-gradient torque density",
        r"$\mathcal{T}_{z,P}$"
        " [code units]",
    ),
    (
        "angular_momentum_flux",
        "Advective angular-momentum flux",
        r"$F_{J,\rm adv}$"
        " [code units]",
    ),
]

saved_files = []


# ============================================================
# Pass 3：描画・保存
# ============================================================
print(
    "\n[INFO] Pass 3/3: "
    "plotting and saving"
)

for frame_index, snapshot in enumerate(
    snapshots
):
    print(
        f"[INFO] Plot "
        f"{frame_index + 1}/{len(snapshots)}: "
        f"step={snapshot['step']:05d}"
    )

    diagnostic = calculate_diagnostics(
        snapshot,
        MAP_RESOLUTION,
    )

    X_AU = (
        diagnostic["X_code"]
        * LENGTH_UNIT_AU
    )

    Y_AU = (
        diagnostic["Y_code"]
        * LENGTH_UNIT_AU
    )

    rho_display = (
        diagnostic["rho_code"]
        * DENSITY_UNIT_MSUN_AU3
    )

    speed = np.hypot(
        diagnostic["vx"],
        diagnostic["vy"],
    )

    finite_speed = speed[
        np.isfinite(speed)
        & (speed > 0.0)
    ]

    if finite_speed.size > 0:
        speed_threshold = (
            np.median(finite_speed)
            * MIN_SPEED_FRACTION
        )
    else:
        speed_threshold = np.inf

    valid_velocity = (
        np.isfinite(speed)
        & (speed > speed_threshold)
    )

    unit_vx = np.full_like(
        speed,
        np.nan,
    )

    unit_vy = np.full_like(
        speed,
        np.nan,
    )

    unit_vx[valid_velocity] = (
        diagnostic["vx"][valid_velocity]
        / speed[valid_velocity]
    )

    unit_vy[valid_velocity] = (
        diagnostic["vy"][valid_velocity]
        / speed[valid_velocity]
    )

    # --------------------------------------------------------
    # フラックスパネル専用：規格化した動径速度ベクトル
    #
    #   v_R > 0 : 外向きの矢印
    #   v_R < 0 : 中心向きの矢印
    #
    # 矢印長は |v_R| に依存させず、方向と符号だけを表示する。
    # v_R がほぼ0の領域は数値ノイズを避けるため非表示にする。
    # --------------------------------------------------------
    R_code = diagnostic["R_code"]
    v_R = diagnostic["v_R"]

    finite_abs_vr = np.abs(v_R)[
        np.isfinite(v_R)
        & (np.abs(v_R) > 0.0)
    ]

    if finite_abs_vr.size > 0:
        vr_threshold = (
            np.median(finite_abs_vr)
            * MIN_SPEED_FRACTION
        )
    else:
        vr_threshold = np.inf

    valid_radial_velocity = (
        np.isfinite(v_R)
        & np.isfinite(R_code)
        & (R_code > 0.0)
        & (np.abs(v_R) > vr_threshold)
    )

    unit_vr_x = np.full_like(v_R, np.nan)
    unit_vr_y = np.full_like(v_R, np.nan)

    radial_sign = np.sign(v_R)

    unit_vr_x[valid_radial_velocity] = (
        radial_sign[valid_radial_velocity]
        * diagnostic["X_code"][valid_radial_velocity]
        / R_code[valid_radial_velocity]
    )

    unit_vr_y[valid_radial_velocity] = (
        radial_sign[valid_radial_velocity]
        * diagnostic["Y_code"][valid_radial_velocity]
        / R_code[valid_radial_velocity]
    )

    quiver_indices = np.linspace(
        0,
        MAP_RESOLUTION - 1,
        QUIVER_N,
        dtype=int,
    )

    selection = np.ix_(
        quiver_indices,
        quiver_indices,
    )

    arrow_length_au = (
        1.15
        * 2.0
        * PLOT_RADIUS_AU
        / max(QUIVER_N - 1, 1)
    )

    fig, axes = plt.subplots(
        1,
        3,
        figsize=FIGSIZE,
        constrained_layout=True,
    )

    for ax, (
        diagnostic_key,
        title,
        colorbar_label,
    ) in zip(
        axes,
        panel_definitions,
    ):
        # 白黒密度背景
        ax.pcolormesh(
            X_AU,
            Y_AU,
            np.clip(
                rho_display,
                rho_vmin,
                rho_vmax,
            ),
            shading="auto",
            cmap=DENSITY_CMAP,
            norm=density_norm,
            rasterized=True,
        )

        torque_limit = (
            diagnostic_limits[
                diagnostic_key
            ]
        )

        torque_norm = SymLogNorm(
            linthresh=(
                torque_limit
                * 1.0e-3
            ),
            linscale=1.0,
            vmin=-torque_limit,
            vmax=torque_limit,
            base=10,
            clip=True,
        )

        torque_image = ax.pcolormesh(
            X_AU,
            Y_AU,
            np.clip(
                diagnostic[diagnostic_key],
                -torque_limit,
                torque_limit,
            ),
            shading="auto",
            cmap=TORQUE_CMAP,
            norm=torque_norm,
            alpha=TORQUE_ALPHA,
            rasterized=True,
        )

        # 自己重力・圧力パネルは従来の速度方向。
        # 移流フラックスパネルだけ規格化したv_R方向を使う。
        if diagnostic_key == "angular_momentum_flux":
            quiver_u = unit_vr_x
            quiver_v = unit_vr_y
        else:
            quiver_u = unit_vx
            quiver_v = unit_vy

        ax.quiver(
            X_AU[selection],
            Y_AU[selection],
            (
                quiver_u[selection]
                * arrow_length_au
            ),
            (
                quiver_v[selection]
                * arrow_length_au
            ),
            color="black",
            angles="xy",
            scale_units="xy",
            scale=1.0,
            width=0.0025,
            pivot="mid",
            headwidth=3.5,
            headlength=4.5,
            headaxislength=4.0,
            alpha=0.85,
        )

        ax.axhline(
            0.0,
            color="gray",
            lw=0.7,
            alpha=0.7,
        )

        ax.axvline(
            0.0,
            color="gray",
            lw=0.7,
            alpha=0.7,
        )

        ax.plot(
            0.0,
            0.0,
            marker="+",
            color="red",
            markersize=11,
            markeredgewidth=2.5,
        )

        ax.set_xlim(
            -PLOT_RADIUS_AU,
            PLOT_RADIUS_AU,
        )

        ax.set_ylim(
            -PLOT_RADIUS_AU,
            PLOT_RADIUS_AU,
        )

        ax.xaxis.set_major_formatter(
            axis_formatter
        )

        ax.yaxis.set_major_formatter(
            axis_formatter
        )

        ax.set_xlabel(
            r"$x\ [10^4\ {\rm AU}]$"
        )

        ax.set_ylabel(
            r"$y\ [10^4\ {\rm AU}]$"
        )

        ax.set_title(title)
        ax.set_aspect("equal")
        ax.grid(False)

        colorbar = fig.colorbar(
            torque_image,
            ax=ax,
            shrink=0.82,
            pad=0.02,
            extend="both",
        )

        colorbar.set_label(
            colorbar_label
        )

    fig.suptitle(
        "Density, velocity direction, "
        "and angular-momentum transport\n"
        f"step={snapshot['step']:05d}, "
        f"t={snapshot['time_myr']:.6e} Myr",
        fontsize=16,
    )

    output_file = output_dir / (
        f"frame_{frame_index:05d}_"
        f"step_{snapshot['step']:05d}_"
        f"time_{snapshot['time_myr']:.6e}Myr.png"
    )

    fig.savefig(
        output_file,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    saved_files.append(output_file)

    print(
        f"[INFO] Saved: "
        f"{output_file.name}"
    )

print(
    f"\n[INFO] Completed: "
    f"{len(saved_files)} PNG files"
)
print(
    f"[INFO] Output: {output_dir}"
)


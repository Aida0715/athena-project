# -*- coding: utf-8 -*-
# Converted from 速度場on密度map.ipynb

# %% cell 1
# ============================================================
# 全タイムステップ：
# x-z平面 Density + normalized (vx, vz) direction
#
# ・計算領域はVTKの実データ境界から自動取得
# ・y=0と交差するAMRブロックのみ使用
# ・PyVistaで幾何学的なy=0面を抽出
# ・nearest補完なし（データ領域外は白抜き）
# ・距離目盛：10^5 AU単位
# ・時間：整数kyr
# ・密度：M_sun AU^-3
# ============================================================

import gc
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata


# ============================================================
# 1. 入出力設定
# ============================================================
vtk_dir = Path(
    os.path.expanduser(
        "~/athena-project/results/〇〇"
    )
).resolve()

output_dir = Path(
    "./xz_density_normalized_velocity"
).resolve()

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

# 使用するAthena++出力ストリーム
STREAM_TOKEN = ".out2."

print(f"[INFO] Input : {vtk_dir}")
print(f"[INFO] Output: {output_dir}")

if not vtk_dir.is_dir():
    raise FileNotFoundError(
        f"VTK directory does not exist: {vtk_dir}"
    )


# ============================================================
# 2. 単位系
# ============================================================
M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15       # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0
KYR_CGS = 1.0e3 * YEAR_CGS

LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
TIME_UNIT_KYR = T_UNIT_CGS / KYR_CGS
MASS_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS

DENSITY_UNIT_MSUN_AU3 = (
    MASS_UNIT_MSUN
    / LENGTH_UNIT_AU**3
)

# 座標を10^5 AU単位で表示
AXIS_SCALE_AU = 1.0e5

print(
    f"[INFO] 1 code length  = "
    f"{LENGTH_UNIT_AU:.6e} AU"
)
print(
    f"[INFO] 1 code time    = "
    f"{TIME_UNIT_KYR:.6e} kyr"
)
print(
    f"[INFO] 1 code density = "
    f"{DENSITY_UNIT_MSUN_AU3:.6e} "
    "M_sun AU^-3"
)


# ============================================================
# 3. 描画設定
# ============================================================

# 補間画像の一辺のピクセル数
MAP_RESOLUTION = 500

# 各方向の矢印本数
QUIVER_N = 31

# 全時刻共通の密度表示範囲を決めるpercentile
DENSITY_PERCENTILES = (1.0, 99.5)

# 密度カラースケールを決める際に除外するシンク半径
SINK_RADIUS_AU = 1000.0

# 代表速度に対してこれより遅い場所では矢印を非表示
MIN_SPEED_FRACTION = 0.01

# 全時刻の密度統計用に、1時刻から最大何点使用するか
MAX_DENSITY_SAMPLES_PER_STEP = 200000

DENSITY_CMAP_NAME = "Greys"
QUIVER_COLOR = "black"

FIGSIZE = (8, 8)
DPI = 200


# ============================================================
# 4. VTK変数名候補
# ============================================================
DENSITY_CANDIDATES = (
    "rho",
    "dens",
    "density",
    "prim_dens",
    "prim_density",
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


# ============================================================
# 5. 補助関数
# ============================================================
def read_vtk_time_code(filename):
    """VTKヘッダからcode timeを取得する。"""

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
            f"VTK header time not found: {filename}"
        )

    return float(match.group(1))


def get_step_number(filename):
    """ファイル名からout2の出力番号を取得する。"""

    name = Path(filename).name

    match = re.search(
        r"(?:prim\.)?out2\.(\d+)",
        name,
    )

    if match is None:
        return None

    return int(match.group(1))


def find_point_field(dataset, candidates):
    """point_dataから候補に一致する変数名を探す。"""

    lower_to_original = {
        name.lower(): name
        for name in dataset.point_data.keys()
    }

    for candidate in candidates:
        key = candidate.lower()

        if key in lower_to_original:
            return lower_to_original[key]

    return None


# ============================================================
# 6. VTKファイル探索
# ============================================================
vtk_files = sorted(
    path
    for path in vtk_dir.glob("*.vtk")
    if "Toyouchi.block" in path.name
)

if STREAM_TOKEN is not None:
    vtk_files = [
        path
        for path in vtk_files
        if STREAM_TOKEN in path.name
    ]

if not vtk_files:
    raise FileNotFoundError(
        f"No VTK files found in {vtk_dir}"
    )

print(
    f"[INFO] VTK file count: {len(vtk_files)}"
)


# ============================================================
# 7. stepごとにVTKファイルを分類
# ============================================================
files_by_step = defaultdict(list)

for filename in vtk_files:
    step = get_step_number(filename)

    if step is not None:
        files_by_step[step].append(filename)

if not files_by_step:
    raise RuntimeError(
        "No out2 timestep numbers were found."
    )


# ============================================================
# 8. スナップショット一覧を作成
# ============================================================
step_table = []

for step in sorted(files_by_step):
    files = sorted(files_by_step[step])

    block_times = np.array(
        [
            read_vtk_time_code(filename)
            for filename in files
        ],
        dtype=float,
    )

    time_code = float(
        np.nanmedian(block_times)
    )

    time_spread = float(
        np.nanmax(block_times)
        - np.nanmin(block_times)
    )

    allowed_spread = max(
        1.0e-10,
        abs(time_code) * 1.0e-10,
    )

    if time_spread > allowed_spread:
        print(
            f"[WARNING] step={step:05d}: "
            f"block times differ by "
            f"{time_spread:.6e} code time"
        )

    time_kyr_float = (
        time_code
        * TIME_UNIT_KYR
    )

    # 表示用の整数kyr
    time_kyr_integer = int(
        np.rint(time_kyr_float)
    )

    step_table.append({
        "step": step,
        "time_code": time_code,
        "time_kyr": time_kyr_float,
        "time_kyr_integer": time_kyr_integer,
        "files": files,
    })

step_table.sort(
    key=lambda item: (
        item["time_code"],
        item["step"],
    )
)

print(
    f"[INFO] Snapshot count: {len(step_table)}"
)
print(
    f"[INFO] First time: "
    f"{step_table[0]['time_kyr_integer']} kyr"
)
print(
    f"[INFO] Last time : "
    f"{step_table[-1]['time_kyr_integer']} kyr"
)


# ============================================================
# 9. 実際のx-z計算領域をVTK境界から自動取得
#
# 最初のstepに属する全ブロックのうち、y=0面と交差する
# ブロックだけを調べる。
# ============================================================
def determine_xz_domain(files):
    x_min_code = np.inf
    x_max_code = -np.inf
    z_min_code = np.inf
    z_max_code = -np.inf

    used_blocks = 0

    for filename in files:
        grid = pv.read(filename)

        xmin, xmax, ymin, ymax, zmin, zmax = (
            grid.bounds
        )

        tolerance = (
            1.0e-12
            * max(
                abs(ymin),
                abs(ymax),
                1.0,
            )
        )

        if (
            ymin - tolerance
            <= 0.0
            <= ymax + tolerance
        ):
            x_min_code = min(
                x_min_code,
                xmin,
            )
            x_max_code = max(
                x_max_code,
                xmax,
            )
            z_min_code = min(
                z_min_code,
                zmin,
            )
            z_max_code = max(
                z_max_code,
                zmax,
            )

            used_blocks += 1

        del grid

    if (
        used_blocks == 0
        or not np.all(
            np.isfinite([
                x_min_code,
                x_max_code,
                z_min_code,
                z_max_code,
            ])
        )
    ):
        raise RuntimeError(
            "Could not determine the x-z domain "
            "intersecting y=0."
        )

    return {
        "x_min_au": (
            x_min_code
            * LENGTH_UNIT_AU
        ),
        "x_max_au": (
            x_max_code
            * LENGTH_UNIT_AU
        ),
        "z_min_au": (
            z_min_code
            * LENGTH_UNIT_AU
        ),
        "z_max_au": (
            z_max_code
            * LENGTH_UNIT_AU
        ),
        "used_blocks": used_blocks,
    }


domain = determine_xz_domain(
    step_table[0]["files"]
)

X_MIN_AU = domain["x_min_au"]
X_MAX_AU = domain["x_max_au"]
Z_MIN_AU = domain["z_min_au"]
Z_MAX_AU = domain["z_max_au"]

DOMAIN_WIDTH_X_AU = X_MAX_AU - X_MIN_AU
DOMAIN_WIDTH_Z_AU = Z_MAX_AU - Z_MIN_AU

if (
    DOMAIN_WIDTH_X_AU <= 0.0
    or DOMAIN_WIDTH_Z_AU <= 0.0
):
    raise RuntimeError(
        "Invalid automatically determined domain."
    )

print(
    "[INFO] Automatically determined x-z domain:"
)
print(
    f"       x = {X_MIN_AU:.6e} -- "
    f"{X_MAX_AU:.6e} AU"
)
print(
    f"       z = {Z_MIN_AU:.6e} -- "
    f"{Z_MAX_AU:.6e} AU"
)
print(
    f"       width = "
    f"{DOMAIN_WIDTH_X_AU:.6e} × "
    f"{DOMAIN_WIDTH_Z_AU:.6e} AU"
)
print(
    f"       y=0 intersecting blocks = "
    f"{domain['used_blocks']}"
)


# ============================================================
# 10. 正確なy=0面をブロックから抽出
# ============================================================
def extract_block_xz_slice(filename):
    """
    cell_dataをpoint_dataへ変換した後、PyVistaのsliceで
    幾何学的なy=0面を抽出する。
    """

    grid = pv.read(filename)

    xmin, xmax, ymin, ymax, zmin, zmax = (
        grid.bounds
    )

    bounds_tolerance = (
        1.0e-12
        * max(
            abs(ymin),
            abs(ymax),
            1.0,
        )
    )

    if not (
        ymin - bounds_tolerance
        <= 0.0
        <= ymax + bounds_tolerance
    ):
        del grid
        return None

    point_grid = grid.cell_data_to_point_data(
        pass_cell_data=False,
    )

    sliced = point_grid.slice(
        normal=(0.0, 1.0, 0.0),
        origin=(0.0, 0.0, 0.0),
    )

    if sliced.n_points == 0:
        del sliced
        del point_grid
        del grid
        return None

    density_name = find_point_field(
        sliced,
        DENSITY_CANDIDATES,
    )

    if density_name is None:
        available = list(
            sliced.point_data.keys()
        )

        raise KeyError(
            f"Density field not found: {filename}\n"
            f"point_data={available}"
        )

    density_code = np.asarray(
        sliced.point_data[density_name]
    ).reshape(-1)

    velocity_name = find_point_field(
        sliced,
        VELOCITY_CANDIDATES,
    )

    if velocity_name is not None:
        velocity = np.asarray(
            sliced.point_data[velocity_name]
        )

    else:
        momentum_name = find_point_field(
            sliced,
            MOMENTUM_CANDIDATES,
        )

        if momentum_name is None:
            available = list(
                sliced.point_data.keys()
            )

            raise KeyError(
                f"Velocity/momentum not found: "
                f"{filename}\n"
                f"point_data={available}"
            )

        momentum = np.asarray(
            sliced.point_data[momentum_name]
        )

        velocity = (
            momentum
            / np.maximum(
                density_code[:, None],
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

    points = np.asarray(
        sliced.points
    )

    x_au = (
        points[:, 0]
        * LENGTH_UNIT_AU
    )

    z_au = (
        points[:, 2]
        * LENGTH_UNIT_AU
    )

    density = (
        density_code
        * DENSITY_UNIT_MSUN_AU3
    )

    vx = velocity[:, 0]
    vz = velocity[:, 2]

    # 浮動小数点誤差だけを考慮した領域判定
    x_tolerance = max(
        DOMAIN_WIDTH_X_AU * 1.0e-12,
        1.0e-8,
    )

    z_tolerance = max(
        DOMAIN_WIDTH_Z_AU * 1.0e-12,
        1.0e-8,
    )

    valid = (
        np.isfinite(x_au)
        & np.isfinite(z_au)
        & np.isfinite(density)
        & (density > 0.0)
        & np.isfinite(vx)
        & np.isfinite(vz)
        & (
            x_au
            >= X_MIN_AU - x_tolerance
        )
        & (
            x_au
            <= X_MAX_AU + x_tolerance
        )
        & (
            z_au
            >= Z_MIN_AU - z_tolerance
        )
        & (
            z_au
            <= Z_MAX_AU + z_tolerance
        )
    )

    result = None

    if np.any(valid):
        result = {
            "x": x_au[valid].copy(),
            "z": z_au[valid].copy(),
            "rho": density[valid].copy(),
            "vx": vx[valid].copy(),
            "vz": vz[valid].copy(),
        }

    del sliced
    del point_grid
    del grid

    return result


# ============================================================
# 11. AMRブロック境界の重複点を統合
# ============================================================
def merge_duplicate_points(data):
    """
    同じ(x,z)に複数ブロックの値が存在する場合に統合する。

    密度：対数平均
    速度：算術平均
    """

    x = np.asarray(data["x"])
    z = np.asarray(data["z"])
    rho = np.asarray(data["rho"])
    vx = np.asarray(data["vx"])
    vz = np.asarray(data["vz"])

    domain_scale = max(
        DOMAIN_WIDTH_X_AU,
        DOMAIN_WIDTH_Z_AU,
    )

    coordinate_tolerance = max(
        domain_scale * 1.0e-10,
        1.0e-8,
    )

    ix = np.rint(
        x / coordinate_tolerance
    ).astype(np.int64)

    iz = np.rint(
        z / coordinate_tolerance
    ).astype(np.int64)

    keys = np.column_stack(
        (ix, iz)
    )

    _, inverse = np.unique(
        keys,
        axis=0,
        return_inverse=True,
    )

    count = np.bincount(
        inverse
    ).astype(float)

    x_merged = (
        np.bincount(
            inverse,
            weights=x,
        )
        / count
    )

    z_merged = (
        np.bincount(
            inverse,
            weights=z,
        )
        / count
    )

    log_rho_merged = (
        np.bincount(
            inverse,
            weights=np.log10(
                np.maximum(
                    rho,
                    1.0e-300,
                )
            ),
        )
        / count
    )

    rho_merged = (
        10.0**log_rho_merged
    )

    vx_merged = (
        np.bincount(
            inverse,
            weights=vx,
        )
        / count
    )

    vz_merged = (
        np.bincount(
            inverse,
            weights=vz,
        )
        / count
    )

    return {
        "x": x_merged,
        "z": z_merged,
        "rho": rho_merged,
        "vx": vx_merged,
        "vz": vz_merged,
    }


# ============================================================
# 12. 1スナップショットのx-z面を抽出
# ============================================================
def extract_snapshot(step_info):
    collected = {
        "x": [],
        "z": [],
        "rho": [],
        "vx": [],
        "vz": [],
    }

    n_used_blocks = 0

    for filename in step_info["files"]:
        try:
            block_data = (
                extract_block_xz_slice(
                    filename
                )
            )

            if block_data is None:
                continue

            n_used_blocks += 1

            for key in collected:
                collected[key].append(
                    block_data[key]
                )

            del block_data

        except Exception as error:
            print(
                f"[WARNING] "
                f"{Path(filename).name}: "
                f"{error}"
            )

    if not collected["x"]:
        raise RuntimeError(
            f"No x-z slice data at "
            f"step={step_info['step']:05d}"
        )

    raw_data = {
        key: np.concatenate(values)
        for key, values in collected.items()
    }

    merged_data = merge_duplicate_points(
        raw_data
    )

    merged_data["step"] = (
        step_info["step"]
    )
    merged_data["time_code"] = (
        step_info["time_code"]
    )
    merged_data["time_kyr"] = (
        step_info["time_kyr"]
    )
    merged_data["time_kyr_integer"] = (
        step_info["time_kyr_integer"]
    )

    print(
        f"[INFO] step="
        f"{step_info['step']:05d}: "
        f"blocks={n_used_blocks}/"
        f"{len(step_info['files'])}, "
        f"raw={len(raw_data['x']):,}, "
        f"merged={len(merged_data['x']):,}"
    )

    del raw_data
    del collected

    return merged_data


# ============================================================
# 13. 線形補間
#
# nearest補完は行わない。
# 凸包の外側はNaNのまま保持する。
# ============================================================
def interpolate_field_linear(
    points_2d,
    values,
    X,
    Z,
):
    values = np.asarray(
        values
    )

    return griddata(
        points_2d,
        values,
        (X, Z),
        method="linear",
        fill_value=np.nan,
    )


# ============================================================
# 14. 全時刻共通の密度カラースケール
#
# 1スナップショットずつ読み込み、統計値だけ保存する。
# 全スナップショットをメモリに保持しない。
# ============================================================
density_samples = []

for index, step_info in enumerate(
    step_table
):
    print(
        f"[INFO] Density scan "
        f"{index + 1}/"
        f"{len(step_table)}: "
        f"step={step_info['step']:05d}"
    )

    try:
        snapshot = extract_snapshot(
            step_info
        )

        radius_au = np.hypot(
            snapshot["x"],
            snapshot["z"],
        )

        mask = (
            np.isfinite(
                snapshot["rho"]
            )
            & (
                snapshot["rho"]
                > 0.0
            )
            & (
                radius_au
                >= SINK_RADIUS_AU
            )
        )

        values = snapshot["rho"][mask]

        if values.size > 0:
            if (
                values.size
                > MAX_DENSITY_SAMPLES_PER_STEP
            ):
                sample_indices = np.linspace(
                    0,
                    values.size - 1,
                    MAX_DENSITY_SAMPLES_PER_STEP,
                    dtype=int,
                )

                values = values[
                    sample_indices
                ]

            density_samples.append(
                values.copy()
            )

        del radius_au
        del mask
        del values
        del snapshot
        gc.collect()

    except Exception as error:
        print(
            f"[WARNING] Density scan failed, "
            f"step={step_info['step']:05d}: "
            f"{error}"
        )

if not density_samples:
    raise RuntimeError(
        "No positive density values were found "
        "outside the sink."
    )

global_density = np.concatenate(
    density_samples
)

rho_vmin, rho_vmax = np.percentile(
    global_density,
    DENSITY_PERCENTILES,
)

del global_density
del density_samples
gc.collect()

if (
    not np.isfinite(rho_vmin)
    or not np.isfinite(rho_vmax)
    or rho_vmin <= 0.0
    or rho_vmax <= rho_vmin
):
    raise RuntimeError(
        "Invalid global density range: "
        f"{rho_vmin}, {rho_vmax}"
    )

density_norm = LogNorm(
    vmin=rho_vmin,
    vmax=rho_vmax,
    clip=True,
)

print(
    f"[INFO] Global density range: "
    f"{rho_vmin:.6e} -- "
    f"{rho_vmax:.6e} "
    "M_sun AU^-3"
)


# ============================================================
# 15. 実データ領域に対応する補間グリッド
# ============================================================
x_axis_au = np.linspace(
    X_MIN_AU,
    X_MAX_AU,
    MAP_RESOLUTION,
)

z_axis_au = np.linspace(
    Z_MIN_AU,
    Z_MAX_AU,
    MAP_RESOLUTION,
)

X_map_au, Z_map_au = np.meshgrid(
    x_axis_au,
    z_axis_au,
)

# 描画座標は10^5 AU単位
X_map_plot = (
    X_map_au
    / AXIS_SCALE_AU
)

Z_map_plot = (
    Z_map_au
    / AXIS_SCALE_AU
)

X_MIN_PLOT = X_MIN_AU / AXIS_SCALE_AU
X_MAX_PLOT = X_MAX_AU / AXIS_SCALE_AU
Z_MIN_PLOT = Z_MIN_AU / AXIS_SCALE_AU
Z_MAX_PLOT = Z_MAX_AU / AXIS_SCALE_AU

# NaN領域を白で表示
density_cmap = plt.get_cmap(
    DENSITY_CMAP_NAME
).copy()

density_cmap.set_bad(
    color="white"
)


# ============================================================
# 16. 全時刻を1枚ずつ描画・保存
# ============================================================
saved_files = []

for frame_index, step_info in enumerate(
    step_table
):
    print(
        f"[INFO] Plotting "
        f"{frame_index + 1}/"
        f"{len(step_table)}: "
        f"step={step_info['step']:05d}"
    )

    try:
        snapshot = extract_snapshot(
            step_info
        )

        points_xz = np.column_stack(
            (
                snapshot["x"],
                snapshot["z"],
            )
        )

        # 密度はlog10(rho)を線形補間
        log_rho_map = (
            interpolate_field_linear(
                points_xz,
                np.log10(
                    np.maximum(
                        snapshot["rho"],
                        1.0e-300,
                    )
                ),
                X_map_au,
                Z_map_au,
            )
        )

        rho_map = np.full_like(
            log_rho_map,
            np.nan,
        )

        finite_log_rho = np.isfinite(
            log_rho_map
        )

        rho_map[finite_log_rho] = (
            10.0
            ** log_rho_map[
                finite_log_rho
            ]
        )

        vx_map = (
            interpolate_field_linear(
                points_xz,
                snapshot["vx"],
                X_map_au,
                Z_map_au,
            )
        )

        vz_map = (
            interpolate_field_linear(
                points_xz,
                snapshot["vz"],
                X_map_au,
                Z_map_au,
            )
        )

        # ----------------------------------------------------
        # 補間後に速度を規格化
        # ----------------------------------------------------
        speed_map = np.hypot(
            vx_map,
            vz_map,
        )

        finite_speed = speed_map[
            np.isfinite(speed_map)
            & (speed_map > 0.0)
        ]

        if finite_speed.size > 0:
            representative_speed = (
                np.median(finite_speed)
            )

            minimum_speed = (
                representative_speed
                * MIN_SPEED_FRACTION
            )
        else:
            minimum_speed = np.inf

        valid_velocity = (
            np.isfinite(vx_map)
            & np.isfinite(vz_map)
            & np.isfinite(speed_map)
            & (
                speed_map
                > minimum_speed
            )
        )

        unit_vx = np.full_like(
            vx_map,
            np.nan,
        )

        unit_vz = np.full_like(
            vz_map,
            np.nan,
        )

        unit_vx[valid_velocity] = (
            vx_map[valid_velocity]
            / speed_map[valid_velocity]
        )

        unit_vz[valid_velocity] = (
            vz_map[valid_velocity]
            / speed_map[valid_velocity]
        )

        # ----------------------------------------------------
        # 矢印を間引く
        # ----------------------------------------------------
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

        X_quiver = (
            X_map_plot[selection]
        )

        Z_quiver = (
            Z_map_plot[selection]
        )

        U_quiver = (
            unit_vx[selection]
        )

        W_quiver = (
            unit_vz[selection]
        )

        # 描画座標単位での矢印長
        arrow_length_plot = (
            1.25
            * min(
                X_MAX_PLOT
                - X_MIN_PLOT,
                Z_MAX_PLOT
                - Z_MIN_PLOT,
            )
            / max(
                QUIVER_N - 1,
                1,
            )
        )

        # ----------------------------------------------------
        # 描画
        # ----------------------------------------------------
        fig, ax = plt.subplots(
            figsize=FIGSIZE
        )

        rho_masked = np.ma.masked_invalid(
            rho_map
        )

        ax.pcolormesh(
            X_map_plot,
            Z_map_plot,
            rho_masked,
            shading="auto",
            cmap=density_cmap,
            norm=density_norm,
            rasterized=True,
        )

        valid_quiver = (
            np.isfinite(U_quiver)
            & np.isfinite(W_quiver)
        )

        ax.quiver(
            X_quiver[valid_quiver],
            Z_quiver[valid_quiver],
            (
                U_quiver[valid_quiver]
                * arrow_length_plot
            ),
            (
                W_quiver[valid_quiver]
                * arrow_length_plot
            ),
            color=QUIVER_COLOR,
            angles="xy",
            scale_units="xy",
            scale=1.0,
            width=0.0025,
            pivot="mid",
            headwidth=3.5,
            headlength=4.5,
            headaxislength=4.0,
            alpha=0.9,
        )

        ax.axhline(
            0.0,
            color="gray",
            linewidth=0.7,
            alpha=0.7,
        )

        ax.axvline(
            0.0,
            color="gray",
            linewidth=0.7,
            alpha=0.7,
        )

        ax.plot(
            0.0,
            0.0,
            marker="+",
            color="red",
            markersize=12,
            markeredgewidth=2.5,
        )

        ax.set_xlim(
            X_MIN_PLOT,
            X_MAX_PLOT,
        )

        ax.set_ylim(
            Z_MIN_PLOT,
            Z_MAX_PLOT,
        )

        # 座標軸名とスケールを明記
        ax.set_xlabel(
            r"$x\ [10^5\ {\rm AU}]$",
            fontsize=13,
        )

        ax.set_ylabel(
            r"$z\ [10^5\ {\rm AU}]$",
            fontsize=13,
        )

        ax.set_title(
            "Density + normalized "
            r"$(v_x,v_z)$ direction"
            "\n"
            f"step={snapshot['step']:05d}, "
            f"t={snapshot['time_kyr_integer']} kyr",
            fontsize=14,
        )

        ax.set_aspect(
            "equal",
            adjustable="box",
        )

        ax.grid(False)

        fig.tight_layout()

        # 整数kyrをファイル名にも使用
        output_file = output_dir / (
            f"frame_{frame_index:05d}_"
            f"step_{snapshot['step']:05d}_"
            f"time_{snapshot['time_kyr_integer']:08d}kyr.png"
        )

        fig.savefig(
            output_file,
            dpi=DPI,
            bbox_inches="tight",
            facecolor="white",
        )

        plt.close(fig)

        saved_files.append(
            output_file
        )

        print(
            f"[INFO] Saved: "
            f"{output_file.name}"
        )

        del snapshot
        del points_xz
        del log_rho_map
        del rho_map
        del vx_map
        del vz_map
        del speed_map
        del finite_speed
        del valid_velocity
        del unit_vx
        del unit_vz
        del rho_masked
        del fig
        del ax

        gc.collect()

    except Exception as error:
        plt.close("all")
        gc.collect()

        print(
            f"[WARNING] Plot failed, "
            f"step={step_info['step']:05d}: "
            f"{error}"
        )


# ============================================================
# 17. 終了表示
# ============================================================
print(
    f"\n[INFO] Completed: "
    f"{len(saved_files)} figures"
)

print(
    f"[INFO] Output directory: "
    f"{output_dir}"
)

if saved_files:
    print(
        f"[INFO] First output: "
        f"{saved_files[0]}"
    )
    print(
        f"[INFO] Last output : "
        f"{saved_files[-1]}"
    )

# %% cell 2
# ============================================================
# 全タイムステップ：
# x-y平面 Density + normalized (vx, vy) direction
#
# ・計算領域はVTKの実データ境界から自動取得
# ・z=0と交差するAMRブロックのみ使用
# ・PyVistaで幾何学的なz=0面を抽出
# ・AMR境界の重複座標を統合
# ・nearest補完なし（実データ領域外は白抜き）
# ・全時刻で共通の密度カラースケール
# ・距離目盛：10^5 AU単位
# ・時間表示：整数kyr
# ・密度：M_sun AU^-3
# ・メモリ節約のため1時刻ずつ処理
# ============================================================

import gc
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata


# ============================================================
# 1. 入出力設定
# ============================================================
vtk_dir = Path(
    os.path.expanduser(
        "~/athena-project/results/〇〇"
    )
).resolve()

output_dir = (
    vtk_dir
    / "xy_density_normalized_velocity"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

# 使用するAthena++出力ストリーム
STREAM_TOKEN = ".out2."

print(f"[INFO] Input : {vtk_dir}")
print(f"[INFO] Output: {output_dir}")

if not vtk_dir.is_dir():
    raise FileNotFoundError(
        f"VTK directory does not exist: {vtk_dir}"
    )


# ============================================================
# 2. 単位系
# ============================================================
M_UNIT_CGS = 4.0e33
L_UNIT_CGS = 7.03e15       # L0 = rb / 2
T_UNIT_CGS = 3.61e10

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0
KYR_CGS = 1.0e3 * YEAR_CGS

LENGTH_UNIT_AU = (
    L_UNIT_CGS
    / AU_CGS
)

TIME_UNIT_KYR = (
    T_UNIT_CGS
    / KYR_CGS
)

MASS_UNIT_MSUN = (
    M_UNIT_CGS
    / MSUN_CGS
)

DENSITY_UNIT_MSUN_AU3 = (
    MASS_UNIT_MSUN
    / LENGTH_UNIT_AU**3
)

# 描画座標を10^5 AU単位にする
AXIS_SCALE_AU = 1.0e5

print(
    f"[INFO] 1 code length  = "
    f"{LENGTH_UNIT_AU:.6e} AU"
)
print(
    f"[INFO] 1 code time    = "
    f"{TIME_UNIT_KYR:.6e} kyr"
)
print(
    f"[INFO] 1 code density = "
    f"{DENSITY_UNIT_MSUN_AU3:.6e} "
    "M_sun AU^-3"
)


# ============================================================
# 3. 解析・描画設定
# ============================================================

# 補間画像の解像度
MAP_RESOLUTION = 700

# 各方向の矢印本数
QUIVER_N = 31

# 全時刻の密度分布から使用するpercentile
DENSITY_PERCENTILES = (
    5.0,
    99.0,
)

# 密度表示範囲の統計から除外するシンク半径
SINK_RADIUS_AU = 1000.0

# 代表速度に対してこれより遅い場所では矢印を表示しない
MIN_SPEED_FRACTION = 0.01

# 1時刻から密度統計に使用する最大点数
MAX_DENSITY_SAMPLES_PER_STEP = 100000

DENSITY_CMAP_NAME = "Greys"
QUIVER_COLOR = "black"

FIGSIZE = (8, 8)
DPI = 200


# ============================================================
# 4. VTKフィールド名候補
# ============================================================
DENSITY_CANDIDATES = (
    "rho",
    "dens",
    "density",
    "prim_dens",
    "prim_density",
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


# ============================================================
# 5. VTKヘッダからcode timeを取得
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


# ============================================================
# 6. ファイル名から出力stepを取得
# ============================================================
def get_step_number(filename):
    match = re.search(
        r"(?:prim\.)?out2\.(\d+)",
        Path(filename).name,
    )

    if match is None:
        return None

    return int(match.group(1))


# ============================================================
# 7. point_dataのフィールド検索
# ============================================================
def find_point_field(dataset, candidates):
    lower_to_original = {
        name.lower(): name
        for name in dataset.point_data.keys()
    }

    for candidate in candidates:
        key = candidate.lower()

        if key in lower_to_original:
            return lower_to_original[key]

    return None


# ============================================================
# 8. VTKファイルを探索
# ============================================================
vtk_files = sorted(
    path
    for path in vtk_dir.glob("*.vtk")
    if "Toyouchi.block" in path.name
)

if STREAM_TOKEN is not None:
    vtk_files = [
        path
        for path in vtk_files
        if STREAM_TOKEN in path.name
    ]

if not vtk_files:
    raise FileNotFoundError(
        f"No VTK files found in {vtk_dir}"
    )

print(
    f"[INFO] VTK file count: "
    f"{len(vtk_files)}"
)


# ============================================================
# 9. out2番号ごとにブロックを分類
# ============================================================
files_by_step = defaultdict(list)

for filename in vtk_files:
    step = get_step_number(filename)

    if step is not None:
        files_by_step[step].append(filename)

if not files_by_step:
    raise RuntimeError(
        "No out2 timestep numbers were found."
    )


# ============================================================
# 10. スナップショット一覧を作成
# ============================================================
step_table = []

for step in sorted(files_by_step):
    files = sorted(
        files_by_step[step]
    )

    block_times = np.array(
        [
            read_vtk_time_code(filename)
            for filename in files
        ],
        dtype=float,
    )

    time_code = float(
        np.nanmedian(block_times)
    )

    time_spread = float(
        np.nanmax(block_times)
        - np.nanmin(block_times)
    )

    allowed_spread = max(
        1.0e-10,
        abs(time_code) * 1.0e-10,
    )

    if time_spread > allowed_spread:
        print(
            f"[WARNING] step={step:05d}: "
            f"block-time spread="
            f"{time_spread:.6e} code time"
        )

    time_kyr = (
        time_code
        * TIME_UNIT_KYR
    )

    step_table.append({
        "step": step,
        "time_code": time_code,
        "time_kyr": time_kyr,
        "time_kyr_integer": int(
            np.rint(time_kyr)
        ),
        "files": files,
    })

# 物理時間順に並べる
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
print(
    f"[INFO] First time: "
    f"{step_table[0]['time_kyr_integer']} kyr"
)
print(
    f"[INFO] Last time : "
    f"{step_table[-1]['time_kyr_integer']} kyr"
)


# ============================================================
# 11. 実際のx-y計算領域をVTK境界から自動取得
#
# 最初のstepに属する全ブロックのうち、z=0面と交差する
# ブロックのx・y境界を取得する。
# ============================================================
def determine_xy_domain(files):
    x_min_code = np.inf
    x_max_code = -np.inf
    y_min_code = np.inf
    y_max_code = -np.inf

    used_blocks = 0

    for filename in files:
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

        if (
            zmin - tolerance
            <= 0.0
            <= zmax + tolerance
        ):
            x_min_code = min(
                x_min_code,
                xmin,
            )
            x_max_code = max(
                x_max_code,
                xmax,
            )
            y_min_code = min(
                y_min_code,
                ymin,
            )
            y_max_code = max(
                y_max_code,
                ymax,
            )

            used_blocks += 1

        del grid

    domain_values = np.array(
        [
            x_min_code,
            x_max_code,
            y_min_code,
            y_max_code,
        ]
    )

    if (
        used_blocks == 0
        or not np.all(
            np.isfinite(domain_values)
        )
    ):
        raise RuntimeError(
            "Could not determine the x-y domain "
            "intersecting z=0."
        )

    return {
        "x_min_au": (
            x_min_code
            * LENGTH_UNIT_AU
        ),
        "x_max_au": (
            x_max_code
            * LENGTH_UNIT_AU
        ),
        "y_min_au": (
            y_min_code
            * LENGTH_UNIT_AU
        ),
        "y_max_au": (
            y_max_code
            * LENGTH_UNIT_AU
        ),
        "used_blocks": used_blocks,
    }


domain = determine_xy_domain(
    step_table[0]["files"]
)

X_MIN_AU = domain["x_min_au"]
X_MAX_AU = domain["x_max_au"]
Y_MIN_AU = domain["y_min_au"]
Y_MAX_AU = domain["y_max_au"]

DOMAIN_WIDTH_X_AU = (
    X_MAX_AU
    - X_MIN_AU
)

DOMAIN_WIDTH_Y_AU = (
    Y_MAX_AU
    - Y_MIN_AU
)

if (
    DOMAIN_WIDTH_X_AU <= 0.0
    or DOMAIN_WIDTH_Y_AU <= 0.0
):
    raise RuntimeError(
        "Invalid automatically determined domain."
    )

print(
    "[INFO] Automatically determined x-y domain:"
)
print(
    f"       x = {X_MIN_AU:.6e} -- "
    f"{X_MAX_AU:.6e} AU"
)
print(
    f"       y = {Y_MIN_AU:.6e} -- "
    f"{Y_MAX_AU:.6e} AU"
)
print(
    f"       width = "
    f"{DOMAIN_WIDTH_X_AU:.6e} × "
    f"{DOMAIN_WIDTH_Y_AU:.6e} AU"
)
print(
    f"       z=0 intersecting blocks = "
    f"{domain['used_blocks']}"
)


# ============================================================
# 12. 各ブロックから正確なz=0面を抽出
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

    # z=0と交差しないブロックは除外
    if not (
        zmin - tolerance
        <= 0.0
        <= zmax + tolerance
    ):
        del grid
        return None

    # cell_dataをpoint_dataへ補間
    point_grid = (
        grid.cell_data_to_point_data(
            pass_cell_data=False,
        )
    )

    # 幾何学的に正確なz=0面を抽出
    sliced = point_grid.slice(
        normal=(0.0, 0.0, 1.0),
        origin=(0.0, 0.0, 0.0),
    )

    if sliced.n_points == 0:
        del sliced
        del point_grid
        del grid
        return None

    density_name = find_point_field(
        sliced,
        DENSITY_CANDIDATES,
    )

    if density_name is None:
        available = list(
            sliced.point_data.keys()
        )

        raise KeyError(
            f"Density field not found: {filename}\n"
            f"point_data={available}"
        )

    density_code = np.asarray(
        sliced.point_data[density_name]
    ).reshape(-1)

    velocity_name = find_point_field(
        sliced,
        VELOCITY_CANDIDATES,
    )

    if velocity_name is not None:
        velocity = np.asarray(
            sliced.point_data[velocity_name]
        )

    else:
        momentum_name = find_point_field(
            sliced,
            MOMENTUM_CANDIDATES,
        )

        if momentum_name is None:
            available = list(
                sliced.point_data.keys()
            )

            raise KeyError(
                f"Velocity/momentum not found: "
                f"{filename}\n"
                f"point_data={available}"
            )

        momentum = np.asarray(
            sliced.point_data[momentum_name]
        )

        velocity = (
            momentum
            / np.maximum(
                density_code[:, None],
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

    points_code = np.asarray(
        sliced.points
    )

    x_au = (
        points_code[:, 0]
        * LENGTH_UNIT_AU
    )

    y_au = (
        points_code[:, 1]
        * LENGTH_UNIT_AU
    )

    density = (
        density_code
        * DENSITY_UNIT_MSUN_AU3
    )

    vx = velocity[:, 0]
    vy = velocity[:, 1]

    x_tolerance = max(
        DOMAIN_WIDTH_X_AU * 1.0e-12,
        1.0e-8,
    )

    y_tolerance = max(
        DOMAIN_WIDTH_Y_AU * 1.0e-12,
        1.0e-8,
    )

    valid = (
        np.isfinite(x_au)
        & np.isfinite(y_au)
        & np.isfinite(density)
        & (density > 0.0)
        & np.isfinite(vx)
        & np.isfinite(vy)
        & (
            x_au
            >= X_MIN_AU - x_tolerance
        )
        & (
            x_au
            <= X_MAX_AU + x_tolerance
        )
        & (
            y_au
            >= Y_MIN_AU - y_tolerance
        )
        & (
            y_au
            <= Y_MAX_AU + y_tolerance
        )
    )

    result = None

    if np.any(valid):
        result = {
            "x": x_au[valid].copy(),
            "y": y_au[valid].copy(),
            "rho": density[valid].copy(),
            "vx": vx[valid].copy(),
            "vy": vy[valid].copy(),
        }

    del sliced
    del point_grid
    del grid

    return result


# ============================================================
# 13. AMRブロック境界の重複座標を統合
# ============================================================
def merge_duplicate_points(data):
    """
    同じ(x,y)に複数ブロックの値が存在する場合に統合する。

    密度：対数平均
    速度：算術平均
    """

    x = np.asarray(data["x"])
    y = np.asarray(data["y"])
    rho = np.asarray(data["rho"])
    vx = np.asarray(data["vx"])
    vy = np.asarray(data["vy"])

    domain_scale = max(
        DOMAIN_WIDTH_X_AU,
        DOMAIN_WIDTH_Y_AU,
    )

    coordinate_tolerance = max(
        domain_scale * 1.0e-10,
        1.0e-8,
    )

    ix = np.rint(
        x / coordinate_tolerance
    ).astype(np.int64)

    iy = np.rint(
        y / coordinate_tolerance
    ).astype(np.int64)

    coordinate_keys = np.column_stack(
        (ix, iy)
    )

    _, inverse = np.unique(
        coordinate_keys,
        axis=0,
        return_inverse=True,
    )

    counts = np.bincount(
        inverse
    ).astype(float)

    x_merged = (
        np.bincount(
            inverse,
            weights=x,
        )
        / counts
    )

    y_merged = (
        np.bincount(
            inverse,
            weights=y,
        )
        / counts
    )

    log_rho_merged = (
        np.bincount(
            inverse,
            weights=np.log10(
                np.maximum(
                    rho,
                    1.0e-300,
                )
            ),
        )
        / counts
    )

    rho_merged = (
        10.0**log_rho_merged
    )

    vx_merged = (
        np.bincount(
            inverse,
            weights=vx,
        )
        / counts
    )

    vy_merged = (
        np.bincount(
            inverse,
            weights=vy,
        )
        / counts
    )

    return {
        "x": x_merged,
        "y": y_merged,
        "rho": rho_merged,
        "vx": vx_merged,
        "vy": vy_merged,
    }


# ============================================================
# 14. 1スナップショット分のx-y面を抽出
# ============================================================
def extract_snapshot(step_info):
    collected = {
        "x": [],
        "y": [],
        "rho": [],
        "vx": [],
        "vy": [],
    }

    used_blocks = 0

    for filename in step_info["files"]:
        try:
            block_data = (
                extract_xy_midplane(
                    filename
                )
            )

            if block_data is None:
                continue

            used_blocks += 1

            for key in collected:
                collected[key].append(
                    block_data[key]
                )

            del block_data

        except Exception as error:
            print(
                f"[WARNING] "
                f"{Path(filename).name}: "
                f"{error}"
            )

    if not collected["x"]:
        raise RuntimeError(
            f"No x-y data at "
            f"step={step_info['step']:05d}"
        )

    raw_data = {
        key: np.concatenate(values)
        for key, values in collected.items()
    }

    merged = merge_duplicate_points(
        raw_data
    )

    merged["step"] = (
        step_info["step"]
    )
    merged["time_code"] = (
        step_info["time_code"]
    )
    merged["time_kyr"] = (
        step_info["time_kyr"]
    )
    merged["time_kyr_integer"] = (
        step_info["time_kyr_integer"]
    )

    print(
        f"[INFO] step="
        f"{step_info['step']:05d}: "
        f"blocks={used_blocks}/"
        f"{len(step_info['files'])}, "
        f"raw={len(raw_data['x']):,}, "
        f"merged={len(merged['x']):,}"
    )

    del raw_data
    del collected

    return merged


# ============================================================
# 15. 線形補間
#
# nearest補完は行わない。
# 実データの凸包外側はNaNのまま残す。
# ============================================================
def interpolate_field_linear(
    points_2d,
    values,
    X,
    Y,
):
    return griddata(
        points_2d,
        np.asarray(values),
        (X, Y),
        method="linear",
        fill_value=np.nan,
    )


# ============================================================
# 16. 全時刻共通の密度カラースケールを計算
#
# 全スナップショットをメモリに保持せず、
# 1時刻ずつ読み込んで密度サンプルだけ保存する。
# ============================================================
density_samples = []

print(
    "\n[INFO] Pass 1/2: "
    "scanning global density range"
)

for index, step_info in enumerate(
    step_table
):
    print(
        f"[INFO] Scan "
        f"{index + 1}/"
        f"{len(step_table)}: "
        f"step={step_info['step']:05d}"
    )

    try:
        snapshot = extract_snapshot(
            step_info
        )

        radius_au = np.hypot(
            snapshot["x"],
            snapshot["y"],
        )

        density_mask = (
            np.isfinite(
                snapshot["rho"]
            )
            & (
                snapshot["rho"]
                > 0.0
            )
            & (
                radius_au
                >= SINK_RADIUS_AU
            )
        )

        values = snapshot["rho"][
            density_mask
        ]

        if values.size > 0:
            if (
                values.size
                > MAX_DENSITY_SAMPLES_PER_STEP
            ):
                sample_indices = np.linspace(
                    0,
                    values.size - 1,
                    MAX_DENSITY_SAMPLES_PER_STEP,
                    dtype=int,
                )

                values = values[
                    sample_indices
                ]

            density_samples.append(
                values.copy()
            )

        del values
        del density_mask
        del radius_au
        del snapshot

        gc.collect()

    except Exception as error:
        print(
            f"[WARNING] Density scan failed, "
            f"step={step_info['step']:05d}: "
            f"{error}"
        )

if not density_samples:
    raise RuntimeError(
        "No positive density values were found "
        "outside the sink."
    )

global_density = np.concatenate(
    density_samples
)

rho_vmin, rho_vmax = np.percentile(
    global_density,
    DENSITY_PERCENTILES,
)

del global_density
del density_samples
gc.collect()

if (
    not np.isfinite(rho_vmin)
    or not np.isfinite(rho_vmax)
    or rho_vmin <= 0.0
    or rho_vmax <= rho_vmin
):
    raise RuntimeError(
        "Invalid global density range: "
        f"{rho_vmin}, {rho_vmax}"
    )

density_norm = LogNorm(
    vmin=rho_vmin,
    vmax=rho_vmax,
    clip=True,
)

print(
    f"\n[INFO] Global density range: "
    f"{rho_vmin:.6e} -- "
    f"{rho_vmax:.6e} "
    "M_sun AU^-3"
)


# ============================================================
# 17. 実データ領域に対応する補間グリッド
# ============================================================
x_axis_au = np.linspace(
    X_MIN_AU,
    X_MAX_AU,
    MAP_RESOLUTION,
)

y_axis_au = np.linspace(
    Y_MIN_AU,
    Y_MAX_AU,
    MAP_RESOLUTION,
)

X_map_au, Y_map_au = np.meshgrid(
    x_axis_au,
    y_axis_au,
)

# 描画座標を10^5 AU単位へ変換
X_map_plot = (
    X_map_au
    / AXIS_SCALE_AU
)

Y_map_plot = (
    Y_map_au
    / AXIS_SCALE_AU
)

X_MIN_PLOT = (
    X_MIN_AU
    / AXIS_SCALE_AU
)

X_MAX_PLOT = (
    X_MAX_AU
    / AXIS_SCALE_AU
)

Y_MIN_PLOT = (
    Y_MIN_AU
    / AXIS_SCALE_AU
)

Y_MAX_PLOT = (
    Y_MAX_AU
    / AXIS_SCALE_AU
)

# NaN領域を白抜き表示
density_cmap = plt.get_cmap(
    DENSITY_CMAP_NAME
).copy()

density_cmap.set_bad(
    color="white"
)


# ============================================================
# 18. 全時刻を1枚ずつ描画・保存
# ============================================================
saved_files = []

print(
    "\n[INFO] Pass 2/2: "
    "plotting and saving PNG files"
)

for frame_index, step_info in enumerate(
    step_table
):
    print(
        f"[INFO] Plot "
        f"{frame_index + 1}/"
        f"{len(step_table)}: "
        f"step={step_info['step']:05d}"
    )

    try:
        snapshot = extract_snapshot(
            step_info
        )

        points_xy = np.column_stack(
            (
                snapshot["x"],
                snapshot["y"],
            )
        )

        # ----------------------------------------------------
        # 密度はlog10(rho)を線形補間
        # ----------------------------------------------------
        log_rho_map = (
            interpolate_field_linear(
                points_xy,
                np.log10(
                    np.maximum(
                        snapshot["rho"],
                        1.0e-300,
                    )
                ),
                X_map_au,
                Y_map_au,
            )
        )

        rho_map = np.full_like(
            log_rho_map,
            np.nan,
        )

        finite_log_rho = np.isfinite(
            log_rho_map
        )

        rho_map[finite_log_rho] = (
            10.0
            ** log_rho_map[
                finite_log_rho
            ]
        )

        # ----------------------------------------------------
        # 速度を線形補間
        # ----------------------------------------------------
        vx_map = (
            interpolate_field_linear(
                points_xy,
                snapshot["vx"],
                X_map_au,
                Y_map_au,
            )
        )

        vy_map = (
            interpolate_field_linear(
                points_xy,
                snapshot["vy"],
                X_map_au,
                Y_map_au,
            )
        )

        # ----------------------------------------------------
        # 補間後に速度を規格化
        # ----------------------------------------------------
        speed_map = np.hypot(
            vx_map,
            vy_map,
        )

        finite_speed = speed_map[
            np.isfinite(speed_map)
            & (speed_map > 0.0)
        ]

        if finite_speed.size > 0:
            minimum_speed = (
                np.median(finite_speed)
                * MIN_SPEED_FRACTION
            )
        else:
            minimum_speed = np.inf

        valid_velocity = (
            np.isfinite(vx_map)
            & np.isfinite(vy_map)
            & np.isfinite(speed_map)
            & (
                speed_map
                > minimum_speed
            )
        )

        unit_vx = np.full_like(
            vx_map,
            np.nan,
        )

        unit_vy = np.full_like(
            vy_map,
            np.nan,
        )

        unit_vx[valid_velocity] = (
            vx_map[valid_velocity]
            / speed_map[valid_velocity]
        )

        unit_vy[valid_velocity] = (
            vy_map[valid_velocity]
            / speed_map[valid_velocity]
        )

        # ----------------------------------------------------
        # 矢印を間引く
        # ----------------------------------------------------
        quiver_indices = np.linspace(
            0,
            MAP_RESOLUTION - 1,
            QUIVER_N,
            dtype=int,
        )

        quiver_selection = np.ix_(
            quiver_indices,
            quiver_indices,
        )

        X_quiver = (
            X_map_plot[
                quiver_selection
            ]
        )

        Y_quiver = (
            Y_map_plot[
                quiver_selection
            ]
        )

        U_quiver = (
            unit_vx[
                quiver_selection
            ]
        )

        V_quiver = (
            unit_vy[
                quiver_selection
            ]
        )

        valid_quiver = (
            np.isfinite(U_quiver)
            & np.isfinite(V_quiver)
        )

        # 描画座標単位での矢印長
        arrow_length_plot = (
            1.25
            * min(
                X_MAX_PLOT
                - X_MIN_PLOT,
                Y_MAX_PLOT
                - Y_MIN_PLOT,
            )
            / max(
                QUIVER_N - 1,
                1,
            )
        )

        # ----------------------------------------------------
        # 描画
        # ----------------------------------------------------
        fig, ax = plt.subplots(
            figsize=FIGSIZE
        )

        rho_masked = np.ma.masked_invalid(
            rho_map
        )

        ax.pcolormesh(
            X_map_plot,
            Y_map_plot,
            rho_masked,
            shading="auto",
            cmap=density_cmap,
            norm=density_norm,
            rasterized=True,
        )

        ax.quiver(
            X_quiver[valid_quiver],
            Y_quiver[valid_quiver],
            (
                U_quiver[valid_quiver]
                * arrow_length_plot
            ),
            (
                V_quiver[valid_quiver]
                * arrow_length_plot
            ),
            color=QUIVER_COLOR,
            angles="xy",
            scale_units="xy",
            scale=1.0,
            width=0.0025,
            pivot="mid",
            headwidth=3.5,
            headlength=4.5,
            headaxislength=4.0,
            alpha=0.9,
        )

        ax.axhline(
            0.0,
            color="gray",
            linewidth=0.7,
            alpha=0.7,
        )

        ax.axvline(
            0.0,
            color="gray",
            linewidth=0.7,
            alpha=0.7,
        )

        ax.plot(
            0.0,
            0.0,
            marker="+",
            color="red",
            markersize=12,
            markeredgewidth=2.5,
        )

        ax.set_xlim(
            X_MIN_PLOT,
            X_MAX_PLOT,
        )

        ax.set_ylim(
            Y_MIN_PLOT,
            Y_MAX_PLOT,
        )

        # 座標軸名と単位を明記
        ax.set_xlabel(
            r"$x\ [10^5\ {\rm AU}]$",
            fontsize=13,
        )

        ax.set_ylabel(
            r"$y\ [10^5\ {\rm AU}]$",
            fontsize=13,
        )

        ax.set_title(
            "Density + normalized "
            r"$(v_x,v_y)$ direction"
            "\n"
            f"step={snapshot['step']:05d}, "
            f"t={snapshot['time_kyr_integer']} kyr",
            fontsize=14,
        )

        ax.set_aspect(
            "equal",
            adjustable="box",
        )

        ax.grid(False)

        fig.tight_layout()

        # frame番号、VTK step、整数kyrの順で保存
        output_file = output_dir / (
            f"frame_{frame_index:05d}_"
            f"step_{snapshot['step']:05d}_"
            f"time_{snapshot['time_kyr_integer']:08d}kyr.png"
        )

        fig.savefig(
            output_file,
            dpi=DPI,
            bbox_inches="tight",
            facecolor="white",
        )

        plt.close(fig)

        saved_files.append(
            output_file
        )

        print(
            f"[INFO] Saved: "
            f"{output_file.name}"
        )

        # ----------------------------------------------------
        # メモリ解放
        # ----------------------------------------------------
        del snapshot
        del points_xy
        del log_rho_map
        del rho_map
        del finite_log_rho
        del vx_map
        del vy_map
        del speed_map
        del finite_speed
        del valid_velocity
        del unit_vx
        del unit_vy
        del rho_masked
        del fig
        del ax

        gc.collect()

    except Exception as error:
        plt.close("all")
        gc.collect()

        print(
            f"[WARNING] Plot failed, "
            f"step={step_info['step']:05d}: "
            f"{error}"
        )


# ============================================================
# 19. 終了表示
# ============================================================
print(
    f"\n[INFO] Completed: "
    f"{len(saved_files)} PNG files"
)

print(
    f"[INFO] Output directory: "
    f"{output_dir}"
)

print(
    f"[INFO] Existing PNG count: "
    f"{len(list(output_dir.glob('*.png')))}"
)

if saved_files:
    print(
        f"[INFO] First output: "
        f"{saved_files[0]}"
    )
    print(
        f"[INFO] Last output : "
        f"{saved_files[-1]}"
    )


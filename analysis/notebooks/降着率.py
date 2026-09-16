# -*- coding: utf-8 -*-
# Converted from 降着率.ipynb

# %% cell 1
# ============================================================
# HSTからシンク表面の降着率 Mdot_flux を読み、
# 時間 [kyr]－降着率 [M_sun yr^-1] を描画する
# ============================================================

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# 入出力設定
# ============================================================
hst_file = Path(
    os.path.expanduser(
        "~/athena-project/results/〇〇/Toyouchi.hst"
    )
).resolve()

output_dir = Path(
    "./stellar_mass_rate_history"
).resolve()

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

output_file = (
    output_dir
    / "mdot_flux_vs_time.png"
)

cleaned_file = (
    output_dir
    / "mdot_flux_cleaned.txt"
)

if not hst_file.is_file():
    raise FileNotFoundError(
        f"HST file does not exist: {hst_file}"
    )


# ============================================================
# コード単位
# ============================================================
M_UNIT_CGS = 4.0e33
T_UNIT_CGS = 3.61e10

M_SUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0
KYR_CGS = 1.0e3 * YEAR_CGS


# code time -> kyr
CODE_TIME_TO_KYR = (
    T_UNIT_CGS
    / KYR_CGS
)

# code mass / code time -> M_sun / yr
CODE_MDOT_TO_MSUN_PER_YR = (
    M_UNIT_CGS
    / M_SUN_CGS
) / (
    T_UNIT_CGS
    / YEAR_CGS
)

print(
    f"[INFO] 1 code time = "
    f"{CODE_TIME_TO_KYR:.6e} kyr"
)

print(
    f"[INFO] 1 code mass-rate = "
    f"{CODE_MDOT_TO_MSUN_PER_YR:.6e} "
    f"M_sun yr^-1"
)


# ============================================================
# HSTヘッダーから列番号を取得
#
# Athena++:
#   [1]=time [2]=dt ...
#
# NumPy:
#   0始まりなので [n] -> n-1
# ============================================================
def read_hst_column_map(filename):
    column_map = {}

    with open(
        filename,
        "r",
        encoding="utf-8",
        errors="ignore",
    ) as handle:

        for line in handle:
            if not line.lstrip().startswith("#"):
                continue

            matches = re.findall(
                r"\[(\d+)\]\s*=\s*([^\s]+)",
                line,
            )

            for number_text, name in matches:
                hst_number = int(
                    number_text
                )

                clean_name = (
                    name.strip()
                    .strip(",")
                )

                column_map[clean_name] = (
                    hst_number - 1
                )

    return column_map


column_map = read_hst_column_map(
    hst_file
)

if not column_map:
    raise RuntimeError(
        "HST header columns could not be parsed."
    )


# ============================================================
# 検出列を表示
# ============================================================
print(
    "\n[INFO] Detected HST columns:"
)

for name, numpy_index in sorted(
    column_map.items(),
    key=lambda item: item[1],
):
    print(
        f"       HST [{numpy_index + 1:2d}] "
        f"-> NumPy [{numpy_index:2d}] "
        f"= {name}"
    )


# ============================================================
# time / Mdot_flux 列の確認
# ============================================================
if "time" not in column_map:
    raise KeyError(
        "HST header does not contain 'time'."
    )

if "Mdot_flux" not in column_map:
    raise KeyError(
        "HST header does not contain "
        "'Mdot_flux'."
    )

time_column = (
    column_map["time"]
)

mdot_column = (
    column_map["Mdot_flux"]
)

print(
    f"\n[INFO] time      : "
    f"HST [{time_column + 1}] "
    f"-> NumPy [{time_column}]"
)

print(
    f"[INFO] Mdot_flux : "
    f"HST [{mdot_column + 1}] "
    f"-> NumPy [{mdot_column}]"
)


# ============================================================
# HSTデータ読み込み
# ============================================================
data = np.loadtxt(
    hst_file,
    comments="#",
    ndmin=2,
)

required_column = max(
    time_column,
    mdot_column,
)

if data.shape[1] <= required_column:
    raise ValueError(
        f"HST data has only "
        f"{data.shape[1]} columns, "
        f"but index {required_column} "
        f"is required."
    )

time_code_raw = np.asarray(
    data[:, time_column],
    dtype=float,
)

mdot_code_raw = np.asarray(
    data[:, mdot_column],
    dtype=float,
)


# ============================================================
# NaN / inf を除外
# ============================================================
finite_mask = (
    np.isfinite(time_code_raw)
    & np.isfinite(mdot_code_raw)
)

time_code_raw = (
    time_code_raw[finite_mask]
)

mdot_code_raw = (
    mdot_code_raw[finite_mask]
)

if time_code_raw.size == 0:
    raise RuntimeError(
        "No finite time/Mdot_flux "
        "data were found."
    )


# ============================================================
# restart重複の整理
#
# 後から書かれたrestart系列を優先する。
# ============================================================
def keep_latest_restart_branch(
    time_values,
    *value_arrays,
):

    time_values = np.asarray(
        time_values,
        dtype=float,
    )

    keep_reversed = []
    latest_accepted_time = np.inf

    for index in range(
        len(time_values) - 1,
        -1,
        -1,
    ):

        current_time = (
            time_values[index]
        )

        scale = max(
            1.0,
            abs(current_time),
            (
                abs(latest_accepted_time)
                if np.isfinite(
                    latest_accepted_time
                )
                else 1.0
            ),
        )

        tolerance = (
            64.0
            * np.finfo(float).eps
            * scale
        )

        if (
            current_time
            < latest_accepted_time
            - tolerance
        ):
            keep_reversed.append(
                index
            )

            latest_accepted_time = (
                current_time
            )

        elif np.isclose(
            current_time,
            latest_accepted_time,
            rtol=0.0,
            atol=tolerance,
        ):
            # 同一時刻なら後に書かれた
            # データを優先
            continue

    keep_indices = np.asarray(
        keep_reversed[::-1],
        dtype=int,
    )

    cleaned_time = (
        time_values[keep_indices]
    )

    cleaned_values = [
        np.asarray(values)[
            keep_indices
        ]
        for values in value_arrays
    ]

    return (
        cleaned_time,
        cleaned_values,
        keep_indices,
    )


(
    time_code,
    cleaned_arrays,
    keep_indices,
) = keep_latest_restart_branch(
    time_code_raw,
    mdot_code_raw,
)

mdot_code = (
    cleaned_arrays[0]
)


print(
    f"\n[INFO] Raw rows     : "
    f"{len(time_code_raw)}"
)

print(
    f"[INFO] Cleaned rows : "
    f"{len(time_code)}"
)

print(
    f"[INFO] Removed rows : "
    f"{len(time_code_raw) - len(time_code)}"
)

if len(time_code) < 2:
    raise RuntimeError(
        "Fewer than two rows remain "
        "after restart cleaning."
    )


# ============================================================
# 物理単位へ変換
# ============================================================
time_kyr = (
    time_code
    * CODE_TIME_TO_KYR
)

mdot_msun_per_yr = (
    mdot_code
    * CODE_MDOT_TO_MSUN_PER_YR
)


# ============================================================
# 平均降着率
# ============================================================

# 単純平均
mdot_sample_mean = np.mean(
    mdot_msun_per_yr
)

# 時間重み付き平均
time_yr = (
    time_kyr
    * 1.0e3
)

elapsed_time_yr = (
    time_yr[-1]
    - time_yr[0]
)

if elapsed_time_yr > 0.0:

    if hasattr(
        np,
        "trapezoid",
    ):
        integrated_mass_msun = (
            np.trapezoid(
                mdot_msun_per_yr,
                x=time_yr,
            )
        )

    else:
        integrated_mass_msun = (
            np.trapz(
                mdot_msun_per_yr,
                x=time_yr,
            )
        )

    mdot_time_mean = (
        integrated_mass_msun
        / elapsed_time_yr
    )

else:
    integrated_mass_msun = np.nan
    mdot_time_mean = np.nan


# ============================================================
# 診断値
# ============================================================
mdot_min = np.min(
    mdot_msun_per_yr
)

mdot_max = np.max(
    mdot_msun_per_yr
)

negative_count = np.count_nonzero(
    mdot_msun_per_yr < 0.0
)


print(
    f"\n[INFO] Time range = "
    f"{time_kyr[0]:.6e} -- "
    f"{time_kyr[-1]:.6e} kyr"
)

print(
    f"[INFO] Mdot range = "
    f"{mdot_min:.6e} -- "
    f"{mdot_max:.6e} "
    f"M_sun yr^-1"
)

print(
    f"[INFO] Time-weighted mean Mdot = "
    f"{mdot_time_mean:.6e} "
    f"M_sun yr^-1"
)

print(
    f"[INFO] Sample mean Mdot        = "
    f"{mdot_sample_mean:.6e} "
    f"M_sun yr^-1"
)

print(
    f"[INFO] Integral of Mdot        = "
    f"{integrated_mass_msun:.6e} "
    f"M_sun"
)

print(
    f"[INFO] Negative Mdot rows      = "
    f"{negative_count}/"
    f"{len(mdot_msun_per_yr)}"
)


# ============================================================
# 整理後データ保存
# ============================================================
np.savetxt(
    cleaned_file,
    np.column_stack(
        (
            time_kyr,
            mdot_msun_per_yr,
        )
    ),
    header=(
        "time_kyr "
        "Mdot_flux_Msun_per_yr"
    ),
)


# ============================================================
# プロット
# ============================================================
fig, ax = plt.subplots(
    figsize=(9, 6),
)

ax.plot(
    time_kyr,
    mdot_msun_per_yr,
    color="tab:blue",
    linewidth=2.0,
    label=r"$\dot{M}_{\rm flux}$",
)

# 時間重み付き平均
ax.axhline(
    mdot_time_mean,
    color="tab:red",
    linewidth=1.5,
    linestyle="--",
    label="Time-weighted mean",
)

# Mdot = 0
ax.axhline(
    0.0,
    color="black",
    linewidth=0.8,
    alpha=0.7,
)

ax.set_xlabel(
    "Time [kyr]"
)

ax.set_ylabel(
    r"$\dot{M}_{\rm flux}$ "
    r"$[M_\odot\,{\rm yr}^{-1}]$"
)

ax.set_title(
    "Accretion rate through the sink surface"
)

ax.grid(
    True,
    linestyle="--",
    alpha=0.45,
)

ax.legend(
    loc="best",
)

# 必要なら有効化
# ax.set_yscale(
#     "symlog",
#     linthresh=1.0e-6,
# )

fig.tight_layout()

fig.savefig(
    output_file,
    dpi=200,
    bbox_inches="tight",
    facecolor="white",
)

plt.show()

plt.close(fig)


# ============================================================
# 保存先
# ============================================================
print(
    f"\n[SAVED] {output_file}"
)

print(
    f"[SAVED] {cleaned_file}"
)


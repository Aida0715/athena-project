# -*- coding: utf-8 -*-
# Converted from 中心星質量変化.ipynb

# %% cell 1
# ============================================================
# 中心星質量 Mstar の時間変化
#
# Athena++ HSTヘッダー：1始まり
# NumPy配列           ：0始まり
#
# ヘッダーの [n] をNumPy列 n-1 へ変換する
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
    "./stellar_mass_history"
).resolve()

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

output_file = (
    output_dir
    / "stellar_mass_vs_time.png"
)

cleaned_file = (
    output_dir
    / "stellar_mass_cleaned.txt"
)

if not hst_file.is_file():
    raise FileNotFoundError(
        f"HST file does not exist: {hst_file}"
    )


# ============================================================
# 単位系
# ============================================================
M_UNIT_CGS = 4.0e33
T_UNIT_CGS = 3.61e10

MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0
KYR_CGS = 1.0e3 * YEAR_CGS

MASS_UNIT_MSUN = (
    M_UNIT_CGS
    / MSUN_CGS
)

TIME_UNIT_KYR = (
    T_UNIT_CGS
    / KYR_CGS
)

print(
    f"[INFO] 1 code mass = "
    f"{MASS_UNIT_MSUN:.6e} M_sun"
)
print(
    f"[INFO] 1 code time = "
    f"{TIME_UNIT_KYR:.6e} kyr"
)


# ============================================================
# HSTヘッダーを解析
#
# 重要：
# HSTの [1], [2], ... は1始まり。
# NumPy列として使うときは必ず1を引く。
# ============================================================
def read_hst_column_map(filename):
    column_map = {}

    with open(
        filename,
        "r",
        encoding="utf-8",
        errors="ignore",
    ) as file:
        for line in file:
            if not line.lstrip().startswith("#"):
                continue

            matches = re.findall(
                r"\[(\d+)\]\s*=\s*([^\s]+)",
                line,
            )

            for hst_index_text, name in matches:
                hst_index = int(
                    hst_index_text
                )

                numpy_index = (
                    hst_index - 1
                )

                clean_name = (
                    name.strip().strip(",")
                )

                column_map[clean_name] = {
                    "hst_index": hst_index,
                    "numpy_index": numpy_index,
                }

    return column_map


column_map = read_hst_column_map(
    hst_file
)

if not column_map:
    raise RuntimeError(
        "HST header columns could not be parsed."
    )

print("[INFO] Detected HST columns:")

for name, information in sorted(
    column_map.items(),
    key=lambda item: item[1]["hst_index"],
):
    print(
        f"       HST [{information['hst_index']:2d}] "
        f"-> NumPy [{information['numpy_index']:2d}] "
        f"= {name}"
    )


# ============================================================
# 名前から列を検索
# ============================================================
def find_column(column_map, candidates):
    lower_map = {
        name.lower(): information
        for name, information
        in column_map.items()
    }

    for candidate in candidates:
        key = candidate.lower()

        if key in lower_map:
            return lower_map[key]

    return None


time_information = find_column(
    column_map,
    ("time",),
)

mstar_information = find_column(
    column_map,
    (
        "Mstar",
        "mstar",
        "M_star",
    ),
)

if time_information is None:
    raise KeyError(
        "time column was not found."
    )

if mstar_information is None:
    raise KeyError(
        "Mstar column was not found."
    )

time_column = (
    time_information["numpy_index"]
)

mstar_column = (
    mstar_information["numpy_index"]
)

print(
    f"[INFO] time : "
    f"HST [{time_information['hst_index']}] "
    f"-> NumPy [{time_column}]"
)

print(
    f"[INFO] Mstar: "
    f"HST [{mstar_information['hst_index']}] "
    f"-> NumPy [{mstar_column}]"
)


# ============================================================
# データ読み込み
# ============================================================
data = np.loadtxt(
    hst_file,
    comments="#",
)

data = np.atleast_2d(data)

if data.shape[1] <= max(
    time_column,
    mstar_column,
):
    raise ValueError(
        f"HST data has only "
        f"{data.shape[1]} columns."
    )

time_code_raw = data[:, time_column]
mstar_code_raw = data[:, mstar_column]

finite = (
    np.isfinite(time_code_raw)
    & np.isfinite(mstar_code_raw)
)

time_code_raw = time_code_raw[finite]
mstar_code_raw = mstar_code_raw[finite]

if time_code_raw.size == 0:
    raise RuntimeError(
        "No finite Mstar data were found."
    )


# ============================================================
# restartにより重複した時間系列を整理
#
# 時刻が巻き戻った場合、その時刻以降の古い系列を
# 新しいrestart系列で置き換える。
# ============================================================
clean_time = []
clean_mstar = []

restart_count = 0
removed_row_count = 0

for time_value, mass_value in zip(
    time_code_raw,
    mstar_code_raw,
):
    tolerance = max(
        1.0e-12,
        abs(time_value) * 1.0e-12,
    )

    if (
        clean_time
        and time_value
        <= clean_time[-1] + tolerance
    ):
        restart_count += 1

        while (
            clean_time
            and clean_time[-1]
            >= time_value - tolerance
        ):
            clean_time.pop()
            clean_mstar.pop()
            removed_row_count += 1

    clean_time.append(
        time_value
    )

    clean_mstar.append(
        mass_value
    )

time_code = np.asarray(
    clean_time,
    dtype=float,
)

mstar_code = np.asarray(
    clean_mstar,
    dtype=float,
)

print(
    f"[INFO] Raw rows       : "
    f"{time_code_raw.size}"
)
print(
    f"[INFO] Cleaned rows   : "
    f"{time_code.size}"
)
print(
    f"[INFO] Restart events : "
    f"{restart_count}"
)
print(
    f"[INFO] Removed rows   : "
    f"{removed_row_count}"
)


# ============================================================
# 物理単位へ変換
# ============================================================
time_kyr = (
    time_code
    * TIME_UNIT_KYR
)

mstar_msun = (
    mstar_code
    * MASS_UNIT_MSUN
)


# ============================================================
# Mstarが減少していないか確認
# ============================================================
mass_difference = np.diff(
    mstar_msun
)

mass_tolerance = max(
    1.0e-12,
    np.max(
        np.abs(mstar_msun)
    ) * 1.0e-10,
)

decrease_indices = np.where(
    mass_difference
    < -mass_tolerance
)[0]

if decrease_indices.size == 0:
    print(
        "[OK] Mstar is monotonically "
        "non-decreasing."
    )

else:
    print(
        f"[WARNING] Mstar decreases at "
        f"{decrease_indices.size} locations."
    )

    for index in decrease_indices[:20]:
        print(
            f"  t={time_kyr[index]:.6e}"
            f" -> {time_kyr[index + 1]:.6e} kyr, "
            f"Mstar={mstar_msun[index]:.6e}"
            f" -> {mstar_msun[index + 1]:.6e} M_sun"
        )


# ============================================================
# 整理後データを保存
# ============================================================
np.savetxt(
    cleaned_file,
    np.column_stack(
        (
            time_kyr,
            mstar_msun,
        )
    ),
    header="time_kyr Mstar_Msun",
)

print(
    f"[INFO] Initial Mstar = "
    f"{mstar_msun[0]:.6e} M_sun"
)
print(
    f"[INFO] Final Mstar   = "
    f"{mstar_msun[-1]:.6e} M_sun"
)


# ============================================================
# 描画
# ============================================================
fig, ax = plt.subplots(
    figsize=(7, 5),
)

ax.plot(
    time_kyr,
    mstar_msun,
    color="black",
    linewidth=1.8,
)

if decrease_indices.size > 0:
    warning_points = (
        decrease_indices + 1
    )

    ax.scatter(
        time_kyr[warning_points],
        mstar_msun[warning_points],
        color="red",
        s=25,
        zorder=5,
        label="Mstar decrease",
    )

    ax.legend()

ax.set_xlabel(
    "Time [kyr]"
)

ax.set_ylabel(
    r""
)

ax.set_title(
    "Central stellar mass evolution"
)

ax.grid(
    True,
    linestyle="--",
    alpha=0.4,
)

# 1点しか残らなかった場合でも点を表示
if time_kyr.size == 1:
    ax.scatter(
        time_kyr,
        mstar_msun,
        color="black",
        s=40,
        zorder=5,
    )

fig.tight_layout()

fig.savefig(
    output_file,
    dpi=200,
    bbox_inches="tight",
    facecolor="white",
)

plt.show()
plt.close(fig)

print(f"[INFO] Saved: {output_file}")
print(f"[INFO] Saved: {cleaned_file}")
     


# -*- coding: utf-8 -*-
"""複数計算の中心星質量 Mstar を1枚の図で比較する。"""

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# 1. ユーザー設定
#
# 2本目以降は、次の2行を番号を増やして追記するだけでよい。
# datadir2 = "~/athena-project/results/○○"
# curve_name2 = r"$B_{z,0}=1\times10^{-19}\ \mu{\rm G}$"
# ============================================================
datadir1 = "~/athena-project/results/○○"
curve_name1 = r"$B_{z0}=3\times10^{-20}\ \mu{\rm G}$"


# 以下は通常変更不要
hst_filename = "Toyouchi.hst"
output_dir = Path("./stellar_mass_history").resolve()
output_file = output_dir / "stellar_mass_vs_time.png"
output_dir.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. 単位系
# ============================================================
M_UNIT_CGS = 4.0e33
T_UNIT_CGS = 3.61e10
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0

MASS_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS
TIME_UNIT_KYR = T_UNIT_CGS / (1.0e3 * YEAR_CGS)

print(f"[INFO] 1 code mass = {MASS_UNIT_MSUN:.6e} M_sun")
print(f"[INFO] 1 code time = {TIME_UNIT_KYR:.6e} kyr")


# ============================================================
# 3. datadir1, datadir2, ... と対応するcurve_nameを自動取得
# ============================================================
def collect_manual_datasets(namespace):
    numbers = sorted(
        int(match.group(1))
        for name in namespace
        if (match := re.fullmatch(r"datadir(\d+)", name))
    )
    if not numbers:
        raise RuntimeError("datadir1, datadir2, ... が設定されていません。")

    datasets = []
    for number in numbers:
        directory_key = f"datadir{number}"
        name_key = f"curve_name{number}"
        if name_key not in namespace:
            raise NameError(
                f"{directory_key} に対応する {name_key} がありません。"
            )

        directory = Path(
            os.path.expanduser(str(namespace[directory_key]))
        ).resolve()
        hst_file = directory if directory.is_file() else directory / hst_filename
        if not hst_file.is_file():
            raise FileNotFoundError(
                f"[{directory_key}] HST file does not exist: {hst_file}"
            )

        datasets.append({
            "number": number,
            "datadir": directory,
            "hst_file": hst_file,
            "curve_name": str(namespace[name_key]),
        })
    return datasets


datasets = collect_manual_datasets(globals())
print("\n[INFO] Datasets:")
for dataset in datasets:
    print(
        f"  [{dataset['number']}] {dataset['curve_name']}\n"
        f"      {dataset['hst_file']}"
    )


# ============================================================
# 4. HSTヘッダー解析
# Athena++ HSTヘッダーは1始まり、NumPy配列は0始まり。
# ============================================================
def read_hst_column_map(filename):
    column_map = {}
    with open(filename, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            if not line.lstrip().startswith("#"):
                continue
            for index_text, name in re.findall(
                r"\[(\d+)\]\s*=\s*([^\s]+)", line
            ):
                hst_index = int(index_text)
                column_map[name.strip().strip(",")] = {
                    "hst_index": hst_index,
                    "numpy_index": hst_index - 1,
                }
    return column_map


def find_column(column_map, candidates):
    lower_map = {name.lower(): info for name, info in column_map.items()}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


# ============================================================
# 5. restartで重複した時間系列を、新しい系列で置換
# ============================================================
def clean_restart_history(time_raw, value_raw):
    clean_time, clean_value = [], []
    restart_count = removed_count = 0

    for time_value, data_value in zip(time_raw, value_raw):
        tolerance = max(1.0e-12, abs(time_value) * 1.0e-12)
        if clean_time and time_value <= clean_time[-1] + tolerance:
            restart_count += 1
            while clean_time and clean_time[-1] >= time_value - tolerance:
                clean_time.pop()
                clean_value.pop()
                removed_count += 1
        clean_time.append(time_value)
        clean_value.append(data_value)

    return (
        np.asarray(clean_time, dtype=float),
        np.asarray(clean_value, dtype=float),
        restart_count,
        removed_count,
    )


# ============================================================
# 6. 1計算分のMstarを読み込む
# ============================================================
def load_stellar_mass_history(dataset):
    hst_file = dataset["hst_file"]
    column_map = read_hst_column_map(hst_file)
    if not column_map:
        raise RuntimeError(f"HST header could not be parsed: {hst_file}")

    time_info = find_column(column_map, ("time",))
    mstar_info = find_column(column_map, ("Mstar", "mstar", "M_star"))
    if time_info is None or mstar_info is None:
        raise KeyError(f"time or Mstar column was not found: {hst_file}")

    time_column = time_info["numpy_index"]
    mstar_column = mstar_info["numpy_index"]
    data = np.atleast_2d(np.loadtxt(hst_file, comments="#"))
    if data.shape[1] <= max(time_column, mstar_column):
        raise ValueError(f"HST has only {data.shape[1]} columns: {hst_file}")

    time_raw = data[:, time_column]
    mstar_raw = data[:, mstar_column]
    finite = np.isfinite(time_raw) & np.isfinite(mstar_raw)
    time_raw, mstar_raw = time_raw[finite], mstar_raw[finite]
    if time_raw.size == 0:
        raise RuntimeError(f"No finite Mstar data were found: {hst_file}")

    time_code, mstar_code, restart_count, removed_count = (
        clean_restart_history(time_raw, mstar_raw)
    )
    time_kyr = time_code * TIME_UNIT_KYR
    mstar_msun = mstar_code * MASS_UNIT_MSUN

    mass_tolerance = max(
        1.0e-12, np.max(np.abs(mstar_msun)) * 1.0e-10
    )
    decrease_indices = np.where(
        np.diff(mstar_msun) < -mass_tolerance
    )[0]

    print(
        f"\n[INFO] Dataset {dataset['number']}: {dataset['curve_name']}\n"
        f"       Raw/Cleaned rows: {time_raw.size}/{time_code.size}\n"
        f"       Restart events  : {restart_count}\n"
        f"       Removed rows    : {removed_count}\n"
        f"       Initial/Final   : {mstar_msun[0]:.6e} / "
        f"{mstar_msun[-1]:.6e} M_sun"
    )
    if decrease_indices.size:
        print(f"[WARNING] Mstar decreases at {decrease_indices.size} locations.")
    else:
        print("[OK] Mstar is monotonically non-decreasing.")

    cleaned_file = output_dir / (
        f"stellar_mass_cleaned_{dataset['number']:02d}.txt"
    )
    np.savetxt(
        cleaned_file,
        np.column_stack((time_kyr, mstar_msun)),
        header="time_kyr Mstar_Msun",
    )
    print(f"[INFO] Saved cleaned data: {cleaned_file}")

    return {
        **dataset,
        "time_kyr": time_kyr,
        "mstar_msun": mstar_msun,
        "decrease_indices": decrease_indices,
    }


histories = [load_stellar_mass_history(dataset) for dataset in datasets]


# ============================================================
# 7. 複数曲線を1枚に描画
# ============================================================
fig, ax = plt.subplots(figsize=(8.5, 6.0))

if len(histories) <= 10:
    colors = plt.get_cmap("tab10").colors[:len(histories)]
else:
    colors = plt.get_cmap("turbo")(
        np.linspace(0.05, 0.95, len(histories))
    )

for history, color in zip(histories, colors):
    time_kyr = history["time_kyr"]
    mstar_msun = history["mstar_msun"]
    ax.plot(
        time_kyr,
        mstar_msun,
        color=color,
        linewidth=1.9,
        label=history["curve_name"],
    )

    if time_kyr.size == 1:
        ax.scatter(time_kyr, mstar_msun, color=color, s=38, zorder=5)

    decrease_indices = history["decrease_indices"]
    if decrease_indices.size:
        warning_points = decrease_indices + 1
        ax.scatter(
            time_kyr[warning_points],
            mstar_msun[warning_points],
            color=color,
            marker="x",
            s=34,
            zorder=6,
            label="_nolegend_",
        )

ax.set_xlabel("Time [kyr]")
ax.set_ylabel(r"$M_{\rm star}\ [M_\odot]$")
ax.set_title("Central stellar mass evolution")
ax.grid(True, linestyle="--", alpha=0.35)
ax.margins(x=0.03, y=0.05)

# 情報ボックスをまず上部へ配置し、実際の描画サイズを測る。
# その高さに応じてy軸上側へ余白を自動追加するため、曲線とは重ならない。
information_box = ax.legend(
    loc="upper center",
    ncol=max(1, int(np.ceil(len(histories) / 6))),
    frameon=True,
    fancybox=True,
    framealpha=0.92,
    facecolor="white",
    edgecolor="0.35",
    title="Models",
)
information_box.set_zorder(20)

fig.tight_layout()
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
legend_bbox = information_box.get_window_extent(renderer=renderer)
axes_bbox = ax.get_window_extent(renderer=renderer)
reserved_fraction = legend_bbox.height / axes_bbox.height + 0.035

if reserved_fraction >= 0.85:
    raise RuntimeError(
        "The information box is too tall for the graph. "
        "Shorten curve names or reduce the number of curves."
    )

current_ymin, current_ymax = ax.get_ylim()
data_ymax = max(float(np.max(item["mstar_msun"])) for item in histories)
new_ymax = current_ymin + (
    (data_ymax - current_ymin) / (1.0 - reserved_fraction)
)
ax.set_ylim(current_ymin, max(current_ymax, new_ymax))

fig.savefig(
    output_file,
    dpi=200,
    bbox_inches="tight",
    facecolor="white",
)
plt.show()
plt.close(fig)

print(f"\n[INFO] Saved graph: {output_file}")

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


# %% cell 2
# ============================================================
# 複数計算の「時間重み付き平均降着率」だけを比較
#
# ここから下は、上の単一計算用コードとは独立して実行可能。
# 瞬間降着率曲線は描画せず、各計算の平均値の水平線だけを描画する。
# ============================================================

import os as _mdot_os
import re as _mdot_re
from pathlib import Path as _MdotPath

import matplotlib.pyplot as _mdot_plt
import numpy as _mdot_np


# ============================================================
# ユーザー設定
#
# 2本目以降は次の2行を番号を増やして追記するだけでよい。
# datadir2 = "~/athena-project/results/○○"
# curve_name2 = r"$B_{z,0}=1\times10^{-19}\ \mu{\rm G}$"
# ============================================================
datadir1 = "~/athena-project/results/○○"
curve_name1 = r"$B_{z0}=〇〇\times10^{〇〇}\ \mu{\rm G}$"


# 以下は通常変更不要
_MDOT_HST_FILENAME = "Toyouchi.hst"
_MDOT_OUTPUT_DIR = _MdotPath("./stellar_mass_rate_comparison").resolve()
_MDOT_OUTPUT_FILE = _MDOT_OUTPUT_DIR / "mdot_time_weighted_means.png"
_MDOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

_MDOT_M_UNIT_CGS = 4.0e33
_MDOT_T_UNIT_CGS = 3.61e10
_MDOT_MSUN_CGS = 1.98847e33
_MDOT_YEAR_CGS = 365.25 * 24.0 * 3600.0

_MDOT_CODE_TIME_TO_KYR = (
    _MDOT_T_UNIT_CGS / (1.0e3 * _MDOT_YEAR_CGS)
)
_MDOT_CODE_RATE_TO_MSUN_PER_YR = (
    (_MDOT_M_UNIT_CGS / _MDOT_MSUN_CGS)
    / (_MDOT_T_UNIT_CGS / _MDOT_YEAR_CGS)
)


def _mdot_collect_datasets(namespace):
    """datadir1, datadir2, ... と同番号のcurve_nameを収集。"""
    numbers = sorted(
        int(match.group(1))
        for name in namespace
        if (match := _mdot_re.fullmatch(r"datadir(\d+)", name))
    )
    if not numbers:
        raise RuntimeError("datadir1, datadir2, ... が設定されていません。")

    datasets = []
    for number in numbers:
        datadir_key = f"datadir{number}"
        curve_name_key = f"curve_name{number}"
        if curve_name_key not in namespace:
            raise NameError(
                f"{datadir_key} に対応する {curve_name_key} がありません。"
            )

        datadir = _MdotPath(
            _mdot_os.path.expanduser(str(namespace[datadir_key]))
        ).resolve()
        hst_file = (
            datadir if datadir.is_file()
            else datadir / _MDOT_HST_FILENAME
        )
        if not hst_file.is_file():
            raise FileNotFoundError(
                f"[{datadir_key}] HST file does not exist: {hst_file}"
            )

        datasets.append({
            "number": number,
            "hst_file": hst_file,
            "curve_name": str(namespace[curve_name_key]),
        })
    return datasets


def _mdot_read_column_map(filename):
    """Athena++ HSTヘッダーの1始まり列を0始まり列へ変換。"""
    column_map = {}
    with open(filename, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.lstrip().startswith("#"):
                continue
            for number_text, name in _mdot_re.findall(
                r"\[(\d+)\]\s*=\s*([^\s]+)", line
            ):
                column_map[name.strip().strip(",")] = int(number_text) - 1
    return column_map


def _mdot_find_column(column_map, candidates):
    lower_map = {name.lower(): index for name, index in column_map.items()}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


def _mdot_clean_restarts(time_values, rate_values):
    """restartで時刻が巻き戻った古い枝を、後から追記された枝で置換。"""
    clean_time = []
    clean_rate = []
    for time_value, rate_value in zip(time_values, rate_values):
        tolerance = max(1.0e-12, abs(time_value) * 1.0e-12)
        if clean_time and time_value <= clean_time[-1] + tolerance:
            while clean_time and clean_time[-1] >= time_value - tolerance:
                clean_time.pop()
                clean_rate.pop()
        clean_time.append(time_value)
        clean_rate.append(rate_value)

    return (
        _mdot_np.asarray(clean_time, dtype=float),
        _mdot_np.asarray(clean_rate, dtype=float),
    )


def _mdot_load_time_weighted_mean(dataset):
    column_map = _mdot_read_column_map(dataset["hst_file"])
    time_column = _mdot_find_column(column_map, ("time",))
    rate_column = _mdot_find_column(
        column_map, ("Mdot_flux", "mdot_flux")
    )
    if time_column is None or rate_column is None:
        raise KeyError(
            f"time or Mdot_flux column was not found: {dataset['hst_file']}"
        )

    data = _mdot_np.atleast_2d(
        _mdot_np.loadtxt(dataset["hst_file"], comments="#")
    )
    if data.shape[1] <= max(time_column, rate_column):
        raise ValueError(
            f"HST has only {data.shape[1]} columns: {dataset['hst_file']}"
        )

    time_raw = _mdot_np.asarray(data[:, time_column], dtype=float)
    rate_raw = _mdot_np.asarray(data[:, rate_column], dtype=float)
    finite = _mdot_np.isfinite(time_raw) & _mdot_np.isfinite(rate_raw)
    time_code, rate_code = _mdot_clean_restarts(
        time_raw[finite], rate_raw[finite]
    )
    if time_code.size < 2:
        raise RuntimeError(
            f"Fewer than two valid rows remain: {dataset['hst_file']}"
        )

    time_kyr = time_code * _MDOT_CODE_TIME_TO_KYR
    rate_msun_per_yr = rate_code * _MDOT_CODE_RATE_TO_MSUN_PER_YR
    elapsed_kyr = time_kyr[-1] - time_kyr[0]
    if elapsed_kyr <= 0.0:
        raise RuntimeError(f"Non-positive elapsed time: {dataset['hst_file']}")

    if hasattr(_mdot_np, "trapezoid"):
        trapezoid = _mdot_np.trapezoid
    else:
        trapezoid = _mdot_np.trapz
    time_weighted_mean = float(
        trapezoid(rate_msun_per_yr, x=time_kyr) / elapsed_kyr
    )

    return {
        **dataset,
        "time_min_kyr": float(time_kyr[0]),
        "time_max_kyr": float(time_kyr[-1]),
        "mean_msun_per_yr": time_weighted_mean,
        "raw_rows": int(time_raw.size),
        "cleaned_rows": int(time_code.size),
    }


_mdot_datasets = _mdot_collect_datasets(globals())
_mdot_results = [
    _mdot_load_time_weighted_mean(dataset)
    for dataset in _mdot_datasets
]

print("\n[INFO] Time-weighted mean accretion rates:")
for result in _mdot_results:
    print(
        f"  [{result['number']}] {result['curve_name']}\n"
        f"      rows = {result['raw_rows']} -> {result['cleaned_rows']}\n"
        f"      time = {result['time_min_kyr']:.6e} -- "
        f"{result['time_max_kyr']:.6e} kyr\n"
        f"      <Mdot>_t = {result['mean_msun_per_yr']:.6e} M_sun yr^-1"
    )


# 瞬間値は描画せず、各モデルの時間範囲に平均線だけを引く。
_mdot_fig, _mdot_ax = _mdot_plt.subplots(figsize=(8.5, 6.0))

if len(_mdot_results) <= 10:
    _mdot_colors = _mdot_plt.get_cmap("tab10").colors[:len(_mdot_results)]
else:
    _mdot_colors = _mdot_plt.get_cmap("turbo")(
        _mdot_np.linspace(0.05, 0.95, len(_mdot_results))
    )

for result, color in zip(_mdot_results, _mdot_colors):
    _mdot_ax.hlines(
        result["mean_msun_per_yr"],
        result["time_min_kyr"],
        result["time_max_kyr"],
        color=color,
        linewidth=2.2,
        linestyle="--",
        label=result["curve_name"],
    )

# Toyouchi et al. (2023) の臨界降着率 [M_sun yr^-1]
_mdot_ax.axhline(
    0.04,
    color="red",
    linewidth=2.0,
    linestyle="-",
    label=r"Toyouchi+23 $\dot{M}_{\rm crit}$",
)

_mdot_ax.set_xlabel("Time [kyr]")
_mdot_ax.set_ylabel(
    r"$\langle\dot{M}_{\rm flux}\rangle_t\ "
    r"[M_\odot\,{\rm yr}^{-1}]$"
)
_mdot_ax.set_title("Time-weighted mean accretion rates")
_mdot_ax.grid(True, linestyle="--", alpha=0.4)
_mdot_ax.margins(x=0.03, y=0.12)

# 情報ボックスを左上に仮配置し、実際の高さを測定する。
# その高さに応じてy軸上側へ空白を自動追加するため、水平線と重ならない。
_mdot_info_box = _mdot_ax.legend(
    loc="upper left",
    ncol=max(1, int(_mdot_np.ceil((len(_mdot_results) + 1) / 6))),
    frameon=True,
    fancybox=True,
    framealpha=0.92,
    facecolor="white",
    edgecolor="0.35",
    title="Models",
)
_mdot_info_box.get_title().set_ha("center")
_mdot_info_box.set_zorder(20)

_mdot_fig.tight_layout()
_mdot_fig.canvas.draw()
_mdot_renderer = _mdot_fig.canvas.get_renderer()
_mdot_legend_bbox = _mdot_info_box.get_window_extent(renderer=_mdot_renderer)
_mdot_axes_bbox = _mdot_ax.get_window_extent(renderer=_mdot_renderer)
_mdot_reserved_fraction = (
    _mdot_legend_bbox.height / _mdot_axes_bbox.height + 0.035
)

if _mdot_reserved_fraction >= 0.85:
    raise RuntimeError(
        "The information box is too tall for the graph. "
        "Shorten curve names or reduce the number of curves."
    )

_mdot_current_ymin, _mdot_current_ymax = _mdot_ax.get_ylim()
_mdot_data_ymax = max(
    [result["mean_msun_per_yr"] for result in _mdot_results]
    + [0.04]
)
_mdot_new_ymax = _mdot_current_ymin + (
    (_mdot_data_ymax - _mdot_current_ymin)
    / (1.0 - _mdot_reserved_fraction)
)
_mdot_ax.set_ylim(
    _mdot_current_ymin,
    max(_mdot_current_ymax, _mdot_new_ymax),
)

_mdot_fig.savefig(
    _MDOT_OUTPUT_FILE,
    dpi=200,
    bbox_inches="tight",
    facecolor="white",
)
_mdot_plt.show()
_mdot_plt.close(_mdot_fig)

print(f"\n[SAVED] {_MDOT_OUTPUT_FILE}")

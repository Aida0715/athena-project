# -*- coding: utf-8 -*-

# %% cell 1: 複数モデル比較（このセル単独で実行可）
"""仮想球表面の正味降着率 Mdot_net = Mdot_in - Mdot_out を比較する。"""
from pathlib import Path
import json
import os
import re

import matplotlib.pyplot as plt
import numpy as np


def load_comparison_metadata(data_path):
    directory = Path(data_path).expanduser().resolve()
    if directory.is_file():
        directory = directory.parent
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(
        directory.glob("*.metadata.json")
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one metadata JSON in {directory}; "
            f"found {len(candidates)}"
        )
    return json.loads(candidates[0].read_text(encoding="utf-8"))


def comparison_run_units(metadata):
    keys = (
        "mass_cgs", "length_cgs", "time_cgs", "velocity_cgs",
        "density_cgs", "pressure_cgs", "magnetic_field_gauss",
        "acceleration_cgs", "gravity_constant_code",
    )
    return {key: float(metadata["units"][key]) for key in keys}

# 2本目以降は datadir2 / curve_name2 を追加。
datadir1 = "~/athena-project/results/○○"
curve_name1 = r"$B_{z0}=○○\ \mu{\rm G}$"

# True: 各モデルの時間重み付き平均降着率を横線で表示
DRAW_TIME_AVERAGES = True

# True: 各JSONの共通fMをタイトルに表示
SHOW_COMMON_FM = True

HST_FILENAME = "Toyouchi.hst"
OUTPUT_DIR = Path("./virtual_sphere_net_accretion_rate").resolve()
OUTPUT_FILE = OUTPUT_DIR / "Mdot_net_model_comparison.png"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MSUN_CGS, YEAR_CGS = 1.98847e33, 365.25*86400.

def collect(namespace):
    numbers = sorted(int(m.group(1)) for key in namespace if (m := re.fullmatch(r"datadir(\d+)",key)))
    result = []
    for number in numbers:
        label_key = f"curve_name{number}"
        if label_key not in namespace: raise NameError(f"{label_key} is required")
        directory = Path(os.path.expanduser(str(namespace[f"datadir{number}"]))).resolve()
        path = directory if directory.is_file() else directory/HST_FILENAME
        if not path.is_file(): raise FileNotFoundError(path)
        metadata = load_comparison_metadata(path.parent)
        result.append((path,str(namespace[label_key]),metadata))
    if not result: raise RuntimeError("datadir1 is not configured")
    return result

def columns(path):
    result = {}
    with path.open("r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if line.lstrip().startswith("#"):
                for number,name in re.findall(r"\[(\d+)\]\s*=\s*([^\s]+)",line): result[name]=int(number)-1
    return result

def curve(path, metadata):
    units = comparison_run_units(metadata)
    time_to_kyr = units["time_cgs"]/(1.e3*YEAR_CGS)
    rate_to_msun_per_yr = ((units["mass_cgs"]/MSUN_CGS) /
                           (units["time_cgs"]/YEAR_CGS))
    cmap = columns(path)
    if "time" not in cmap or "Mdot_net" not in cmap:
        raise KeyError(f"time/Mdot_net column missing: {path}")
    table = np.loadtxt(path,comments="#",ndmin=2)
    time,rate = table[:,cmap["time"]],table[:,cmap["Mdot_net"]]
    good=np.isfinite(time)&np.isfinite(rate); time,rate=time[good],rate[good]
    order=np.argsort(time,kind="stable"); time,rate=time[order],rate[order]
    _,reverse_index=np.unique(time[::-1],return_index=True); keep=np.sort(len(time)-1-reverse_index)
    time,rate=time[keep],rate[keep]
    if len(time)<1: raise RuntimeError(f"At least one valid time is required: {path}")
    return time*time_to_kyr,rate*rate_to_msun_per_yr

def time_weighted_mean(time, rate):
    finite = np.isfinite(time) & np.isfinite(rate)
    time, rate = time[finite], rate[finite]
    if len(time) == 0:
        return np.nan
    if len(time) == 1 or not time[-1] > time[0]:
        return float(rate[-1])
    dt = np.diff(time)
    return float(np.sum(0.5*(rate[:-1]+rate[1:])*dt)/(time[-1]-time[0]))

datasets = collect(globals())
title = "Net accretion rate through the virtual-sphere surface"
if SHOW_COMMON_FM:
    fm_values = np.array([
        float(metadata["initial_conditions"]["fM"])
        for _, _, metadata in datasets
    ])
    if not np.allclose(fm_values, fm_values[0], rtol=1.0e-12, atol=0.0):
        raise ValueError(
            "SHOW_COMMON_FM=True, but fM is not common to all models: "
            f"{fm_values.tolist()}"
        )
    title += rf", $f_M={fm_values[0]:g}$"

fig,ax=plt.subplots(figsize=(9,5.5))
for path,label,metadata in datasets:
    time,rate=curve(path,metadata)
    line, = ax.plot(time,rate,lw=1.4,label=label)
    average = time_weighted_mean(time,rate)
    if DRAW_TIME_AVERAGES and np.isfinite(average):
        ax.hlines(
            average,time[0],time[-1],colors=line.get_color(),
            linestyles="--",lw=1.7,label=f"{label} time average",
        )
    print(f"[INFO] {label}: time-weighted mean(Mdot_net)={average:.6e} M_sun yr^-1")
ax.axhline(0.,color="0.35",lw=.8)
# loc="best"により、凡例ボックはデータとの重なりが最小の位置へ自動配置される。
ax.set(xlabel="Time [kyr]",ylabel=r"$\dot{M}_{\rm net}\ [M_\odot\,{\rm yr}^{-1}]$",title=title)
ax.grid(alpha=.3); ax.legend(title="Models",loc="best")
fig.tight_layout(); fig.savefig(OUTPUT_FILE,dpi=200,bbox_inches="tight")
plt.show(); plt.close(fig)
print(f"[SAVED] {OUTPUT_FILE}")


# %% cell 2: 1モデル用（このセル単独で実行可）
from pathlib import Path
import json
import os
import re

import matplotlib.pyplot as plt
import numpy as np


# -------------------------- ユーザー設定 --------------------------
SINGLE_DATADIR = Path(
    os.path.expanduser("~/athena-project/results/○○")
).resolve()
SINGLE_MODEL_LABEL = "Model"
SINGLE_HST_FILENAME = "Toyouchi.hst"
SINGLE_OUTPUT_DIR = SINGLE_DATADIR / "virtual_sphere_net_accretion_rate"
SINGLE_OUTPUT_FILE = SINGLE_OUTPUT_DIR / "Mdot_net_single_model.png"

MDOT_CRIT_MSUN_PER_YR = 0.04

SINGLE_MSUN_CGS = 1.98847e33
SINGLE_YEAR_CGS = 365.25 * 86400.0


def load_single_metadata(directory):
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(
        directory.glob("*.metadata.json")
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one metadata JSON in {directory}; "
            f"found {len(candidates)}"
        )
    return json.loads(candidates[0].read_text(encoding="utf-8"))


def single_hst_columns(path):
    result = {}
    with path.open("r",encoding="utf-8",errors="ignore") as handle:
        for line in handle:
            if line.lstrip().startswith("#"):
                for number,name in re.findall(
                    r"\[(\d+)\]\s*=\s*([^\s]+)",line
                ):
                    result[name] = int(number)-1
    return result


def read_single_mdot_net(directory):
    metadata = load_single_metadata(directory)
    hst_path = directory / SINGLE_HST_FILENAME
    if not hst_path.is_file():
        raise FileNotFoundError(hst_path)
    columns = single_hst_columns(hst_path)
    if "time" not in columns or "Mdot_net" not in columns:
        raise KeyError(f"time/Mdot_net column missing: {hst_path}")
    table = np.loadtxt(hst_path,comments="#",ndmin=2)
    time = table[:,columns["time"]]
    rate = table[:,columns["Mdot_net"]]
    valid = np.isfinite(time) & np.isfinite(rate)
    time,rate = time[valid],rate[valid]
    order = np.argsort(time,kind="stable")
    time,rate = time[order],rate[order]
    _,reverse_index = np.unique(time[::-1],return_index=True)
    keep = np.sort(len(time)-1-reverse_index)
    time,rate = time[keep],rate[keep]
    if len(time) < 1:
        raise RuntimeError(f"No valid time/Mdot_net rows: {hst_path}")
    units = metadata["units"]
    time_kyr = time*float(units["time_cgs"])/(1.0e3*SINGLE_YEAR_CGS)
    rate_msun_yr = rate*float(units["mass_cgs"])/SINGLE_MSUN_CGS/(
        float(units["time_cgs"])/SINGLE_YEAR_CGS
    )
    return time_kyr,rate_msun_yr,metadata


def single_time_weighted_mean(time,rate):
    if len(time) == 1 or not time[-1] > time[0]:
        return float(rate[-1])
    dt = np.diff(time)
    return float(np.sum(0.5*(rate[:-1]+rate[1:])*dt)/(time[-1]-time[0]))


single_time,single_rate,single_metadata = read_single_mdot_net(SINGLE_DATADIR)
single_average = single_time_weighted_mean(single_time,single_rate)
single_fm = float(single_metadata["initial_conditions"]["fM"])
single_bz0_microgauss = float(
    single_metadata["magnetic_field"]["initial_bz_microgauss"]
)

SINGLE_OUTPUT_DIR.mkdir(parents=True,exist_ok=True)
fig,ax = plt.subplots(figsize=(9,5.5))
instantaneous_line, = ax.plot(
    single_time,single_rate,lw=1.4,label=SINGLE_MODEL_LABEL
)
ax.hlines(
    single_average,single_time[0],single_time[-1],
    colors=instantaneous_line.get_color(),linestyles="--",lw=1.8,
    label=rf"Time average: ${single_average:.3g}\ M_\odot\,{{\rm yr}}^{{-1}}$",
)
ax.axhline(
    MDOT_CRIT_MSUN_PER_YR,color="red",linestyle="-",lw=1.8,
    label=rf"Toyouchi+23 $\dot{{M}}_{{\rm crit}}={MDOT_CRIT_MSUN_PER_YR:g}\ M_\odot\,{{\rm yr}}^{{-1}}$",
)
ax.axhline(0.0,color="0.35",lw=0.8)
ax.set(
    xlabel="Time [kyr]",
    ylabel=r"$\dot{M}_{\rm net}\ [M_\odot\,{\rm yr}^{-1}]$",
    title=(
        r"Net accretion rate through the virtual-sphere surface, "
        rf"$f_M={single_fm:g}$, "
        rf"$B_{{z0}}={single_bz0_microgauss:g}\,\mu{{\rm G}}$"
    ),
)
ax.grid(alpha=0.3)
ax.legend(loc="best")
fig.tight_layout()
fig.savefig(SINGLE_OUTPUT_FILE,dpi=200,bbox_inches="tight")
plt.show()
plt.close(fig)
print(f"[INFO] Time-weighted mean(Mdot_net)={single_average:.6e} M_sun yr^-1")
print(f"[SAVED] {SINGLE_OUTPUT_FILE}")

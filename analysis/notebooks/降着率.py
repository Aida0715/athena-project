# -*- coding: utf-8 -*-
"""シンクなし計算の全域ガス質量変化率を比較する。

これは中心星降着率ではなく HST の -dM_gas/dt であり、
外部境界の流出入や数値フロアの寄与も含み得る。
"""
from pathlib import Path
import os, re
import matplotlib.pyplot as plt
import numpy as np
try:
    from toyouchi_metadata import load_toyouchi_metadata, run_units
except ModuleNotFoundError:
    import json
    def load_toyouchi_metadata(data_path):
        directory = Path(data_path).expanduser().resolve(); directory = directory if directory.is_dir() else directory.parent
        preferred = directory / "Toyouchi.metadata.json"; candidates = [preferred] if preferred.is_file() else sorted(directory.glob("*.metadata.json"))
        if len(candidates) != 1: raise FileNotFoundError(f"Expected one metadata JSON in {directory}; found {len(candidates)}")
        metadata = json.loads(candidates[0].read_text(encoding="utf-8")); metadata["_path"] = str(candidates[0]); return metadata
    def run_units(metadata):
        keys = ("mass_cgs", "length_cgs", "time_cgs", "velocity_cgs", "density_cgs", "pressure_cgs", "magnetic_field_gauss", "acceleration_cgs", "gravity_constant_code")
        return {key: float(metadata["units"][key]) for key in keys}

# 2本目以降は datadir2 / curve_name2 を追加。
datadir1 = "~/athena-project/results/○○"
curve_name1 = r"$B_{z0}=○○\ \mu{\rm G}$"

HST_FILENAME = "Toyouchi.hst"
OUTPUT_DIR = Path("./gas_mass_change_rate").resolve()
OUTPUT_FILE = OUTPUT_DIR / "gas_mass_change_rate.png"
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
        result.append((path,str(namespace[label_key])))
    if not result: raise RuntimeError("datadir1 is not configured")
    return result

def columns(path):
    result = {}
    with path.open("r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if line.lstrip().startswith("#"):
                for number,name in re.findall(r"\[(\d+)\]\s*=\s*([^\s]+)",line): result[name]=int(number)-1
    return result

def curve(path):
    units = run_units(load_toyouchi_metadata(path.parent))
    time_to_kyr = units["time_cgs"]/(1.e3*YEAR_CGS)
    rate_to_msun_per_yr = ((units["mass_cgs"]/MSUN_CGS) /
                           (units["time_cgs"]/YEAR_CGS))
    cmap = columns(path)
    if "time" not in cmap or "mass" not in cmap: raise KeyError(f"time/mass column missing: {path}")
    table = np.loadtxt(path,comments="#",ndmin=2)
    time,mass = table[:,cmap["time"]],table[:,cmap["mass"]]
    good=np.isfinite(time)&np.isfinite(mass); time,mass=time[good],mass[good]
    order=np.argsort(time,kind="stable"); time,mass=time[order],mass[order]
    _,reverse_index=np.unique(time[::-1],return_index=True); keep=np.sort(len(time)-1-reverse_index)
    time,mass=time[keep],mass[keep]
    if len(time)<2: raise RuntimeError(f"At least two unique times required: {path}")
    return time*time_to_kyr,-np.gradient(mass,time)*rate_to_msun_per_yr

fig,ax=plt.subplots(figsize=(9,5.5))
for path,label in collect(globals()):
    time,rate=curve(path); ax.plot(time,rate,lw=1.4,label=label)
    finite=np.isfinite(rate)
    if np.any(finite): print(f"[INFO] {label}: mean(-dMgas/dt)={np.mean(rate[finite]):.6e} M_sun yr^-1")
ax.axhline(0.,color="0.35",lw=.8)
ax.set(xlabel="Time [kyr]",ylabel=r"$-dM_{\rm gas}/dt\ [M_\odot\,{\rm yr}^{-1}]$",title="Global gas-mass change rate (not stellar accretion)")
ax.grid(alpha=.3); ax.legend(title="Models",loc="best")
fig.tight_layout(); fig.savefig(OUTPUT_FILE,dpi=200,bbox_inches="tight")
plt.show(); plt.close(fig)
print(f"[SAVED] {OUTPUT_FILE}")

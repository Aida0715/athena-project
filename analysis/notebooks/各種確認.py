# -*- coding: utf-8 -*-
"""Toyouchi: シンクなしの現行計算向け HST/VTK 健全性確認。"""

from pathlib import Path
from collections import defaultdict
import csv, gc, re
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from toyouchi_metadata import load_toyouchi_metadata, run_units

# ======================== ユーザー設定 ========================
DATA_DIR = Path("~/athena-project/results/○○").expanduser()
HST_FILENAME = "Toyouchi.hst"
VTK_GLOB = "**/*.block*.out2.*.vtk"
OUTPUT_DIR = DATA_DIR / "checks"
SHOW_PLOTS = True

# 現行の基準単位
RUN_UNITS = run_units(load_toyouchi_metadata(DATA_DIR))
M_UNIT_G = RUN_UNITS["mass_cgs"]
L_UNIT_CM = RUN_UNITS["length_cgs"]
T_UNIT_S = RUN_UNITS["time_cgs"]
MSUN_G = 1.98847e33
TIME_KYR = T_UNIT_S / (1.e3 * 365.25 * 86400.)
MASS_MSUN = M_UNIT_G / MSUN_G
V_UNIT = L_UNIT_CM / T_UNIT_S
RHO_UNIT = M_UNIT_G / L_UNIT_CM**3
B_UNIT = RUN_UNITS["magnetic_field_gauss"]
E_UNIT = M_UNIT_G * V_UNIT**2
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def hst_columns(path):
    result = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.lstrip().startswith("#"):
                for number, name in re.findall(r"\[(\d+)\]\s*=\s*([^\s]+)", line):
                    result[name] = int(number)-1
    return result

def load_hst(path):
    columns = hst_columns(path)
    if "time" not in columns:
        raise KeyError("HST header has no time column")
    data = np.loadtxt(path, comments="#", ndmin=2)
    data = data[np.isfinite(data[:, columns["time"]])]
    data = data[np.argsort(data[:, columns["time"]], kind="stable")]
    _, reverse_index = np.unique(data[::-1, columns["time"]], return_index=True)
    return columns, data[np.sort(len(data)-1-reverse_index)]

def plot_hst(columns, data):
    time = data[:, columns["time"]] * TIME_KYR
    series = []
    if "mass" in columns:
        series.append(("Gas mass", data[:, columns["mass"]]*MASS_MSUN, r"$M_{\rm gas}\ [M_\odot]$"))
    if "Mstar" in columns:  # 固定中心星は現行コードに存在する
        series.append(("Central-star mass", data[:, columns["Mstar"]]*MASS_MSUN, r"$M_\star\ [M_\odot]$"))
    if not series:
        return
    fig, axes = plt.subplots(len(series), 1, figsize=(8, 3.2*len(series)), sharex=True, squeeze=False)
    for ax, (title, values, ylabel) in zip(axes[:, 0], series):
        ax.plot(time, values, lw=1.5); ax.set(title=title, ylabel=ylabel); ax.grid(alpha=.3)
    axes[-1, 0].set_xlabel("Time [kyr]")
    fig.tight_layout(); fig.savefig(OUTPUT_DIR/"history_overview.png", dpi=180)
    if SHOW_PLOTS: plt.show()
    plt.close(fig)

def vtk_step(path):
    match = re.search(r"\.out\d+\.(\d+)\.vtk$", path.name)
    return int(match.group(1)) if match else None

def grids(dataset):
    if isinstance(dataset, pv.MultiBlock):
        for block in dataset:
            if block is not None: yield from grids(block)
    else:
        yield dataset

def analyze_step(step, paths):
    out = defaultdict(float)
    out.update(step=step, blocks=0, cells=0, nonfinite=0, nonpositive_rho=0,
               rho_min=np.inf, rho_max=-np.inf, b_max=0., time=np.nan)
    times = []
    for path in paths:
        dataset = pv.read(path)
        for grid in grids(dataset):
            sized = grid.compute_cell_sizes(length=False, area=False, volume=True)
            vol = np.asarray(sized.cell_data["Volume"], dtype=float)
            rho = np.asarray(grid.cell_data["rho"], dtype=float).reshape(-1)
            out["blocks"] += 1; out["cells"] += rho.size
            out["nonfinite"] += np.count_nonzero(~np.isfinite(rho))
            out["nonpositive_rho"] += np.count_nonzero(rho <= 0.)
            good = np.isfinite(rho) & np.isfinite(vol) & (vol > 0.)
            if np.any(good):
                out["mass"] += np.sum(rho[good]*vol[good]); out["volume"] += np.sum(vol[good])
                out["rho_min"] = min(out["rho_min"], np.min(rho[good])); out["rho_max"] = max(out["rho_max"], np.max(rho[good]))
            if "Bcc" in grid.cell_data:
                b = np.asarray(grid.cell_data["Bcc"], dtype=float).reshape(-1, 3)
                b2 = np.einsum("ij,ij->i", b, b); bgood = good & np.isfinite(b2)
                if np.any(bgood):
                    out["magnetic_energy"] += np.sum(.5*b2[bgood]*vol[bgood]); out["b_max"] = max(out["b_max"], np.sqrt(np.max(b2[bgood])))
            if "TIME" in grid.field_data: times.extend(np.ravel(grid.field_data["TIME"]).astype(float))
            del sized
        del dataset; gc.collect()
    if times: out["time"] = float(np.median(times))
    return dict(out)

def save_vtk(rows):
    if not rows:
        print("[WARNING] VTK files not found; HST check only."); return
    fields = ["step","time_kyr","blocks","cells","volume_cm3","mass_msun","magnetic_energy_erg","rho_min_gcm3","rho_max_gcm3","b_max_gauss","nonfinite","nonpositive_rho"]
    with (OUTPUT_DIR/"vtk_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for r in rows:
            writer.writerow(dict(step=r["step"], time_kyr=r["time"]*TIME_KYR, blocks=int(r["blocks"]), cells=int(r["cells"]), volume_cm3=r["volume"]*L_UNIT_CM**3, mass_msun=r["mass"]*MASS_MSUN, magnetic_energy_erg=r["magnetic_energy"]*E_UNIT, rho_min_gcm3=r["rho_min"]*RHO_UNIT, rho_max_gcm3=r["rho_max"]*RHO_UNIT, b_max_gauss=r["b_max"]*B_UNIT, nonfinite=int(r["nonfinite"]), nonpositive_rho=int(r["nonpositive_rho"])))
    t = np.array([r["time"] for r in rows])*TIME_KYR
    values = [(np.array([r["mass"] for r in rows])*MASS_MSUN, r"Gas mass [$M_\odot$]", False), (np.array([r["magnetic_energy"] for r in rows])*E_UNIT, "Magnetic energy [erg]", True), (np.array([r["rho_min"] for r in rows])*RHO_UNIT, r"Minimum density [g cm$^{-3}$]", True), (np.array([r["rho_max"] for r in rows])*RHO_UNIT, r"Maximum density [g cm$^{-3}$]", True)]
    fig, axes = plt.subplots(2,2,figsize=(10,7),sharex=True)
    for ax,(y,label,log) in zip(axes.flat, values):
        ax.plot(t,y,marker="o",ms=2); ax.set_ylabel(label); ax.grid(alpha=.3)
        if log and np.any(y>0): ax.set_yscale("log")
    axes[1,0].set_xlabel("Time [kyr]"); axes[1,1].set_xlabel("Time [kyr]")
    fig.tight_layout(); fig.savefig(OUTPUT_DIR/"vtk_overview.png",dpi=180)
    if SHOW_PLOTS: plt.show()
    plt.close(fig)

hst_path = DATA_DIR/HST_FILENAME
if hst_path.is_file():
    cmap, table = load_hst(hst_path); plot_hst(cmap, table); print(f"[OK] HST rows: {len(table)}")
else: print(f"[WARNING] HST file not found: {hst_path}")
groups = defaultdict(list)
for path in sorted(DATA_DIR.glob(VTK_GLOB)):
    step = vtk_step(path)
    if step is not None: groups[step].append(path)
rows = []
for index,(step,paths) in enumerate(sorted(groups.items()),1):
    print(f"[VTK {index}/{len(groups)}] step={step:05d}, files={len(paths)}")
    rows.append(analyze_step(step,paths))
save_vtk(rows)
print(f"[SAVED] {OUTPUT_DIR}")

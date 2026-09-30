# -*- coding: utf-8 -*-
"""複数計算について、計算領域全体の磁場エネルギー時間変化を描く。"""

import gc
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


# ============================================================
# 1. ユーザー設定
#
# 2本目以降は、次の2行を番号を増やして追記するだけでよい。
# datadir2 = "~/athena-project/results/○○"
# curve_name2 = r"$B_{z,0}=1\times10^{-1}\ \mu{\rm G}$"
# ============================================================
datadir1 = "~/athena-project/results/○○"
curve_name1 = r"$B_{z0}=3\times10^{-2}\ \mu{\rm G}$"


# 以下は通常変更不要
vtk_output_number = 2
output_dir = Path("./magnetic_energy_history").resolve()
output_file = output_dir / "total_magnetic_energy_vs_time.png"
output_dir.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. Toyouchi.cppの新単位系（L0 = rb / 2）
# ============================================================
M_UNIT_CGS = 4.0e33       # g
L_UNIT_CGS = 7.03e15      # cm
T_UNIT_CGS = 3.61e10      # s

YEAR_CGS = 365.25 * 24.0 * 3600.0
TIME_UNIT_KYR = T_UNIT_CGS / (1.0e3 * YEAR_CGS)

V_UNIT_CGS = L_UNIT_CGS / T_UNIT_CGS
RHO_UNIT_CGS = M_UNIT_CGS / L_UNIT_CGS**3
P_UNIT_CGS = RHO_UNIT_CGS * V_UNIT_CGS**2
B_UNIT_GAUSS = np.sqrt(4.0 * np.pi * P_UNIT_CGS)

# code magnetic energy (B_code^2/2 * V_code) -> erg
MAGNETIC_ENERGY_UNIT_ERG = P_UNIT_CGS * L_UNIT_CGS**3

print(f"[INFO] 1 code time = {TIME_UNIT_KYR:.6e} kyr")
print(f"[INFO] 1 code magnetic field = {B_UNIT_GAUSS:.6e} G")
print(
    f"[INFO] 1 code energy = {MAGNETIC_ENERGY_UNIT_ERG:.6e} erg"
)


# ============================================================
# 3. datadir1, datadir2, ... と対応するcurve_nameを自動取得
# ============================================================
def collect_datasets(namespace):
    numbers = sorted(
        int(match.group(1))
        for name in namespace
        if (match := re.fullmatch(r"datadir(\d+)", name))
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

        directory = Path(
            os.path.expanduser(str(namespace[datadir_key]))
        ).resolve()
        if not directory.is_dir():
            raise NotADirectoryError(
                f"[{datadir_key}] Directory does not exist: {directory}"
            )

        datasets.append({
            "number": number,
            "datadir": directory,
            "curve_name": str(namespace[curve_name_key]),
        })
    return datasets


# ============================================================
# 4. Athena++ legacy VTKの時刻を取得
# ============================================================
def read_athena_vtk_time(filename):
    with open(filename, "rb") as handle:
        handle.readline()
        title = handle.readline().decode("ascii", errors="ignore")

    match = re.search(
        r"time\s*=\s*"
        r"([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        title,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTK time was not found: {filename}")
    return float(match.group(1))


def find_cell_array_name(grid, candidates):
    for name in candidates:
        if name in grid.cell_data:
            return name
    raise KeyError(
        "Magnetic-field cell data were not found.\n"
        f"Candidates: {candidates}\n"
        f"Available cell_data: {list(grid.cell_data.keys())}"
    )


# ============================================================
# 5. 1 MeshBlockの磁場エネルギー
#
# Athena++の磁場規約では
#   u_B,code = (Bx^2 + By^2 + Bz^2) / 2
# CGSでは
#   u_B = B_phys^2 / (8*pi)
# ============================================================
def block_magnetic_energy_erg(filename):
    grid = pv.read(filename)
    try:
        magnetic_name = find_cell_array_name(
            grid,
            ("Bcc", "bcc", "B", "magnetic_field"),
        )
        magnetic_field = np.asarray(
            grid.cell_data[magnetic_name], dtype=np.float64
        ).reshape(-1, 3)

        if magnetic_field.shape[0] != grid.n_cells:
            raise ValueError(
                f"Bcc cell count does not match grid.n_cells: {filename}"
            )

        finite = np.all(np.isfinite(magnetic_field), axis=1)
        if not np.all(finite):
            raise ValueError(f"Bcc contains NaN or inf: {filename}")

        # Toyouchi計算はCartesianで、各MeshBlock内ではセル体積が一様。
        xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
        block_volume_code = (
            (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        )
        if grid.n_cells <= 0 or block_volume_code <= 0.0:
            raise ValueError(f"Invalid MeshBlock volume: {filename}")
        cell_volume_code = block_volume_code / float(grid.n_cells)

        b_squared_code = np.einsum(
            "ij,ij->i", magnetic_field, magnetic_field
        )
        energy_code = 0.5 * cell_volume_code * np.sum(
            b_squared_code, dtype=np.float64
        )
        return float(energy_code * MAGNETIC_ENERGY_UNIT_ERG)
    finally:
        del grid


# ============================================================
# 6. 1計算分を逐次解析
# ============================================================
def analyze_dataset(dataset):
    directory = dataset["datadir"]
    vtk_pattern = (
        f"*.block*.out{vtk_output_number}.*.vtk"
    )
    vtk_files = sorted(directory.glob(vtk_pattern))
    if not vtk_files:
        raise FileNotFoundError(
            "MeshBlock VTK files were not found.\n"
            f"Checked: {directory / vtk_pattern}"
        )

    step_pattern = re.compile(
        rf"\.out{vtk_output_number}\.(\d+)\.vtk$"
    )
    files_by_step = defaultdict(list)
    for filename in vtk_files:
        match = step_pattern.search(filename.name)
        if match:
            files_by_step[int(match.group(1))].append(filename)
    if not files_by_step:
        raise RuntimeError(f"No VTK output steps were found: {directory}")

    steps = sorted(files_by_step)
    time_kyr_values = []
    energy_erg_values = []
    block_count_values = []

    print(
        f"\n[DATASET {dataset['number']}] {dataset['curve_name']}\n"
        f"  directory = {directory}\n"
        f"  steps     = {len(steps)}"
    )

    for index, step in enumerate(steps, start=1):
        step_files = files_by_step[step]
        time_code = read_athena_vtk_time(step_files[0])
        total_energy_erg = 0.0

        # 1ブロックずつ読み、加算後に破棄するため全VTKを保持しない。
        for block_index, filename in enumerate(step_files, start=1):
            total_energy_erg += block_magnetic_energy_erg(filename)
            if block_index % 128 == 0:
                gc.collect()

        time_kyr_values.append(time_code * TIME_UNIT_KYR)
        energy_erg_values.append(total_energy_erg)
        block_count_values.append(len(step_files))

        print(
            f"  [{index:3d}/{len(steps):3d}] "
            f"step={step:05d}, "
            f"t={time_kyr_values[-1]:.3f} kyr, "
            f"blocks={len(step_files)}, "
            f"E_B={total_energy_erg:.6e} erg"
        )
        gc.collect()

    time_kyr = np.asarray(time_kyr_values, dtype=float)
    energy_erg = np.asarray(energy_erg_values, dtype=float)
    block_counts = np.asarray(block_count_values, dtype=int)
    steps_array = np.asarray(steps, dtype=int)

    # restart後などでVTK番号と物理時刻の順が異なる場合にも時刻順で描く。
    order = np.argsort(time_kyr, kind="stable")
    time_kyr = time_kyr[order]
    energy_erg = energy_erg[order]
    block_counts = block_counts[order]
    steps_array = steps_array[order]

    numeric_file = output_dir / (
        f"magnetic_energy_{dataset['number']:02d}.txt"
    )
    np.savetxt(
        numeric_file,
        np.column_stack((steps_array, time_kyr, energy_erg, block_counts)),
        fmt=["%d", "%.10e", "%.10e", "%d"],
        header="output_step time_kyr total_magnetic_energy_erg meshblock_count",
    )
    print(f"  [SAVED] {numeric_file}")

    return {
        **dataset,
        "time_kyr": time_kyr,
        "energy_erg": energy_erg,
        "block_counts": block_counts,
    }


datasets = collect_datasets(globals())
results = [analyze_dataset(dataset) for dataset in datasets]


# ============================================================
# 7. 描画
# ============================================================
fig, ax = plt.subplots(figsize=(8.5, 6.0))

if len(results) <= 10:
    colors = plt.get_cmap("tab10").colors[:len(results)]
else:
    colors = plt.get_cmap("turbo")(
        np.linspace(0.05, 0.95, len(results))
    )

all_energies = []
for result, color in zip(results, colors):
    ax.plot(
        result["time_kyr"],
        result["energy_erg"],
        color=color,
        linewidth=1.9,
        label=result["curve_name"],
    )
    all_energies.append(result["energy_erg"])

all_energies = np.concatenate(all_energies)
finite_energies = all_energies[np.isfinite(all_energies)]
if finite_energies.size == 0:
    raise RuntimeError("No finite magnetic-energy values were found.")

# 全値が正なら、桁変化を比較しやすい対数軸を使用する。
use_log_scale = np.all(finite_energies > 0.0)
if use_log_scale:
    ax.set_yscale("log")
else:
    print(
        "[WARNING] Zero magnetic energy is present; "
        "the y axis is kept linear."
    )

ax.set_xlabel("Time [kyr]")
ax.set_ylabel(r"Total magnetic energy $E_B$ [erg]")
ax.set_title("Total magnetic energy in the computational domain")
ax.grid(True, which="both", linestyle="--", alpha=0.35)
ax.margins(x=0.03, y=0.08)

# 曲線との重なりが最小になる位置をMatplotlibが自動選択する。
information_box = ax.legend(
    loc="best",
    frameon=True,
    fancybox=True,
    framealpha=0.92,
    facecolor="white",
    edgecolor="0.35",
    title="Models",
)
information_box.set_zorder(20)

fig.tight_layout()
fig.savefig(
    output_file,
    dpi=200,
    bbox_inches="tight",
    facecolor="white",
)
plt.show()
plt.close(fig)

print(f"\n[SAVED] {output_file}")

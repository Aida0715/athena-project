# -*- coding: utf-8 -*-
"""複数計算の指定領域内における体積積分圧力比の時間発展を描く。

    beta_global = integral(P_gas dV) / integral(P_mag dV)
    P_gas,code = rho * cs^2
    P_mag,code = (Bx^2 + By^2 + Bz^2) / 2

Athena++の磁場規格化を使うと、圧力と体積の物理単位変換係数は
分子と分母で打ち消し合う。数値出力には確認用に両積分値 [erg]も保存する。
"""

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
# 最新run/Toyouchi-test71.shのOUTDIRと初期磁場を初期値にしている。
# 2本目以降は次の2行を、番号を増やして追記するだけでよい。
# datadir2 = "~/athena-project/results/〇〇"
# curve_name2 = r"$B_{z,0}=0.1\ \mu{\rm G}$"
# ============================================================
datadir1 = "~/athena-project/results/〇〇"
curve_name1 = r"$B_{z0}=20\ \mu{\rm G}$"


# 解析領域：原点中心の直方体
# analysis_fraction=0.1 なら、x/y/z各軸の幅を実データ全域幅の1/10にする。
analysis_fraction = 0.1
analysis_center_code = (0.0, 0.0, 0.0)

# 範囲をAUで直接指定する場合はanalysis_bounds_auを設定する。
# 例: analysis_bounds_au = ((-1e4, 1e4), (-1e4, 1e4), (-5e3, 5e3))
analysis_bounds_au = None

# コード単位で指定したい場合の予備設定。
# analysis_bounds_auとの同時指定はできない。
# 例: analysis_bounds_code = ((-5.0, 5.0), (-5.0, 5.0), (-2.0, 2.0))
analysis_bounds_code = None

# athinput.Toyouchi_61の新単位系（L0=rb/2）に対応する等温音速。
ISO_SOUND_SPEED_CODE = 0.707
vtk_output_number = 2

output_dir = Path("./plasma_beta_history").resolve()
output_file = output_dir / "volume_integrated_plasma_beta_vs_time.png"
output_dir.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. Toyouchi.cppの新単位系（L0 = rb / 2）
# ============================================================
M_UNIT_CGS = 4.0e33       # g
L_UNIT_CGS = 7.03e15      # cm
T_UNIT_CGS = 3.61e10      # s

YEAR_CGS = 365.25 * 24.0 * 3600.0
AU_CGS = 1.495978707e13
TIME_UNIT_KYR = T_UNIT_CGS / (1.0e3 * YEAR_CGS)
V_UNIT_CGS = L_UNIT_CGS / T_UNIT_CGS
RHO_UNIT_CGS = M_UNIT_CGS / L_UNIT_CGS**3
P_UNIT_CGS = RHO_UNIT_CGS * V_UNIT_CGS**2
PV_UNIT_ERG = P_UNIT_CGS * L_UNIT_CGS**3

print(f"[INFO] 1 code time = {TIME_UNIT_KYR:.6e} kyr")
print(f"[INFO] cs(code)    = {ISO_SOUND_SPEED_CODE:.6g}")


# ============================================================
# 3. datadir1, datadir2, ... とcurve_nameを自動取得
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
        curve_key = f"curve_name{number}"
        if curve_key not in namespace:
            raise NameError(f"{datadir_key} に対応する {curve_key} がありません。")

        directory = Path(
            os.path.expanduser(str(namespace[datadir_key]))
        ).resolve()
        if not directory.is_dir():
            raise NotADirectoryError(f"[{datadir_key}] Directory does not exist: {directory}")

        datasets.append({
            "number": number,
            "datadir": directory,
            "curve_name": str(namespace[curve_key]),
        })
    return datasets


# ============================================================
# 4. VTKヘッダとセル配列
# ============================================================
def read_athena_vtk_time(filename):
    with open(filename, "rb") as handle:
        handle.readline()
        title = handle.readline().decode("ascii", errors="ignore")

    match = re.search(
        r"time\s*=\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        title,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"VTK time was not found: {filename}")
    return float(match.group(1))


def find_cell_array_name(grid, candidates, quantity):
    lower_names = {name.lower(): name for name in grid.cell_data.keys()}
    for candidate in candidates:
        if candidate.lower() in lower_names:
            return lower_names[candidate.lower()]
    raise KeyError(
        f"{quantity} cell data were not found.\n"
        f"Candidates: {candidates}\n"
        f"Available: {list(grid.cell_data.keys())}"
    )


def group_vtk_files(directory):
    pattern = f"*.block*.out{vtk_output_number}.*.vtk"
    step_pattern = re.compile(rf"\.out{vtk_output_number}\.(\d+)\.vtk$")
    files_by_step = defaultdict(list)
    for filename in sorted(directory.glob(pattern)):
        match = step_pattern.search(filename.name)
        if match:
            files_by_step[int(match.group(1))].append(filename)
    if not files_by_step:
        raise FileNotFoundError(f"No MeshBlock VTK files were found: {directory / pattern}")
    return files_by_step


# ============================================================
# 5. 実VTKから全計算領域と解析領域を決める
# ============================================================
def detect_domain_bounds(step_files):
    domain = np.array([np.inf, -np.inf, np.inf, -np.inf, np.inf, -np.inf])
    for filename in step_files:
        grid = pv.read(filename)
        try:
            bounds = np.asarray(grid.bounds, dtype=float)
            domain[0] = min(domain[0], bounds[0])
            domain[1] = max(domain[1], bounds[1])
            domain[2] = min(domain[2], bounds[2])
            domain[3] = max(domain[3], bounds[3])
            domain[4] = min(domain[4], bounds[4])
            domain[5] = max(domain[5], bounds[5])
        finally:
            del grid
    if not np.all(np.isfinite(domain)):
        raise RuntimeError("The computational-domain bounds could not be determined.")
    return domain


def make_analysis_bounds(domain):
    if analysis_bounds_au is not None and analysis_bounds_code is not None:
        raise ValueError(
            "Specify only one of analysis_bounds_au and analysis_bounds_code."
        )

    if analysis_bounds_au is not None:
        bounds_array = np.asarray(analysis_bounds_au, dtype=float)
        if bounds_array.shape != (3, 2):
            raise ValueError(
                "analysis_bounds_au must have shape "
                "((xmin,xmax),(ymin,ymax),(zmin,zmax))."
            )
        bounds = (bounds_array * AU_CGS / L_UNIT_CGS).reshape(-1)
    elif analysis_bounds_code is not None:
        bounds_array = np.asarray(analysis_bounds_code, dtype=float)
        if bounds_array.shape != (3, 2):
            raise ValueError("analysis_bounds_code must have shape ((xmin,xmax),(ymin,ymax),(zmin,zmax)).")
        bounds = bounds_array.reshape(-1)
    else:
        fraction = float(analysis_fraction)
        if not (0.0 < fraction <= 1.0):
            raise ValueError("analysis_fraction must satisfy 0 < analysis_fraction <= 1.")
        center = np.asarray(analysis_center_code, dtype=float)
        if center.shape != (3,):
            raise ValueError("analysis_center_code must contain three values.")
        full_width = np.array([
            domain[1] - domain[0],
            domain[3] - domain[2],
            domain[5] - domain[4],
        ])
        half_width = 0.5 * fraction * full_width
        bounds = np.array([
            center[0] - half_width[0], center[0] + half_width[0],
            center[1] - half_width[1], center[1] + half_width[1],
            center[2] - half_width[2], center[2] + half_width[2],
        ])

    # 計算領域外へはみ出した部分は切り詰める。
    bounds[[0, 2, 4]] = np.maximum(bounds[[0, 2, 4]], domain[[0, 2, 4]])
    bounds[[1, 3, 5]] = np.minimum(bounds[[1, 3, 5]], domain[[1, 3, 5]])
    if np.any(bounds[[1, 3, 5]] <= bounds[[0, 2, 4]]):
        raise ValueError(f"Analysis region does not overlap the domain: {bounds}")
    return bounds


# ============================================================
# 6. 1 MeshBlockの圧力体積を積分
# ============================================================
def integrate_block(filename, region):
    grid = pv.read(filename)
    centers = None
    try:
        block_bounds = np.asarray(grid.bounds, dtype=float)
        if (
            block_bounds[1] <= region[0] or block_bounds[0] >= region[1]
            or block_bounds[3] <= region[2] or block_bounds[2] >= region[3]
            or block_bounds[5] <= region[4] or block_bounds[4] >= region[5]
        ):
            return 0.0, 0.0, 0.0, 0

        rho_name = find_cell_array_name(
            grid, ("rho", "dens", "density", "prim_dens"), "density"
        )
        magnetic_name = find_cell_array_name(
            grid, ("Bcc", "bcc", "B", "magnetic_field"), "magnetic field"
        )
        rho = np.asarray(grid.cell_data[rho_name], dtype=np.float64).reshape(-1)
        magnetic = np.asarray(
            grid.cell_data[magnetic_name], dtype=np.float64
        ).reshape(-1, 3)
        if rho.size != grid.n_cells or magnetic.shape[0] != grid.n_cells:
            raise ValueError(f"Cell-data size mismatch: {filename}")

        cell_counts = np.asarray(grid.dimensions, dtype=int) - 1
        if np.any(cell_counts <= 0) or int(np.prod(cell_counts)) != grid.n_cells:
            raise ValueError(
                f"Unexpected VTK dimensions={grid.dimensions}, n_cells={grid.n_cells}: {filename}"
            )
        cell_width = np.array([
            (block_bounds[1] - block_bounds[0]) / cell_counts[0],
            (block_bounds[3] - block_bounds[2]) / cell_counts[1],
            (block_bounds[5] - block_bounds[4]) / cell_counts[2],
        ])

        centers = np.asarray(grid.cell_centers().points, dtype=np.float64)
        overlap_x = np.maximum(
            0.0,
            np.minimum(centers[:, 0] + 0.5 * cell_width[0], region[1])
            - np.maximum(centers[:, 0] - 0.5 * cell_width[0], region[0]),
        )
        overlap_y = np.maximum(
            0.0,
            np.minimum(centers[:, 1] + 0.5 * cell_width[1], region[3])
            - np.maximum(centers[:, 1] - 0.5 * cell_width[1], region[2]),
        )
        overlap_z = np.maximum(
            0.0,
            np.minimum(centers[:, 2] + 0.5 * cell_width[2], region[5])
            - np.maximum(centers[:, 2] - 0.5 * cell_width[2], region[4]),
        )
        overlap_volume = overlap_x * overlap_y * overlap_z
        selected = overlap_volume > 0.0
        finite = (
            selected
            & np.isfinite(rho)
            & np.all(np.isfinite(magnetic), axis=1)
            & (rho >= 0.0)
        )
        if not np.any(finite):
            return 0.0, 0.0, 0.0, 0

        volume = overlap_volume[finite]
        gas_pressure = rho[finite] * ISO_SOUND_SPEED_CODE**2
        magnetic_pressure = 0.5 * np.einsum(
            "ij,ij->i", magnetic[finite], magnetic[finite]
        )
        gas_pv_code = np.sum(gas_pressure * volume, dtype=np.float64)
        magnetic_pv_code = np.sum(magnetic_pressure * volume, dtype=np.float64)
        return (
            float(gas_pv_code),
            float(magnetic_pv_code),
            float(np.sum(volume, dtype=np.float64)),
            int(np.count_nonzero(finite)),
        )
    finally:
        del centers, grid


# ============================================================
# 7. 1計算分を低メモリで逐次解析
# ============================================================
def analyze_dataset(dataset):
    files_by_step = group_vtk_files(dataset["datadir"])
    steps = sorted(files_by_step)
    domain = detect_domain_bounds(files_by_step[steps[0]])
    region = make_analysis_bounds(domain)

    print(
        f"\n[DATASET {dataset['number']}] {dataset['curve_name']}\n"
        f"  directory = {dataset['datadir']}\n"
        f"  steps     = {len(steps)}\n"
        f"  domain    = x[{domain[0]:.6g},{domain[1]:.6g}] "
        f"y[{domain[2]:.6g},{domain[3]:.6g}] "
        f"z[{domain[4]:.6g},{domain[5]:.6g}] code\n"
        f"  analysis  = x[{region[0]:.6g},{region[1]:.6g}] "
        f"y[{region[2]:.6g},{region[3]:.6g}] "
        f"z[{region[4]:.6g},{region[5]:.6g}] code"
    )

    rows = []
    for index, step in enumerate(steps, start=1):
        step_files = files_by_step[step]
        time_code = read_athena_vtk_time(step_files[0])
        gas_pv_code = magnetic_pv_code = integrated_volume_code = 0.0
        selected_cells = 0

        for block_index, filename in enumerate(step_files, start=1):
            gas_part, magnetic_part, volume_part, cells_part = integrate_block(
                filename, region
            )
            gas_pv_code += gas_part
            magnetic_pv_code += magnetic_part
            integrated_volume_code += volume_part
            selected_cells += cells_part
            if block_index % 128 == 0:
                gc.collect()

        if integrated_volume_code <= 0.0:
            raise RuntimeError(f"No cell volume was integrated at step={step:05d}.")
        beta_global = (
            gas_pv_code / magnetic_pv_code
            if magnetic_pv_code > 0.0
            else np.inf
        )
        rows.append((
            step,
            time_code * TIME_UNIT_KYR,
            beta_global,
            gas_pv_code * PV_UNIT_ERG,
            magnetic_pv_code * PV_UNIT_ERG,
            integrated_volume_code,
            selected_cells,
            len(step_files),
        ))
        print(
            f"  [{index:3d}/{len(steps):3d}] step={step:05d}, "
            f"t={rows[-1][1]:.3f} kyr, beta={beta_global:.6e}, "
            f"cells={selected_cells}, blocks={len(step_files)}"
        )
        gc.collect()

    values = np.asarray(rows, dtype=float)
    order = np.argsort(values[:, 1], kind="stable")
    values = values[order]

    numeric_file = output_dir / f"volume_integrated_plasma_beta_{dataset['number']:02d}.txt"
    np.savetxt(
        numeric_file,
        values,
        fmt=["%d", "%.10e", "%.10e", "%.10e", "%.10e", "%.10e", "%d", "%d"],
        header=(
            "output_step time_kyr beta_global "
            "integral_Pgas_dV_erg integral_Pmag_dV_erg "
            "integrated_volume_code selected_cell_count meshblock_count"
        ),
    )
    print(f"  [SAVED] {numeric_file}")
    return {
        **dataset,
        "time_kyr": values[:, 1],
        "beta_global": values[:, 2],
        "domain": domain,
        "region": region,
    }


# ============================================================
# 8. 複数計算を1枚に描画
# ============================================================
datasets = collect_datasets(globals())
results = [analyze_dataset(dataset) for dataset in datasets]

fig, ax = plt.subplots(figsize=(8.5, 6.0))
if len(results) <= 10:
    colors = plt.get_cmap("tab10").colors[:len(results)]
else:
    colors = plt.get_cmap("turbo")(np.linspace(0.05, 0.95, len(results)))

positive_values = []
for result, color in zip(results, colors):
    valid = (
        np.isfinite(result["time_kyr"])
        & np.isfinite(result["beta_global"])
        & (result["beta_global"] > 0.0)
    )
    if not np.any(valid):
        print(f"[WARNING] No finite positive beta values: {result['curve_name']}")
        continue
    beta_values = result["beta_global"][valid]
    positive_values.append(beta_values)
    ax.plot(
        result["time_kyr"][valid],
        beta_values,
        color=color,
        linewidth=2.0,
        label=result["curve_name"],
    )

ax.axhline(1.0, color="black", linestyle="--", linewidth=1.3, label=r"$\beta=1$")
ax.set_yscale("log")
ax.set_xlabel("Time [kyr]")
ax.set_ylabel(r"$\beta_{\rm global}=\int P_{\rm gas}\,dV\,/\,\int P_{\rm mag}\,dV$")
if analysis_bounds_au is None and analysis_bounds_code is None:
    region_text = f"{analysis_fraction:g} of full domain"
else:
    region_text = "manually specified box"
ax.set_title("Volume-integrated plasma beta evolution\n" + region_text)
ax.grid(True, which="both", alpha=0.3, linestyle="--")
ax.legend(loc="best", framealpha=0.9)
fig.tight_layout()
fig.savefig(output_file, dpi=200, bbox_inches="tight", facecolor="white")
plt.show()
plt.close(fig)

print(f"\n[SAVED] {output_file}")
print("[DONE] Volume-integrated plasma-beta analysis finished.")

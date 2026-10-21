# -*- coding: utf-8 -*-
"""Toyouchi VTKの3次元球殻平均密度プロファイルを全時刻分作成する。"""

import gc
import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
try:
    from toyouchi_metadata import load_toyouchi_metadata, run_units
except ModuleNotFoundError:
    import json
    from pathlib import Path
    def load_toyouchi_metadata(data_path):
        directory = Path(data_path).expanduser().resolve(); directory = directory if directory.is_dir() else directory.parent
        preferred = directory / "Toyouchi.metadata.json"; candidates = [preferred] if preferred.is_file() else sorted(directory.glob("*.metadata.json"))
        if len(candidates) != 1: raise FileNotFoundError(f"Expected one metadata JSON in {directory}; found {len(candidates)}")
        metadata = json.loads(candidates[0].read_text(encoding="utf-8")); metadata["_path"] = str(candidates[0]); return metadata
    def run_units(metadata):
        keys = ("mass_cgs", "length_cgs", "time_cgs", "velocity_cgs", "density_cgs", "pressure_cgs", "magnetic_field_gauss", "acceleration_cgs", "gravity_constant_code")
        return {key: float(metadata["units"][key]) for key in keys}


# ============================================================
# 1. ユーザー設定
# ============================================================

# Toyouchi.block*.out2.*.vtk があるディレクトリ
vtk_dir = os.path.expanduser(
    "~/athena-project/results/〇〇"
)

# 画像の保存先。絶対パスでも指定可能。
output_dir = os.path.expanduser(
    "./radial_density_profiles"
)

# 中心座標 [code length]
center_code = np.array([0.0, 0.0, 0.0], dtype=float)

# 対数球殻ビン数
n_bins = 100

figsize = (7.2, 5.2)
dpi = 200
overwrite_existing = False

os.makedirs(output_dir, exist_ok=True)


# ============================================================
# 2. Toyouchi.cppが出力した実行時メタデータ
# ============================================================

METADATA = load_toyouchi_metadata(vtk_dir)
RUN_UNITS = run_units(METADATA)
M_UNIT_CGS = RUN_UNITS["mass_cgs"]
L_UNIT_CGS = RUN_UNITS["length_cgs"]
T_UNIT_CGS = RUN_UNITS["time_cgs"]

AU_CGS = 1.495978707e13
YEAR_CGS = 365.25 * 24.0 * 3600.0

RHO_UNIT_CGS = M_UNIT_CGS / L_UNIT_CGS**3
LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
TIME_UNIT_KYR = T_UNIT_CGS / (1.0e3 * YEAR_CGS)

print(f"=== Units from {METADATA['_path']} ===")
print(f"1 code length  = {LENGTH_UNIT_AU:.6f} AU")
print(f"1 code time    = {TIME_UNIT_KYR:.6f} kyr")
print(f"1 code density = {RHO_UNIT_CGS:.6e} g cm^-3")
print(f"Output         = {os.path.abspath(output_dir)}")


# ============================================================
# 3. VTKファイルを出力時刻ごとに整理
# ============================================================

vtk_pattern = os.path.join(vtk_dir, "Toyouchi.block*.out2.*.vtk")
vtk_files = glob.glob(vtk_pattern)

if not vtk_files:
    raise FileNotFoundError(
        "VTK files were not found.\n"
        f"Checked pattern: {vtk_pattern}"
    )


def block_number(path):
    match = re.search(r"block(\d+)", os.path.basename(path))
    return int(match.group(1)) if match else -1


step_dict = defaultdict(list)
for path in vtk_files:
    match = re.search(r"\.out2\.(\d+)\.vtk$", path)
    if match:
        step_dict[int(match.group(1))].append(path)

for step in step_dict:
    step_dict[step].sort(key=block_number)

timesteps = sorted(step_dict)
if not timesteps:
    raise RuntimeError("No output timestep could be extracted from VTK names.")

print(f"Found {len(vtk_files)} VTK files in {len(timesteps)} timesteps.")


# ============================================================
# 4. ヘッダ時刻と密度配列の取得
# ============================================================

def read_vtk_time(path):
    with open(path, "rb") as stream:
        stream.readline()
        title = stream.readline().decode("ascii", errors="ignore")

    match = re.search(
        r"time\s*=\s*([+\-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+\-]?\d+)?)",
        title,
        flags=re.IGNORECASE,
    )
    return float(match.group(1)) if match else np.nan


def density_name(grid):
    for name in ("rho", "dens", "density", "prim_dens", "prim_density"):
        if name in grid.cell_data:
            return name
    raise KeyError(
        "Density array was not found. "
        f"Available cell arrays: {list(grid.cell_data.keys())}"
    )


def cell_radius_code(grid):
    centers = np.asarray(grid.cell_centers().points, dtype=float)
    displacement = centers - center_code[None, :]
    return np.sqrt(np.sum(displacement * displacement, axis=1))


def uniform_cell_volume_code(grid):
    """Cartesian Athena++ MeshBlock内の1セルの体積を返す。"""
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    block_volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
    if grid.n_cells <= 0 or not np.isfinite(block_volume) or block_volume <= 0.0:
        raise ValueError("Invalid MeshBlock volume or cell count.")
    return block_volume / float(grid.n_cells)


# ============================================================
# 5. 1時刻の球殻平均密度
# ============================================================

def radial_density_profile(files):
    # 1回目：実在セル中心の最小・最大半径を取得
    radius_min_code = np.inf
    radius_max_code = -np.inf

    for path in files:
        grid = pv.read(path)
        radius = cell_radius_code(grid)
        valid = np.isfinite(radius) & (radius > 0.0)
        if np.any(valid):
            radius_min_code = min(radius_min_code, float(np.min(radius[valid])))
            radius_max_code = max(radius_max_code, float(np.max(radius[valid])))
        del radius, grid

    if not (
        np.isfinite(radius_min_code)
        and np.isfinite(radius_max_code)
        and radius_min_code > 0.0
        and radius_max_code > radius_min_code
    ):
        raise RuntimeError("Could not determine the radial range from VTK cells.")

    bin_edges_code = np.geomspace(
        radius_min_code,
        np.nextafter(radius_max_code, np.inf),
        n_bins + 1,
    )

    density_volume_sum = np.zeros(n_bins, dtype=float)
    volume_sum = np.zeros(n_bins, dtype=float)
    cell_count = np.zeros(n_bins, dtype=np.int64)

    # 2回目：AMRのセル数で偏らないよう、体積重み付きで球殻積分
    for path in files:
        grid = pv.read(path)
        radius = cell_radius_code(grid)
        rho_code = np.asarray(grid.cell_data[density_name(grid)]).reshape(-1)
        cell_volume = uniform_cell_volume_code(grid)

        valid = (
            np.isfinite(radius)
            & (radius >= radius_min_code)
            & (radius <= radius_max_code)
            & np.isfinite(rho_code)
            & (rho_code > 0.0)
        )

        if np.any(valid):
            indices = np.searchsorted(
                bin_edges_code,
                radius[valid],
                side="right",
            ) - 1
            indices = np.clip(indices, 0, n_bins - 1)

            density_volume_sum += np.bincount(
                indices,
                weights=rho_code[valid] * cell_volume,
                minlength=n_bins,
            )
            volume_sum += np.bincount(
                indices,
                weights=np.full(indices.size, cell_volume),
                minlength=n_bins,
            )
            cell_count += np.bincount(indices, minlength=n_bins)

        del radius, rho_code, grid

    rho_shell_code = np.divide(
        density_volume_sum,
        volume_sum,
        out=np.full(n_bins, np.nan, dtype=float),
        where=volume_sum > 0.0,
    )

    radius_centers_au = (
        np.sqrt(bin_edges_code[:-1] * bin_edges_code[1:])
        * LENGTH_UNIT_AU
    )

    return {
        "radius_au": radius_centers_au,
        "rho_cgs": rho_shell_code * RHO_UNIT_CGS,
        "cell_count": cell_count,
        "radius_min_au": radius_min_code * LENGTH_UNIT_AU,
        "radius_max_au": radius_max_code * LENGTH_UNIT_AU,
    }


# ============================================================
# 6. 全タイムステップを描画・保存
# ============================================================

saved_files = []

for number, step in enumerate(timesteps, start=1):
    files = step_dict[step]
    output_file = os.path.join(output_dir, f"radial_density_{step:05d}.png")

    if os.path.isfile(output_file) and not overwrite_existing:
        print(
            f"[SKIP {number:3d}/{len(timesteps):3d}] "
            f"step={step:05d}: {output_file}"
        )
        saved_files.append(output_file)
        continue

    print(
        f"[PROCESS {number:3d}/{len(timesteps):3d}] "
        f"step={step:05d}"
    )

    time_code = read_vtk_time(files[0])
    time_kyr = time_code * TIME_UNIT_KYR
    profile = radial_density_profile(files)

    valid = (
        np.isfinite(profile["radius_au"])
        & (profile["radius_au"] > 0.0)
        & np.isfinite(profile["rho_cgs"])
        & (profile["rho_cgs"] > 0.0)
        & (profile["cell_count"] > 0)
    )

    if not np.any(valid):
        raise RuntimeError(f"No valid radial bins at step={step:05d}.")

    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(
        profile["radius_au"][valid],
        profile["rho_cgs"][valid],
        color="black",
        linewidth=1.6,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(profile["radius_min_au"], profile["radius_max_au"])
    ax.set_xlabel(r"Radius $r$ [AU]")
    ax.set_ylabel(r"Density $\rho$ [g cm$^{-3}$]")
    ax.set_title(
        "Spherically averaged density profile\n"
        f"step={step:05d}, t={time_kyr:.1f} kyr"
    )

    ax.grid(
        True,
        which="both",
        linestyle=":",
        linewidth=0.6,
        alpha=0.35,
    )
    fig.tight_layout()
    fig.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    saved_files.append(output_file)
    print(
        f"[SAVED] step={step:05d}, t={time_kyr:.1f} kyr\n"
        f"        radial range="
        f"{profile['radius_min_au']:.3e} -- "
        f"{profile['radius_max_au']:.3e} AU\n"
        f"        {output_file}"
    )

    del profile, valid
    gc.collect()


print("\n=== Finished ===")
print(f"Saved profiles : {len(saved_files)}")
print(f"Output directory: {os.path.abspath(output_dir)}")

# -*- coding: utf-8 -*-
"""SMR VTKから各半径の累積回転数を求める。

各時刻の円環で質量加重平均した角速度
    Omega(R,t) = sum(m_i v_phi,i/R_i) / sum(m_i)
を求め、時間方向に
    N(R) = integral |Omega(R,t)| dt / (2 pi)
を台形積分する。これは「その半径のガスが持つ代表的な
角速度を時間積分したEuler的回転数」であり、個々のガス粒子を
追跡した回転数ではない。
"""

from collections import defaultdict
from pathlib import Path
import csv
import gc
import json
import re

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


# ============================================================
# 1. ユーザー設定
# ============================================================
DATADIR = Path("~/athena-project/results/○○").expanduser().resolve()
VTK_GLOB = "Toyouchi.block*.out2.*.vtk"
STEPS = None                     # None: 全時刻、例: [0, 10, 20]

# 回転円盤の厚み。Noneならz方向の全計算領域を使う。
Z_HALF_WIDTH_AU = 2500.0

# 半径範囲。Noneなら計算領域から自動決定する。
R_MIN_AU = None
R_MAX_AU = None
N_RADIAL_BINS = 100
MIN_CELLS_PER_ANNULUS = 8

# True: 回転方向に関係なく回転した回数を積算。
# False: 逆回転を負とする符号付き累積回転数。
USE_ABSOLUTE_OMEGA = True

# Toyouchi.cppの初期回転
# v_phi=0.5*sqrt(G*Menc/R)を解析時間中保った場合の基準線。
DRAW_INITIAL_ROTATION_REFERENCE = True

OUTPUT_DIR = DATADIR / "checks"
OUTPUT_PNG = OUTPUT_DIR / "gas_total_rotations_vs_radius.png"
OUTPUT_CSV = OUTPUT_DIR / "gas_total_rotations_vs_radius.csv"

FIGSIZE = (8.5, 6.0)
DPI = 180
SHOW_IMAGE = False

AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0
KYR_CGS = 1.0e3 * YEAR_CGS


# ============================================================
# 2. metadataとVTK整理
# ============================================================
def load_metadata(directory):
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(
        directory.glob("*.metadata.json")
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one metadata JSON in {directory}; "
            f"found {len(candidates)}"
        )
    metadata = json.loads(candidates[0].read_text(encoding="utf-8"))
    if metadata.get("format") != "Toyouchi runtime metadata":
        raise ValueError(f"Unsupported metadata: {candidates[0]}")
    metadata["_path"] = str(candidates[0])
    return metadata


def collect_steps():
    groups = defaultdict(list)
    for path in DATADIR.glob(VTK_GLOB):
        match = re.search(r"\.(\d+)\.vtk$", path.name)
        if match:
            step = int(match.group(1))
            if STEPS is None or step in STEPS:
                groups[step].append(path)
    if not groups:
        raise FileNotFoundError(f"No VTK files: {DATADIR / VTK_GLOB}")
    return {step: sorted(paths) for step, paths in sorted(groups.items())}


def read_vtk_time(path):
    with path.open("rb") as handle:
        header = handle.read(4096).decode("ascii", errors="ignore")
    match = re.search(
        r"\btime\s*=\s*([-+0-9.eEdD]+)", header, flags=re.IGNORECASE
    )
    if not match:
        raise ValueError(f"VTK time was not found in {path}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def block_geometry(grid, path):
    bounds = np.asarray(grid.bounds, dtype=float).reshape(3, 2)
    dimensions = np.asarray(grid.dimensions, dtype=int) - 1
    if np.any(dimensions <= 0) or int(np.prod(dimensions)) != grid.n_cells:
        raise ValueError(
            f"Unexpected dimensions: {path}; dimensions={grid.dimensions}, "
            f"n_cells={grid.n_cells}"
        )
    spacing = (bounds[:, 1] - bounds[:, 0]) / dimensions
    if np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError(f"Invalid cell spacing: {path}; {spacing}")
    return bounds, spacing


# ============================================================
# 3. 共通半径ビンの決定
# ============================================================
def make_radial_edges(first_paths, full_bounds, length_au):
    finest_xy = np.inf
    for path in first_paths:
        grid = pv.read(path)
        try:
            _, spacing = block_geometry(grid, path)
            finest_xy = min(finest_xy, spacing[0], spacing[1])
        finally:
            del grid

    # 一部の方位だけでなく完全な円環が入る最大半径。
    auto_rmax = min(
        full_bounds[0, 1], -full_bounds[0, 0],
        full_bounds[1, 1], -full_bounds[1, 0],
    )
    rmin = (
        float(R_MIN_AU) / length_au
        if R_MIN_AU is not None else 2.0 * finest_xy
    )
    rmax = (
        float(R_MAX_AU) / length_au
        if R_MAX_AU is not None else auto_rmax
    )
    if not (np.isfinite(rmin) and np.isfinite(rmax) and 0.0 < rmin < rmax):
        raise ValueError(f"Invalid radial range: {rmin=} {rmax=} [code]")
    return np.geomspace(rmin, rmax, int(N_RADIAL_BINS) + 1)


# ============================================================
# 4. 1時刻の円環平均Omega
# ============================================================
def omega_profile(paths, radial_edges, z_half_width_code):
    nbins = radial_edges.size - 1
    ann_mass = np.zeros(nbins, dtype=np.float64)
    ann_momega = np.zeros(nbins, dtype=np.float64)
    ann_count = np.zeros(nbins, dtype=np.int64)

    for block_number, path in enumerate(paths, start=1):
        grid = pv.read(path)
        try:
            if "rho" not in grid.cell_data or "vel" not in grid.cell_data:
                raise KeyError(f"{path}: rho or vel is absent")
            rho = np.asarray(grid.cell_data["rho"], dtype=np.float64).reshape(-1)
            velocity = np.asarray(grid.cell_data["vel"], dtype=np.float64)
            xyz = np.asarray(grid.cell_centers().points, dtype=np.float64)
            _, spacing = block_geometry(grid, path)
            dx, dy, dz = spacing

            valid = (
                np.isfinite(rho) & (rho > 0.0)
                & np.all(np.isfinite(velocity), axis=1)
            )
            if z_half_width_code is None:
                dz_overlap = np.full(rho.shape, dz, dtype=np.float64)
            else:
                zlo = np.maximum(xyz[:, 2] - 0.5 * dz, -z_half_width_code)
                zhi = np.minimum(xyz[:, 2] + 0.5 * dz, z_half_width_code)
                dz_overlap = np.maximum(zhi - zlo, 0.0)
                valid &= dz_overlap > 0.0
            if not np.any(valid):
                continue

            x = xyz[valid, 0]
            y = xyz[valid, 1]
            radius = np.hypot(x, y)
            vel = velocity[valid]
            mass = rho[valid] * dx * dy * dz_overlap[valid]
            vphi = np.divide(
                -y * vel[:, 0] + x * vel[:, 1], radius,
                out=np.zeros_like(radius), where=radius > 0.0,
            )
            omega = np.divide(
                vphi, radius, out=np.zeros_like(radius), where=radius > 0.0
            )
            ibin = np.searchsorted(radial_edges, radius, side="right") - 1
            use = (
                (ibin >= 0) & (ibin < nbins)
                & np.isfinite(omega) & np.isfinite(mass) & (mass > 0.0)
            )
            np.add.at(ann_mass, ibin[use], mass[use])
            np.add.at(ann_momega, ibin[use], mass[use] * omega[use])
            np.add.at(ann_count, ibin[use], 1)
        finally:
            del grid
        if block_number % 128 == 0:
            gc.collect()

    omega = np.divide(
        ann_momega, ann_mass, out=np.full(nbins, np.nan),
        where=(ann_mass > 0.0) & (ann_count >= MIN_CELLS_PER_ANNULUS),
    )
    return omega, ann_count


# ============================================================
# 5. 全時刻の積分と出力
# ============================================================
def integrate_rotations(times, omega_profiles):
    order = np.argsort(times, kind="stable")
    times = np.asarray(times, dtype=float)[order]
    omega_profiles = np.asarray(omega_profiles, dtype=float)[order]
    if times.size < 2:
        raise RuntimeError("At least two VTK times are required")
    if np.any(np.diff(times) <= 0.0):
        raise RuntimeError("VTK times must be unique and strictly increasing")

    omega_used = np.abs(omega_profiles) if USE_ABSOLUTE_OMEGA else omega_profiles
    rotations = np.zeros(omega_profiles.shape[1], dtype=np.float64)
    covered_time = np.zeros_like(rotations)
    for index, dt in enumerate(np.diff(times)):
        pair_valid = (
            np.isfinite(omega_used[index])
            & np.isfinite(omega_used[index + 1])
        )
        rotations[pair_valid] += (
            0.5 * (omega_used[index, pair_valid] + omega_used[index + 1, pair_valid])
            * dt / (2.0 * np.pi)
        )
        covered_time[pair_valid] += dt
    total_time = times[-1] - times[0]
    coverage = covered_time / total_time
    rotations[coverage == 0.0] = np.nan
    return times, omega_profiles, rotations, coverage


def initial_menc_msun(radius_au, metadata):
    """Toyouchi.cppのMencMMHMsun()と同じ初期包有質量。"""
    initial = metadata["initial_conditions"]
    fit_rmin = float(initial["fit_rmin_au"])
    fit_rmax = float(initial["fit_rmax_au"])
    mass_coeff = np.asarray(initial["mass_fit_coefficients"], dtype=float)
    core_mass = float(initial["core_mass_msun"])
    core_a, core_b, core_c, core_d = np.asarray(
        initial["core_coefficients"], dtype=float
    )

    radius_au = np.asarray(radius_au, dtype=float)
    menc = np.full_like(radius_au, np.nan)
    inner = (radius_au > 0.0) & (radius_au < fit_rmin)
    if np.any(inner):
        t = radius_au[inner] / fit_rmin
        t2 = t * t
        menc[inner] = core_mass * t**3 * (
            core_a + t2 * (core_b + t2 * (core_c + core_d * t2))
        )
    outer = radius_au >= fit_rmin
    if np.any(outer):
        x = np.log10(np.minimum(radius_au[outer], fit_rmax))
        # C++のPolynomial()は係数を高次からHorner法で評価する。
        menc[outer] = 10.0 ** np.polyval(mass_coeff, x)
    return menc


def initial_rotation_reference(radius_au, duration_code, metadata):
    """Toyouchi.cppの初期v_phiをduration_code維持した回転数。"""
    units = metadata["units"]
    initial = metadata["initial_conditions"]
    length_au = float(units["length_cgs"]) / AU_CGS
    mass_code_per_msun = MSUN_CGS / float(units["mass_cgs"])
    gravity_code = float(units["gravity_constant_code"])
    fit_rmax_au = float(initial["fit_rmax_au"])

    radius_au = np.asarray(radius_au, dtype=float)
    radius_code = radius_au / length_au
    # Toyouchi.cppではfit範囲外でv_phiを端点値に固定する。
    rotation_radius_au = np.minimum(radius_au, fit_rmax_au)
    rotation_radius_code = rotation_radius_au / length_au
    menc_code = initial_menc_msun(rotation_radius_au, metadata) * mass_code_per_msun
    vphi_code = 0.5 * np.sqrt(
        gravity_code * menc_code / rotation_radius_code
    )
    omega_code = np.divide(
        vphi_code, radius_code,
        out=np.full_like(radius_code, np.nan), where=radius_code > 0.0,
    )
    return omega_code * duration_code / (2.0 * np.pi)


def save_outputs(
    radius_au, rotations, theory_rotations, coverage,
    times_code, time_cgs, z_label,
):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    time_kyr = times_code * time_cgs / KYR_CGS

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "radius_AU", "measured_total_rotations",
            "initial_rotation_reference", "measured_over_reference",
            "time_coverage_fraction",
            "start_time_kyr", "end_time_kyr",
        ])
        for radius, count, reference, fraction in zip(
            radius_au, rotations, theory_rotations, coverage
        ):
            ratio = (
                count / reference
                if np.isfinite(count) and np.isfinite(reference) and reference > 0.0
                else np.nan
            )
            writer.writerow([
                f"{radius:.10e}", f"{count:.10e}", f"{reference:.10e}",
                f"{ratio:.10e}", f"{fraction:.10e}",
                f"{time_kyr[0]:.10e}", f"{time_kyr[-1]:.10e}",
            ])

    fig, ax = plt.subplots(figsize=FIGSIZE)
    valid = np.isfinite(rotations) & (coverage > 0.0)
    ax.semilogx(
        radius_au[valid], rotations[valid], color="tab:blue", lw=2.0,
        label="Measured (mass-weighted)",
    )
    valid_theory = np.isfinite(theory_rotations) & (theory_rotations >= 0.0)
    if DRAW_INITIAL_ROTATION_REFERENCE and np.any(valid_theory):
        ax.semilogx(
            radius_au[valid_theory], theory_rotations[valid_theory],
            color="black", ls="--", lw=1.7,
            label=r"Initial $0.5\sqrt{GM_{\rm enc}/R}$ reference",
        )
    ax.set_xlabel("Radius [AU]")
    ax.set_ylabel("Total number of rotations")
    ax.set_title(
        "Cumulative gas rotations\n"
        rf"$t={time_kyr[0]:.1f}$--${time_kyr[-1]:.1f}$ kyr, {z_label}"
    )
    ax.grid(True, which="both", alpha=0.28)
    ax.axhline(0.0, color="black", lw=0.8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=DPI)
    if SHOW_IMAGE:
        plt.show()
    plt.close(fig)


def main():
    if not DATADIR.is_dir():
        raise FileNotFoundError(DATADIR)
    metadata = load_metadata(DATADIR)
    units = metadata["units"]
    length_au = float(units["length_cgs"]) / AU_CGS
    time_cgs = float(units["time_cgs"])
    mesh = metadata["mesh"]
    full_bounds = np.array([
        [mesh["x1min_code"], mesh["x1max_code"]],
        [mesh["x2min_code"], mesh["x2max_code"]],
        [mesh["x3min_code"], mesh["x3max_code"]],
    ], dtype=float)
    z_half_width_code = (
        None if Z_HALF_WIDTH_AU is None else float(Z_HALF_WIDTH_AU) / length_au
    )
    z_label = (
        "full z domain" if Z_HALF_WIDTH_AU is None
        else rf"$|z|\leq {Z_HALF_WIDTH_AU:g}$ AU"
    )

    groups = collect_steps()
    radial_edges = make_radial_edges(next(iter(groups.values())), full_bounds, length_au)
    radial_centers = np.sqrt(radial_edges[:-1] * radial_edges[1:])
    times = []
    profiles = []

    print(f"[INFO] Metadata : {metadata['_path']}")
    print(f"[INFO] Input    : {DATADIR}")
    print(f"[INFO] Output   : {OUTPUT_DIR}")
    print(f"[INFO] Steps    : {len(groups)}")
    print(f"[INFO] z range  : {z_label}")
    for number, (step, paths) in enumerate(groups.items(), start=1):
        time_code = read_vtk_time(paths[0])
        print(
            f"[RUN {number:3d}/{len(groups):3d}] step={step:05d}, "
            f"t={time_code * time_cgs / KYR_CGS:.1f} kyr, blocks={len(paths)}"
        )
        omega, _ = omega_profile(paths, radial_edges, z_half_width_code)
        times.append(time_code)
        profiles.append(omega)
        gc.collect()

    times, profiles, rotations, coverage = integrate_rotations(times, profiles)
    theory_rotations = initial_rotation_reference(
        radial_centers * length_au, times[-1] - times[0], metadata
    )
    save_outputs(
        radial_centers * length_au, rotations, theory_rotations, coverage,
        times, time_cgs, z_label,
    )
    print(f"[SAVED] {OUTPUT_PNG}")
    print(f"[SAVED] {OUTPUT_CSV}")


if __name__ == "__main__":
    main()

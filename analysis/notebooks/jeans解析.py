# -*- coding: utf-8 -*-
"""Toyouchi SMR VTKのJeans不安定性とJeans解像度を解析する。

球内部の不安定度:
    J(r) = Menc(r) / MJ[rho_bar(<r)]
    rho_bar(<r) = 3 Menc(r)/(4 pi r^3)
    MJ = pi^(5/2)/6 * cs^3/(G^(3/2) rho_bar^(1/2))

Jeans解像度:
    lambdaJ = cs * sqrt(pi/(G rho))
    NJ = lambdaJ / Delta_x

Jマップは球対称な1次元J(r)をx-y面へ射影したもの。
NJマップはz=0面を横切る実SMRセルの局所値を描画する。
"""

from collections import defaultdict
from pathlib import Path
import csv
import gc
import json
import re

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pyvista as pv


# ============================================================
# 1. ユーザー設定
# ============================================================
DATADIR = Path("~/athena-project/results/○○").expanduser().resolve()
VTK_GLOB = "Toyouchi.block*.out2.*.vtk"
STEPS = None                         # None: 全VTK時刻

OUTPUT_ROOT = DATADIR / "jeans_analysis"
OUTPUT_DIRS = {
    "instability": {
        "full": OUTPUT_ROOT / "spherical_instability_maps" / "full",
        "zoom_1over10": (
            OUTPUT_ROOT / "spherical_instability_maps" / "zoom_1over10"
        ),
    },
    "resolution": {
        "full": OUTPUT_ROOT / "jeans_resolution_maps" / "full",
        "zoom_1over10": (
            OUTPUT_ROOT / "jeans_resolution_maps" / "zoom_1over10"
        ),
    },
    "profiles": OUTPUT_ROOT / "radial_profiles",
    "tables": OUTPUT_ROOT / "tables",
}

# fullと各軸幅の1/10 zoom。同じVTK走査中に両方を集計する。
MAP_VIEWS = {"full": 1.0, "zoom_1over10": 0.1}
MAP_RESOLUTION = 700
MAP_AXIS_SCALE_AU = 1.0e4
N_RADIAL_BINS = 100
MIN_CELLS_PER_SHELL = 8

# 1を白とし、不安定側(J>1)を赤で表示。
INSTABILITY_DISPLAY_MIN = 0.1
INSTABILITY_DISPLAY_MAX = 10.0

# Jeans解像度は4セルを白とし、不足側を赤で表示。
JEANS_CENTER_CELLS = 4.0
JEANS_RESOLUTION_MIN = 1.0
JEANS_RESOLUTION_MAX = 32.0

FIGSIZE = (8.3, 7.0)
PROFILE_FIGSIZE = (8.3, 6.3)
DPI = 180
SHOW_IMAGES = False
OVERWRITE = False

AU_CGS = 1.495978707e13
KYR_CGS = 1.0e3 * 365.25 * 24.0 * 3600.0


# ============================================================
# 2. metadata・VTK整理
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
        raise ValueError(f"VTK time was not found: {path}")
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
        raise ValueError(f"Invalid spacing: {path}; {spacing}")
    return bounds, spacing


def map_bounds(full_bounds, fraction):
    width = full_bounds[:2, 1] - full_bounds[:2, 0]
    half = 0.5 * fraction * width
    return np.array([[-half[0], half[0]], [-half[1], half[1]]])


def inspect_finest_spacing(paths):
    finest = np.inf
    for path in paths:
        grid = pv.read(path)
        try:
            _, spacing = block_geometry(grid, path)
            finest = min(finest, *spacing)
        finally:
            del grid
    if not np.isfinite(finest):
        raise RuntimeError("Could not determine the finest cell width")
    return finest


def radial_edges(first_paths, full_bounds):
    finest = inspect_finest_spacing(first_paths)
    rmax = min(
        full_bounds[0, 1], -full_bounds[0, 0],
        full_bounds[1, 1], -full_bounds[1, 0],
        full_bounds[2, 1], -full_bounds[2, 0],
    )
    rmin = finest
    if not 0.0 < rmin < rmax:
        raise ValueError(f"Invalid radial interval: {rmin=} {rmax=}")
    # Mencの累積和から中心セルの質量を落とさないよう、
    # 最初の球殻だけは0 <= r < finestとする。
    return np.r_[0.0, np.geomspace(rmin, rmax, N_RADIAL_BINS)]


def add_rectangles(diff, ix0, ix1, iy0, iy1, values):
    np.add.at(diff, (iy0, ix0), values)
    np.add.at(diff, (iy1, ix0), -values)
    np.add.at(diff, (iy0, ix1), -values)
    np.add.at(diff, (iy1, ix1), values)


# ============================================================
# 3. 1時刻の集計
# ============================================================
def aggregate_step(paths, edges, bounds_by_view, cs_code, g_code):
    nbins = edges.size - 1
    shell_mass = np.zeros(nbins)
    shell_volume = np.zeros(nbins)
    shell_nj_mass = np.zeros(nbins)
    shell_nj_min = np.full(nbins, np.inf)
    shell_count = np.zeros(nbins, dtype=np.int64)

    nx = ny = int(MAP_RESOLUTION)
    map_buffers = {}
    for view_name, bounds_xy in bounds_by_view.items():
        x_edges = np.linspace(bounds_xy[0, 0], bounds_xy[0, 1], nx + 1)
        y_edges = np.linspace(bounds_xy[1, 0], bounds_xy[1, 1], ny + 1)
        map_buffers[view_name] = {
            "bounds": bounds_xy,
            "x_edges": x_edges,
            "y_edges": y_edges,
            "pixel_dx": x_edges[1] - x_edges[0],
            "pixel_dy": y_edges[1] - y_edges[0],
            "nj_diff": np.zeros((ny + 1, nx + 1)),
            "weight_diff": np.zeros((ny + 1, nx + 1)),
        }

    for block_number, path in enumerate(paths, start=1):
        grid = pv.read(path)
        try:
            if "rho" not in grid.cell_data:
                raise KeyError(f"{path}: rho is absent")
            rho = np.asarray(grid.cell_data["rho"], dtype=float).reshape(-1)
            xyz = np.asarray(grid.cell_centers().points, dtype=float)
            _, spacing = block_geometry(grid, path)
            dx, dy, dz = spacing
            cell_width = max(dx, dy, dz)  # 保守的に最大辺を使用
            volume = dx * dy * dz
            finite = np.isfinite(rho) & (rho > 0.0)
            if not np.any(finite):
                continue

            radius3d = np.linalg.norm(xyz, axis=1)
            lambda_j = np.full_like(rho, np.nan)
            lambda_j[finite] = cs_code * np.sqrt(np.pi / (g_code * rho[finite]))
            nj = lambda_j / cell_width
            ibin = np.searchsorted(edges, radius3d, side="right") - 1
            use = finite & np.isfinite(nj) & (ibin >= 0) & (ibin < nbins)
            if np.any(use):
                bins = ibin[use]
                mass = rho[use] * volume
                shell_mass += np.bincount(bins, weights=mass, minlength=nbins)
                shell_volume += np.bincount(
                    bins, weights=np.full(bins.size, volume), minlength=nbins
                )
                shell_nj_mass += np.bincount(
                    bins, weights=mass * nj[use], minlength=nbins
                )
                shell_count += np.bincount(bins, minlength=nbins)
                np.minimum.at(shell_nj_min, bins, nj[use])

            # z=0面を横切るセルだけをx-yマップへ射影。
            zlo = xyz[:, 2] - 0.5 * dz
            zhi = xyz[:, 2] + 0.5 * dz
            on_plane = finite & np.isfinite(nj) & (zlo <= 0.0) & (zhi > 0.0)
            if np.any(on_plane):
                x = xyz[on_plane, 0]
                y = xyz[on_plane, 1]
                xlo, xhi = x - 0.5*dx, x + 0.5*dx
                ylo, yhi = y - 0.5*dy, y + 0.5*dy
                plane_nj = nj[on_plane]
                for buffer in map_buffers.values():
                    bounds_xy = buffer["bounds"]
                    in_map = (
                        (xhi > bounds_xy[0, 0]) & (xlo < bounds_xy[0, 1])
                        & (yhi > bounds_xy[1, 0]) & (ylo < bounds_xy[1, 1])
                    )
                    if not np.any(in_map):
                        continue
                    ix0 = np.floor(
                        (xlo[in_map]-bounds_xy[0, 0])/buffer["pixel_dx"]
                    ).astype(int)
                    ix1 = np.ceil(
                        (xhi[in_map]-bounds_xy[0, 0])/buffer["pixel_dx"]
                    ).astype(int)
                    iy0 = np.floor(
                        (ylo[in_map]-bounds_xy[1, 0])/buffer["pixel_dy"]
                    ).astype(int)
                    iy1 = np.ceil(
                        (yhi[in_map]-bounds_xy[1, 0])/buffer["pixel_dy"]
                    ).astype(int)
                    ix0, ix1 = np.clip(ix0, 0, nx), np.clip(ix1, 0, nx)
                    iy0, iy1 = np.clip(iy0, 0, ny), np.clip(iy1, 0, ny)
                    good = (ix1 > ix0) & (iy1 > iy0)
                    values = plane_nj[in_map][good]
                    add_rectangles(
                        buffer["nj_diff"], ix0[good], ix1[good],
                        iy0[good], iy1[good], values,
                    )
                    add_rectangles(
                        buffer["weight_diff"], ix0[good], ix1[good],
                        iy0[good], iy1[good], np.ones(values.size),
                    )
        finally:
            del grid
        if block_number % 128 == 0:
            gc.collect()

    centers = np.sqrt(edges[:-1] * edges[1:])
    centers[0] = 0.5 * edges[1]
    enclosed_radius = edges[1:]
    menc = np.cumsum(shell_mass)
    # Mencは各球殻の外縁までの質量なので、平均密度も
    # 外縁半径の球体積で割る。
    rho_mean = np.divide(3.0*menc, 4.0*np.pi*enclosed_radius**3,
                         out=np.full_like(menc, np.nan), where=menc > 0.0)
    mj_mean = (np.pi**2.5/6.0) * cs_code**3 / (g_code**1.5 * np.sqrt(rho_mean))
    instability = np.divide(menc, mj_mean, out=np.full_like(menc, np.nan), where=mj_mean > 0.0)
    nj_mean = np.divide(shell_nj_mass, shell_mass, out=np.full_like(shell_mass, np.nan),
                        where=(shell_mass > 0.0) & (shell_count >= MIN_CELLS_PER_SHELL))
    shell_nj_min[~np.isfinite(shell_nj_min)] = np.nan
    shell_nj_min[shell_count < MIN_CELLS_PER_SHELL] = np.nan

    valid_j = np.isfinite(instability) & (instability > 0.0)
    maps = {}
    for view_name, buffer in map_buffers.items():
        nj_sum = np.cumsum(
            np.cumsum(buffer["nj_diff"][:-1, :-1], axis=0), axis=1
        )
        map_weight = np.cumsum(
            np.cumsum(buffer["weight_diff"][:-1, :-1], axis=0), axis=1
        )
        nj_map = np.divide(
            nj_sum, map_weight, out=np.full_like(nj_sum, np.nan),
            where=map_weight > 0,
        )
        xcent = 0.5 * (buffer["x_edges"][:-1] + buffer["x_edges"][1:])
        ycent = 0.5 * (buffer["y_edges"][:-1] + buffer["y_edges"][1:])
        xx, yy = np.meshgrid(xcent, ycent)
        rr = np.hypot(xx, yy)
        instability_map = np.full_like(rr, np.nan)
        if np.count_nonzero(valid_j) >= 2:
            inside = (
                (rr >= enclosed_radius[valid_j].min())
                & (rr <= enclosed_radius[valid_j].max())
            )
            instability_map[inside] = 10.0**np.interp(
                np.log10(rr[inside]), np.log10(enclosed_radius[valid_j]),
                np.log10(instability[valid_j]),
            )
        maps[view_name] = {
            "x_edges": buffer["x_edges"],
            "y_edges": buffer["y_edges"],
            "nj_map": nj_map,
            "instability_map": instability_map,
        }

    return {
        "radius": centers, "enclosed_radius": enclosed_radius,
        "menc": menc, "rho_mean": rho_mean, "mj_mean": mj_mean,
        "instability": instability, "maps": maps,
        "nj_mean": nj_mean, "nj_min": shell_nj_min,
        "shell_count": shell_count,
    }


# ============================================================
# 4. 描画・保存
# ============================================================
def extent_scaled(map_data, length_au):
    return [map_data["x_edges"][0]*length_au/MAP_AXIS_SCALE_AU,
            map_data["x_edges"][-1]*length_au/MAP_AXIS_SCALE_AU,
            map_data["y_edges"][0]*length_au/MAP_AXIS_SCALE_AU,
            map_data["y_edges"][-1]*length_au/MAP_AXIS_SCALE_AU]


def format_map(ax):
    exponent = int(np.log10(MAP_AXIS_SCALE_AU))
    ax.set_xlabel(rf"$x$ [$10^{{{exponent}}}$ AU]")
    ax.set_ylabel(rf"$y$ [$10^{{{exponent}}}$ AU]")


def centered_log_image(ax, values, extent, vmin, center, vmax, cmap):
    log_values = np.log10(values)
    return ax.imshow(
        log_values, origin="lower", extent=extent, cmap=cmap,
        norm=TwoSlopeNorm(vmin=np.log10(vmin), vcenter=np.log10(center),
                          vmax=np.log10(vmax)),
        interpolation="nearest", aspect="equal",
    )


def save_outputs(step, time_kyr, result, length_au, mass_msun):
    suffix = f"{step:05d}"
    for view_name, map_data in result["maps"].items():
        extent = extent_scaled(map_data, length_au)
        view_label = "full" if view_name == "full" else "1/10 zoom"

        fig, ax = plt.subplots(figsize=FIGSIZE)
        image = centered_log_image(
            ax, map_data["instability_map"], extent,
            INSTABILITY_DISPLAY_MIN, 1.0, INSTABILITY_DISPLAY_MAX, "RdBu_r"
        )
        format_map(ax)
        ax.set_title(
            rf"Spherical Jeans instability $M_{{\rm enc}}/M_J$, "
            rf"$t={time_kyr:.1f}$ kyr, {view_label}"
        )
        cbar = fig.colorbar(image, ax=ax, label=r"$M_{\rm enc}/M_J$")
        ticks = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
        cbar.set_ticks(np.log10(ticks))
        cbar.set_ticklabels([f"{x:g}" for x in ticks])
        fig.tight_layout()
        fig.savefig(
            OUTPUT_DIRS["instability"][view_name]
            / f"jeans_instability_{suffix}_{view_name}.png",
            dpi=DPI,
        )
        if SHOW_IMAGES: plt.show()
        plt.close(fig)

        fig, ax = plt.subplots(figsize=FIGSIZE)
        image = centered_log_image(
            ax, map_data["nj_map"], extent,
            JEANS_RESOLUTION_MIN, JEANS_CENTER_CELLS,
            JEANS_RESOLUTION_MAX, "RdBu"
        )
        format_map(ax)
        ax.set_title(
            rf"Jeans resolution $\lambda_J/\Delta x$ at $z=0$, "
            rf"$t={time_kyr:.1f}$ kyr, {view_label}"
        )
        cbar = fig.colorbar(image, ax=ax, label=r"$N_J=\lambda_J/\Delta x$")
        ticks = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
        cbar.set_ticks(np.log10(ticks))
        cbar.set_ticklabels([f"{x:g}" for x in ticks])
        fig.tight_layout()
        fig.savefig(
            OUTPUT_DIRS["resolution"][view_name]
            / f"jeans_resolution_{suffix}_{view_name}.png",
            dpi=DPI,
        )
        if SHOW_IMAGES: plt.show()
        plt.close(fig)

    radius_au = result["radius"] * length_au
    enclosed_radius_au = result["enclosed_radius"] * length_au
    fig, ax = plt.subplots(figsize=PROFILE_FIGSIZE)
    ax.loglog(radius_au, result["nj_mean"], lw=1.8, label="Mass-weighted mean")
    ax.loglog(radius_au, result["nj_min"], lw=1.2, label="Minimum in shell")
    for level, style in ((4.0, "-"), (8.0, "--"), (16.0, ":")):
        ax.axhline(level, color="black", ls=style, lw=1.0, label=rf"$N_J={level:g}$")
    ax.set_xlabel("Spherical radius r [AU]")
    ax.set_ylabel(r"Jeans resolution $N_J=\lambda_J/\Delta x$")
    ax.set_title(rf"Jeans-resolution profile, $t={time_kyr:.1f}$ kyr")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIRS["profiles"] / f"jeans_resolution_profile_{suffix}.png", dpi=DPI)
    if SHOW_IMAGES: plt.show()
    plt.close(fig)

    with (OUTPUT_DIRS["tables"] / f"jeans_profile_{suffix}.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "shell_center_radius_AU", "enclosed_outer_radius_AU",
            "Menc_Msun", "mean_density_code", "MJ_mean_Msun",
            "Menc_over_MJ", "NJ_mass_weighted_mean", "NJ_shell_min", "cell_count",
        ])
        for i in range(radius_au.size):
            writer.writerow([
                f"{radius_au[i]:.10e}", f"{enclosed_radius_au[i]:.10e}",
                f"{result['menc'][i]*mass_msun:.10e}",
                f"{result['rho_mean'][i]:.10e}", f"{result['mj_mean'][i]*mass_msun:.10e}",
                f"{result['instability'][i]:.10e}", f"{result['nj_mean'][i]:.10e}",
                f"{result['nj_min'][i]:.10e}", int(result["shell_count"][i]),
            ])


# ============================================================
# 5. 全時刻実行
# ============================================================
def main():
    if not DATADIR.is_dir():
        raise FileNotFoundError(DATADIR)
    if any(not (0.0 < fraction <= 1.0) for fraction in MAP_VIEWS.values()):
        raise ValueError("All MAP_VIEWS fractions must be in (0, 1]")
    metadata = load_metadata(DATADIR)
    units = metadata["units"]
    length_au = float(units["length_cgs"]) / AU_CGS
    mass_msun = float(units["mass_cgs"]) / 1.98847e33
    time_cgs = float(units["time_cgs"])
    cs_code = float(metadata["hydro"]["sound_speed_code"])
    g_code = float(units["gravity_constant_code"])
    mesh = metadata["mesh"]
    full_bounds = np.array([
        [mesh["x1min_code"], mesh["x1max_code"]],
        [mesh["x2min_code"], mesh["x2max_code"]],
        [mesh["x3min_code"], mesh["x3max_code"]],
    ], dtype=float)
    bounds_by_view = {
        name: map_bounds(full_bounds, fraction)
        for name, fraction in MAP_VIEWS.items()
    }
    groups = collect_steps()
    edges = radial_edges(next(iter(groups.values())), full_bounds)
    for value in OUTPUT_DIRS.values():
        if isinstance(value, dict):
            for directory in value.values():
                directory.mkdir(parents=True, exist_ok=True)
        else:
            value.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Metadata : {metadata['_path']}")
    print(f"[INFO] Input    : {DATADIR}")
    print(f"[INFO] Output   : {OUTPUT_ROOT}")
    print(f"[INFO] Steps    : {len(groups)}")
    print(f"[INFO] cs(code) : {cs_code:.8e}; G(code)={g_code:.8e}")

    for number, (step, paths) in enumerate(groups.items(), start=1):
        suffix = f"{step:05d}"
        expected_outputs = [
            OUTPUT_DIRS["profiles"] / f"jeans_resolution_profile_{suffix}.png"
        ]
        for view_name in MAP_VIEWS:
            expected_outputs.extend([
                OUTPUT_DIRS["instability"][view_name]
                / f"jeans_instability_{suffix}_{view_name}.png",
                OUTPUT_DIRS["resolution"][view_name]
                / f"jeans_resolution_{suffix}_{view_name}.png",
            ])
        if all(path.exists() for path in expected_outputs) and not OVERWRITE:
            print(f"[SKIP {number:3d}/{len(groups):3d}] step={step:05d}")
            continue
        time_code = read_vtk_time(paths[0])
        time_kyr = time_code * time_cgs / KYR_CGS
        print(f"[RUN  {number:3d}/{len(groups):3d}] step={step:05d}, t={time_kyr:.1f} kyr, blocks={len(paths)}")
        result = aggregate_step(paths, edges, bounds_by_view, cs_code, g_code)
        save_outputs(step, time_kyr, result, length_au, mass_msun)
        finite_j = result["instability"][np.isfinite(result["instability"])]
        zoom_nj = result["maps"]["zoom_1over10"]["nj_map"]
        finite_nj = zoom_nj[np.isfinite(zoom_nj)]
        if finite_j.size and finite_nj.size:
            print(f"       max(Menc/MJ)={finite_j.max():.3e}, min(NJ,z=0)={finite_nj.min():.3e}")
        del result
        gc.collect()
    print(f"[DONE] {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
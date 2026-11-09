# -*- coding: utf-8 -*-

# %% cell 1: x-z平面（このセル単独で実行可）
"""x-z平面のプラズマβマップ。磁力線の有無は手動で切り替える。"""

from collections import defaultdict
from pathlib import Path
import gc
import json
import os
import re

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, LogLocator
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata


# -------------------------- ユーザー設定 --------------------------
XZ_VTK_DIR = Path(os.path.expanduser("~/athena-project/results/○○")).resolve()

# True: プラズマβに細い黒色の投影磁力線を重ねる
XZ_DRAW_FIELDLINES = True
# True: metadataのracc内のβマップだけを白抜きにする
XZ_MASK_RACC = True

XZ_OUTPUT_ROOT = XZ_VTK_DIR / "xz_plasma_beta_maps"
XZ_OUTPUT_DIRS = {
    "full": XZ_OUTPUT_ROOT / "full",
    "zoom": XZ_OUTPUT_ROOT / "zoom_1over10",
}
XZ_BETA_VMIN = None
XZ_BETA_VMAX = None
XZ_BETA_CMAP = "plasma"
XZ_MAP_N = 700
XZ_ZOOM_FRACTION = 0.1
XZ_Y_SLICE_CODE = 0.0
XZ_STREAM_DENSITY = 1.1
XZ_STREAM_LINEWIDTH = 0.45
XZ_STREAM_ARROWSIZE = 0.55
XZ_FIGSIZE = (8.5, 8.0)
XZ_SAVE_DPI = 200
XZ_SAVE_IMAGES = True
XZ_SHOW_IMAGES = False
XZ_AU_CGS = 1.495978707e13
XZ_YEAR_CGS = 365.25 * 24.0 * 3600.0


def xz_load_metadata(directory):
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(directory.glob("*.metadata.json"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one metadata JSON in {directory}; found {len(candidates)}"
        )
    return json.loads(candidates[0].read_text(encoding="utf-8"))


XZ_METADATA = xz_load_metadata(XZ_VTK_DIR)
XZ_UNITS = XZ_METADATA["units"]
XZ_LENGTH_UNIT_AU = float(XZ_UNITS["length_cgs"]) / XZ_AU_CGS
XZ_TIME_UNIT_KYR = float(XZ_UNITS["time_cgs"]) / XZ_YEAR_CGS / 1.0e3
XZ_CS_CODE = float(XZ_METADATA["hydro"]["sound_speed_code"])
XZ_RACC_AU = float(XZ_METADATA["virtual_sphere"]["radius_au"]) if XZ_MASK_RACC else 0.0
XZ_BZ0 = float(XZ_METADATA["magnetic_field"]["initial_bz_microgauss"])
XZ_FM = float(XZ_METADATA["initial_conditions"]["fM"])
XZ_RUN_LABEL = rf"$B_{{z0}}={XZ_BZ0:g}\,\mu{{\rm G}}$, $f_M={XZ_FM:g}$"

for xz_directory in XZ_OUTPUT_DIRS.values():
    xz_directory.mkdir(parents=True, exist_ok=True)
if not XZ_VTK_DIR.is_dir():
    raise FileNotFoundError(XZ_VTK_DIR)


def xz_step_number(path):
    match = re.search(r"\.(\d+)\.vtk$", path.name)
    return int(match.group(1)) if match else None


def xz_read_time(path):
    with path.open("rb") as handle:
        header = handle.read(4096).decode("ascii", errors="ignore")
    match = re.search(r"\btime\s*=\s*([-+0-9.eEdD]+)", header, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"VTK time was not found: {path}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


XZ_FILES_BY_STEP = defaultdict(list)
for xz_path in sorted(XZ_VTK_DIR.glob("Toyouchi.block*.out2.*.vtk")):
    xz_step = xz_step_number(xz_path)
    if xz_step is not None:
        XZ_FILES_BY_STEP[xz_step].append(xz_path)
if not XZ_FILES_BY_STEP:
    raise FileNotFoundError(f"No out2 VTK files in {XZ_VTK_DIR}")
XZ_STEPS = sorted(XZ_FILES_BY_STEP)


def xz_find_array(grid, candidates):
    lower = {name.lower(): name for name in grid.cell_data.keys()}
    for candidate in candidates:
        if candidate in grid.cell_data:
            return candidate
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def xz_load_slice(step):
    x_parts, z_parts, rho_parts, b_parts = [], [], [], []
    xmin, xmax, zmin, zmax = np.inf, -np.inf, np.inf, -np.inf
    files = XZ_FILES_BY_STEP[step]
    for path in files:
        grid = pv.read(path)
        bounds = grid.bounds
        ymin, ymax = bounds[2], bounds[3]
        tolerance = 1.0e-12 * max(abs(ymin), abs(ymax), abs(XZ_Y_SLICE_CODE), 1.0)
        if not ymin-tolerance <= XZ_Y_SLICE_CODE <= ymax+tolerance:
            del grid
            continue
        rho_name = xz_find_array(grid, ("rho", "density"))
        b_name = xz_find_array(grid, ("Bcc", "bcc", "B", "magnetic_field"))
        if rho_name is None or b_name is None:
            raise KeyError(
                f"rho/Bcc is missing: {path.name}; cell_data={list(grid.cell_data.keys())}"
            )
        centers = np.asarray(grid.cell_centers().points)
        rho = np.asarray(grid.cell_data[rho_name]).reshape(-1)
        magnetic = np.asarray(grid.cell_data[b_name])[:, :3]
        distance = np.abs(centers[:, 1]-XZ_Y_SLICE_CODE)
        nearest = np.min(distance)
        layer = np.abs(distance-nearest) <= max(1.0e-12, nearest*1.0e-6)
        valid = layer & np.isfinite(rho) & (rho > 0.0) & np.all(np.isfinite(magnetic), axis=1)
        if np.any(valid):
            x_parts.append(centers[valid, 0])
            z_parts.append(centers[valid, 2])
            rho_parts.append(rho[valid])
            b_parts.append(magnetic[valid])
            xmin, xmax = min(xmin, bounds[0]), max(xmax, bounds[1])
            zmin, zmax = min(zmin, bounds[4]), max(zmax, bounds[5])
        del centers, rho, magnetic, grid
    if not x_parts:
        raise RuntimeError(f"No x-z cells at step={step:05d}")
    coordinates = np.column_stack((np.concatenate(x_parts), np.concatenate(z_parts)))
    rho = np.concatenate(rho_parts)
    magnetic = np.vstack(b_parts)
    unique, inverse = np.unique(coordinates, axis=0, return_inverse=True)
    count = np.bincount(inverse).astype(float)
    rho_average = np.bincount(inverse, weights=rho) / count
    b_average = np.column_stack([
        np.bincount(inverse, weights=magnetic[:, component]) / count
        for component in range(3)
    ])
    return {
        "step": step, "time_code": xz_read_time(files[0]),
        "x": unique[:, 0], "z": unique[:, 1], "rho": rho_average, "B": b_average,
        "bounds": np.array([xmin, xmax, zmin, zmax]),
    }


def xz_beta(rho, magnetic):
    b2 = np.einsum("ij,ij->i", magnetic, magnetic)
    result = np.full(rho.shape, np.nan)
    valid = np.isfinite(rho) & (rho > 0.0) & np.isfinite(b2) & (b2 > 0.0)
    result[valid] = 2.0*rho[valid]*XZ_CS_CODE**2/b2[valid]
    return result


def xz_beta_limits():
    if XZ_BETA_VMIN is not None and XZ_BETA_VMAX is not None:
        return float(XZ_BETA_VMIN), float(XZ_BETA_VMAX)
    samples = []
    for index, step in enumerate(XZ_STEPS, 1):
        data = xz_load_slice(step)
        values = xz_beta(data["rho"], data["B"])
        values = values[np.isfinite(values) & (values > 0.0)]
        if values.size:
            stride = max(1, values.size//20000)
            samples.append(values[::stride])
        print(f"[XZ SCALE {index:3d}/{len(XZ_STEPS):3d}] step={step:05d}")
        del data, values
        gc.collect()
    if not samples:
        raise RuntimeError("No finite positive plasma beta values")
    all_values = np.concatenate(samples)
    automatic_min, automatic_max = np.percentile(all_values, [1.0, 99.0])
    vmin = automatic_min if XZ_BETA_VMIN is None else float(XZ_BETA_VMIN)
    vmax = automatic_max if XZ_BETA_VMAX is None else float(XZ_BETA_VMAX)
    if not 0.0 < vmin < vmax:
        raise ValueError(f"Invalid beta limits: {vmin}, {vmax}")
    return vmin, vmax


def xz_axis_format(ax, bounds):
    largest = max(abs(value) for value in bounds)
    exponent = int(np.floor(np.log10(largest))) if largest > 0.0 else 0
    scale = 10.0**exponent
    formatter = FuncFormatter(lambda value, position: f"{value/scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\ {{\rm AU}}]$")
    ax.set_ylabel(rf"$z\ [10^{{{exponent}}}\ {{\rm AU}}]$")


def xz_interpolate(points, values, X, Z):
    try:
        result = griddata(points, values, (X, Z), method="linear")
    except Exception:
        result = np.full_like(X, np.nan, dtype=float)
    if np.any(~np.isfinite(result)):
        nearest = griddata(points, values, (X, Z), method="nearest")
        result = np.where(np.isfinite(result), result, nearest)
    return result


def xz_plot(data, fraction, output_dir, suffix, beta_vmin, beta_vmax):
    step = data["step"]
    bounds = data["bounds"]*XZ_LENGTH_UNIT_AU
    full_xmin, full_xmax, full_zmin, full_zmax = bounds
    if fraction is None:
        xmin, xmax, zmin, zmax = bounds
        range_name = "full domain"
    else:
        xhalf = 0.5*(full_xmax-full_xmin)*fraction
        zhalf = 0.5*(full_zmax-full_zmin)*fraction
        xmin, xmax, zmin, zmax = -xhalf, xhalf, -zhalf, zhalf
        range_name = "1/10 zoom"
    x_au, z_au = data["x"]*XZ_LENGTH_UNIT_AU, data["z"]*XZ_LENGTH_UNIT_AU
    beta = xz_beta(data["rho"], data["B"])
    margin_x, margin_z = 0.03*(xmax-xmin), 0.03*(zmax-zmin)
    source = (
        (x_au >= xmin-margin_x) & (x_au <= xmax+margin_x)
        & (z_au >= zmin-margin_z) & (z_au <= zmax+margin_z)
        & np.isfinite(beta) & (beta > 0.0)
    )
    if np.count_nonzero(source) < 4:
        raise RuntimeError(f"Too few x-z source cells: step={step:05d}")
    x_axis = np.linspace(xmin, xmax, XZ_MAP_N)
    z_axis = np.linspace(zmin, zmax, XZ_MAP_N)
    X, Z = np.meshgrid(x_axis, z_axis)
    points = np.column_stack((x_au[source], z_au[source]))
    log_beta = xz_interpolate(points, np.log10(beta[source]), X, Z)
    beta_grid = 10.0**log_beta
    if XZ_DRAW_FIELDLINES:
        bx_grid = np.nan_to_num(xz_interpolate(points, data["B"][source, 0], X, Z), nan=0.0)
        bz_grid = np.nan_to_num(xz_interpolate(points, data["B"][source, 2], X, Z), nan=0.0)
        has_field = np.any(np.hypot(bx_grid, bz_grid) > 0.0)
    else:
        bx_grid = bz_grid = None
        has_field = False
    if XZ_MASK_RACC:
        beta_grid = np.ma.masked_where(X*X+Z*Z < XZ_RACC_AU**2, beta_grid)
    fig, ax = plt.subplots(figsize=XZ_FIGSIZE)
    image = ax.pcolormesh(
        X, Z, beta_grid, shading="auto", cmap=XZ_BETA_CMAP,
        norm=LogNorm(vmin=beta_vmin, vmax=beta_vmax), rasterized=True,
    )
    if XZ_DRAW_FIELDLINES and has_field:
        ax.streamplot(
            x_axis, z_axis, bx_grid, bz_grid, color="black",
            density=XZ_STREAM_DENSITY, linewidth=XZ_STREAM_LINEWIDTH,
            arrowsize=XZ_STREAM_ARROWSIZE, integration_direction="both",
            broken_streamlines=False, minlength=0.03, maxlength=20.0, zorder=3,
        )
    cbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    cbar.set_label(r"Plasma $\beta=P_{\rm gas}/P_{\rm mag}$")
    cbar.locator = LogLocator(base=10.0)
    cbar.update_ticks()
    finite = beta_grid[np.isfinite(beta_grid)]
    if finite.size and finite.min() <= 1.0 <= finite.max():
        ax.contour(
            X, Z, beta_grid, levels=[1.0],
            colors="white" if XZ_DRAW_FIELDLINES else "black",
            linewidths=0.75, linestyles="--" if XZ_DRAW_FIELDLINES else "-", zorder=4,
        )
    ax.set(xlim=(xmin, xmax), ylim=(zmin, zmax))
    ax.set_aspect("equal", adjustable="box")
    xz_axis_format(ax, (xmin, xmax, zmin, zmax))
    title = (
        r"Plasma $\beta$ and projected magnetic field lines: x-z"
        if XZ_DRAW_FIELDLINES else r"Plasma $\beta$: x-z"
    )
    time_kyr = data["time_code"]*XZ_TIME_UNIT_KYR
    ax.set_title(
        f"{title}\nstep={step:05d}, t={time_kyr:.1f} kyr, "
        f"{range_name}, {XZ_RUN_LABEL}"
    )
    fig.tight_layout()
    tag = "with_Blines" if XZ_DRAW_FIELDLINES else "no_Blines"
    output_path = output_dir/f"xz_plasma_beta_{tag}_{step:05d}_{suffix}.png"
    if XZ_SAVE_IMAGES:
        fig.savefig(output_path, dpi=XZ_SAVE_DPI, bbox_inches="tight", facecolor="white")
        print(f"[SAVED] {output_path}")
    if XZ_SHOW_IMAGES:
        plt.show()
    plt.close(fig)
    del X, Z, beta_grid, log_beta, beta
    if bx_grid is not None:
        del bx_grid, bz_grid
    return output_path


XZ_BETA_LIMITS = xz_beta_limits()
print(f"[XZ] beta range = {XZ_BETA_LIMITS[0]:.6e} -- {XZ_BETA_LIMITS[1]:.6e}")
for xz_index, xz_step in enumerate(XZ_STEPS, 1):
    print(f"[XZ {xz_index:3d}/{len(XZ_STEPS):3d}] step={xz_step:05d}")
    xz_data = xz_load_slice(xz_step)
    xz_plot(xz_data, None, XZ_OUTPUT_DIRS["full"], "full", *XZ_BETA_LIMITS)
    xz_plot(xz_data, XZ_ZOOM_FRACTION, XZ_OUTPUT_DIRS["zoom"], "zoom_1over10", *XZ_BETA_LIMITS)
    del xz_data
    gc.collect()


# %% cell 2: x-y平面（このセル単独で実行可）
"""x-y平面のプラズマβマップ。full/zoomと磁力線の有無を手動指定する。"""

from collections import defaultdict
from pathlib import Path
import gc
import json
import os
import re

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, LogLocator
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata


# -------------------------- ユーザー設定 --------------------------
XY_VTK_DIR = Path(os.path.expanduser("~/athena-project/results/○○")).resolve()
# "full" または "zoom"。メモリ負荷を抑えるため一度に片方だけ実行する。
XY_PLOT_MODE = "zoom"
# True: プラズマβに細い黒色の投影磁力線を重ねる
XY_DRAW_FIELDLINES = True
# True: metadataのracc内のβマップだけを白抜きにする
XY_MASK_RACC = True

XY_OUTPUT_ROOT = XY_VTK_DIR / "xy_plasma_beta_maps"
XY_OUTPUT_DIRS = {"full": XY_OUTPUT_ROOT/"full", "zoom": XY_OUTPUT_ROOT/"zoom_1over10"}
XY_BETA_VMIN = None
XY_BETA_VMAX = None
XY_BETA_CMAP = "plasma"
XY_MAP_N = 700
XY_ZOOM_FRACTION = 0.1
XY_Z_SLICE_CODE = 0.0
XY_STREAM_DENSITY = 1.1
XY_STREAM_LINEWIDTH = 0.45
XY_STREAM_ARROWSIZE = 0.55
XY_FIGSIZE = (8.5, 8.0)
XY_SAVE_DPI = 200
XY_SAVE_IMAGES = True
XY_SHOW_IMAGES = False
XY_AU_CGS = 1.495978707e13
XY_YEAR_CGS = 365.25 * 24.0 * 3600.0

if XY_PLOT_MODE not in {"full", "zoom"}:
    raise ValueError('XY_PLOT_MODE must be "full" or "zoom".')


def xy_load_metadata(directory):
    preferred = directory / "Toyouchi.metadata.json"
    candidates = [preferred] if preferred.is_file() else sorted(directory.glob("*.metadata.json"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one metadata JSON in {directory}; found {len(candidates)}"
        )
    return json.loads(candidates[0].read_text(encoding="utf-8"))


XY_METADATA = xy_load_metadata(XY_VTK_DIR)
XY_UNITS = XY_METADATA["units"]
XY_LENGTH_UNIT_AU = float(XY_UNITS["length_cgs"]) / XY_AU_CGS
XY_TIME_UNIT_KYR = float(XY_UNITS["time_cgs"]) / XY_YEAR_CGS / 1.0e3
XY_CS_CODE = float(XY_METADATA["hydro"]["sound_speed_code"])
XY_RACC_AU = float(XY_METADATA["virtual_sphere"]["radius_au"]) if XY_MASK_RACC else 0.0
XY_BZ0 = float(XY_METADATA["magnetic_field"]["initial_bz_microgauss"])
XY_FM = float(XY_METADATA["initial_conditions"]["fM"])
XY_RUN_LABEL = rf"$B_{{z0}}={XY_BZ0:g}\,\mu{{\rm G}}$, $f_M={XY_FM:g}$"
for xy_directory in XY_OUTPUT_DIRS.values():
    xy_directory.mkdir(parents=True, exist_ok=True)
if not XY_VTK_DIR.is_dir():
    raise FileNotFoundError(XY_VTK_DIR)


def xy_step_number(path):
    match = re.search(r"\.(\d+)\.vtk$", path.name)
    return int(match.group(1)) if match else None


def xy_read_time(path):
    with path.open("rb") as handle:
        header = handle.read(4096).decode("ascii", errors="ignore")
    match = re.search(r"\btime\s*=\s*([-+0-9.eEdD]+)", header, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"VTK time was not found: {path}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


XY_FILES_BY_STEP = defaultdict(list)
for xy_path in sorted(XY_VTK_DIR.glob("Toyouchi.block*.out2.*.vtk")):
    xy_step = xy_step_number(xy_path)
    if xy_step is not None:
        XY_FILES_BY_STEP[xy_step].append(xy_path)
if not XY_FILES_BY_STEP:
    raise FileNotFoundError(f"No out2 VTK files in {XY_VTK_DIR}")
XY_STEPS = sorted(XY_FILES_BY_STEP)


def xy_find_array(grid, candidates):
    lower = {name.lower(): name for name in grid.cell_data.keys()}
    for candidate in candidates:
        if candidate in grid.cell_data:
            return candidate
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def xy_load_slice(step):
    x_parts, y_parts, rho_parts, b_parts = [], [], [], []
    xmin, xmax, ymin, ymax = np.inf, -np.inf, np.inf, -np.inf
    files = XY_FILES_BY_STEP[step]
    for path in files:
        grid = pv.read(path)
        bounds = grid.bounds
        zmin, zmax = bounds[4], bounds[5]
        tolerance = 1.0e-12 * max(abs(zmin), abs(zmax), abs(XY_Z_SLICE_CODE), 1.0)
        if not zmin-tolerance <= XY_Z_SLICE_CODE <= zmax+tolerance:
            del grid
            continue
        rho_name = xy_find_array(grid, ("rho", "density"))
        b_name = xy_find_array(grid, ("Bcc", "bcc", "B", "magnetic_field"))
        if rho_name is None or b_name is None:
            raise KeyError(
                f"rho/Bcc is missing: {path.name}; cell_data={list(grid.cell_data.keys())}"
            )
        centers = np.asarray(grid.cell_centers().points)
        rho = np.asarray(grid.cell_data[rho_name]).reshape(-1)
        magnetic = np.asarray(grid.cell_data[b_name])[:, :3]
        distance = np.abs(centers[:, 2]-XY_Z_SLICE_CODE)
        nearest = np.min(distance)
        layer = np.abs(distance-nearest) <= max(1.0e-12, nearest*1.0e-6)
        valid = layer & np.isfinite(rho) & (rho > 0.0) & np.all(np.isfinite(magnetic), axis=1)
        if np.any(valid):
            x_parts.append(centers[valid, 0])
            y_parts.append(centers[valid, 1])
            rho_parts.append(rho[valid])
            b_parts.append(magnetic[valid])
            xmin, xmax = min(xmin, bounds[0]), max(xmax, bounds[1])
            ymin, ymax = min(ymin, bounds[2]), max(ymax, bounds[3])
        del centers, rho, magnetic, grid
    if not x_parts:
        raise RuntimeError(f"No x-y cells at step={step:05d}")
    coordinates = np.column_stack((np.concatenate(x_parts), np.concatenate(y_parts)))
    rho = np.concatenate(rho_parts)
    magnetic = np.vstack(b_parts)
    unique, inverse = np.unique(coordinates, axis=0, return_inverse=True)
    count = np.bincount(inverse).astype(float)
    rho_average = np.bincount(inverse, weights=rho) / count
    b_average = np.column_stack([
        np.bincount(inverse, weights=magnetic[:, component]) / count
        for component in range(3)
    ])
    return {
        "step": step, "time_code": xy_read_time(files[0]),
        "x": unique[:, 0], "y": unique[:, 1], "rho": rho_average, "B": b_average,
        "bounds": np.array([xmin, xmax, ymin, ymax]),
    }


def xy_beta(rho, magnetic):
    b2 = np.einsum("ij,ij->i", magnetic, magnetic)
    result = np.full(rho.shape, np.nan)
    valid = np.isfinite(rho) & (rho > 0.0) & np.isfinite(b2) & (b2 > 0.0)
    result[valid] = 2.0*rho[valid]*XY_CS_CODE**2/b2[valid]
    return result


def xy_beta_limits():
    if XY_BETA_VMIN is not None and XY_BETA_VMAX is not None:
        return float(XY_BETA_VMIN), float(XY_BETA_VMAX)
    samples = []
    for index, step in enumerate(XY_STEPS, 1):
        data = xy_load_slice(step)
        values = xy_beta(data["rho"], data["B"])
        values = values[np.isfinite(values) & (values > 0.0)]
        if values.size:
            stride = max(1, values.size//20000)
            samples.append(values[::stride])
        print(f"[XY SCALE {index:3d}/{len(XY_STEPS):3d}] step={step:05d}")
        del data, values
        gc.collect()
    if not samples:
        raise RuntimeError("No finite positive plasma beta values")
    all_values = np.concatenate(samples)
    automatic_min, automatic_max = np.percentile(all_values, [1.0, 99.0])
    vmin = automatic_min if XY_BETA_VMIN is None else float(XY_BETA_VMIN)
    vmax = automatic_max if XY_BETA_VMAX is None else float(XY_BETA_VMAX)
    if not 0.0 < vmin < vmax:
        raise ValueError(f"Invalid beta limits: {vmin}, {vmax}")
    return vmin, vmax


def xy_axis_format(ax, bounds):
    largest = max(abs(value) for value in bounds)
    exponent = int(np.floor(np.log10(largest))) if largest > 0.0 else 0
    scale = 10.0**exponent
    formatter = FuncFormatter(lambda value, position: f"{value/scale:g}")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf"$x\ [10^{{{exponent}}}\ {{\rm AU}}]$")
    ax.set_ylabel(rf"$y\ [10^{{{exponent}}}\ {{\rm AU}}]$")


def xy_interpolate(points, values, X, Y):
    try:
        result = griddata(points, values, (X, Y), method="linear")
    except Exception:
        result = np.full_like(X, np.nan, dtype=float)
    if np.any(~np.isfinite(result)):
        nearest = griddata(points, values, (X, Y), method="nearest")
        result = np.where(np.isfinite(result), result, nearest)
    return result


def xy_plot(data, fraction, output_dir, suffix, beta_vmin, beta_vmax):
    step = data["step"]
    bounds = data["bounds"]*XY_LENGTH_UNIT_AU
    full_xmin, full_xmax, full_ymin, full_ymax = bounds
    if fraction is None:
        xmin, xmax, ymin, ymax = bounds
        range_name = "full domain"
    else:
        xhalf = 0.5*(full_xmax-full_xmin)*fraction
        yhalf = 0.5*(full_ymax-full_ymin)*fraction
        xmin, xmax, ymin, ymax = -xhalf, xhalf, -yhalf, yhalf
        range_name = "1/10 zoom"
    x_au, y_au = data["x"]*XY_LENGTH_UNIT_AU, data["y"]*XY_LENGTH_UNIT_AU
    beta = xy_beta(data["rho"], data["B"])
    margin_x, margin_y = 0.03*(xmax-xmin), 0.03*(ymax-ymin)
    source = (
        (x_au >= xmin-margin_x) & (x_au <= xmax+margin_x)
        & (y_au >= ymin-margin_y) & (y_au <= ymax+margin_y)
        & np.isfinite(beta) & (beta > 0.0)
    )
    if np.count_nonzero(source) < 4:
        raise RuntimeError(f"Too few x-y source cells: step={step:05d}")
    x_axis = np.linspace(xmin, xmax, XY_MAP_N)
    y_axis = np.linspace(ymin, ymax, XY_MAP_N)
    X, Y = np.meshgrid(x_axis, y_axis)
    points = np.column_stack((x_au[source], y_au[source]))
    log_beta = xy_interpolate(points, np.log10(beta[source]), X, Y)
    beta_grid = 10.0**log_beta
    if XY_DRAW_FIELDLINES:
        bx_grid = np.nan_to_num(xy_interpolate(points, data["B"][source, 0], X, Y), nan=0.0)
        by_grid = np.nan_to_num(xy_interpolate(points, data["B"][source, 1], X, Y), nan=0.0)
        has_field = np.any(np.hypot(bx_grid, by_grid) > 0.0)
    else:
        bx_grid = by_grid = None
        has_field = False
    if XY_MASK_RACC:
        beta_grid = np.ma.masked_where(X*X+Y*Y < XY_RACC_AU**2, beta_grid)
    fig, ax = plt.subplots(figsize=XY_FIGSIZE)
    image = ax.pcolormesh(
        X, Y, beta_grid, shading="auto", cmap=XY_BETA_CMAP,
        norm=LogNorm(vmin=beta_vmin, vmax=beta_vmax), rasterized=True,
    )
    if XY_DRAW_FIELDLINES and has_field:
        ax.streamplot(
            x_axis, y_axis, bx_grid, by_grid, color="black",
            density=XY_STREAM_DENSITY, linewidth=XY_STREAM_LINEWIDTH,
            arrowsize=XY_STREAM_ARROWSIZE, integration_direction="both",
            broken_streamlines=False, minlength=0.03, maxlength=20.0, zorder=3,
        )
    cbar = fig.colorbar(image, ax=ax, pad=0.025, extend="both")
    cbar.set_label(r"Plasma $\beta=P_{\rm gas}/P_{\rm mag}$")
    cbar.locator = LogLocator(base=10.0)
    cbar.update_ticks()
    finite = beta_grid[np.isfinite(beta_grid)]
    if finite.size and finite.min() <= 1.0 <= finite.max():
        ax.contour(
            X, Y, beta_grid, levels=[1.0],
            colors="white" if XY_DRAW_FIELDLINES else "black",
            linewidths=0.75, linestyles="--" if XY_DRAW_FIELDLINES else "-", zorder=4,
        )
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.set_aspect("equal", adjustable="box")
    xy_axis_format(ax, (xmin, xmax, ymin, ymax))
    title = (
        r"Plasma $\beta$ and projected magnetic field lines: x-y"
        if XY_DRAW_FIELDLINES else r"Plasma $\beta$: x-y"
    )
    time_kyr = data["time_code"]*XY_TIME_UNIT_KYR
    ax.set_title(
        f"{title}\nstep={step:05d}, t={time_kyr:.1f} kyr, "
        f"{range_name}, {XY_RUN_LABEL}"
    )
    fig.tight_layout()
    tag = "with_Blines" if XY_DRAW_FIELDLINES else "no_Blines"
    output_path = output_dir/f"xy_plasma_beta_{tag}_{step:05d}_{suffix}.png"
    if XY_SAVE_IMAGES:
        fig.savefig(output_path, dpi=XY_SAVE_DPI, bbox_inches="tight", facecolor="white")
        print(f"[SAVED] {output_path}")
    if XY_SHOW_IMAGES:
        plt.show()
    plt.close(fig)
    del X, Y, beta_grid, log_beta, beta
    if bx_grid is not None:
        del bx_grid, by_grid
    return output_path


XY_BETA_LIMITS = xy_beta_limits()
print(f"[XY] beta range = {XY_BETA_LIMITS[0]:.6e} -- {XY_BETA_LIMITS[1]:.6e}")
XY_FRACTION = None if XY_PLOT_MODE == "full" else XY_ZOOM_FRACTION
XY_SUFFIX = "full" if XY_PLOT_MODE == "full" else "zoom_1over10"
XY_OUTPUT_DIR = XY_OUTPUT_DIRS[XY_PLOT_MODE]
for xy_index, xy_step in enumerate(XY_STEPS, 1):
    print(f"[XY {xy_index:3d}/{len(XY_STEPS):3d}] step={xy_step:05d}")
    xy_data = xy_load_slice(xy_step)
    xy_plot(xy_data, XY_FRACTION, XY_OUTPUT_DIR, XY_SUFFIX, *XY_BETA_LIMITS)
    del xy_data
    gc.collect()

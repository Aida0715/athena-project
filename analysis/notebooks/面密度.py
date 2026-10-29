# -*- coding: utf-8 -*-
"""全VTK時刻の面密度 Sigma=integral rho dz と円環平均プロファイル。

Q値.pyと同じく既定では中心から |z|<=2500 AU を積分（Noneなら全高さ）。
SMRの各ブロックを順次読み、投影セルと画素の重なり面積で保存的に投影。
円環プロファイルは投影セルを小矩形に分割して質量を配分する近似。
マップ解像度・full/zoom設定には依存しない。完全な円環のみ表示する。
横軸表示範囲は半径×密度グラフ.pyと同じ球半径セル中心の最小・最大値。
ただし面密度の物理的な半径は円筒半径R。単位はmetadata JSONから取得。
磁力線.pyと同じくマップ座標はAUを10のべきで表示。格子線は描画しない。

DATADIRを変更して実行。必要環境: numpy, matplotlib, pyvista。
"""
from pathlib import Path
from collections import defaultdict
import json
import re

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter

# %% ユーザー設定
DATADIR = Path('~/athena-project/results/〇〇').expanduser()
VTK_GLOB = 'Toyouchi.block*.out2.*.vtk'
MAP_MODE = 'both'             # 'full', 'zoom', 'both'
ZOOM_FRACTION = 0.1           # 各軸の表示幅を全領域の1/10にする
MAP_RESOLUTION = 700          # 各軸の画素数（最細格子への全域展開はしない）
Z_HALF_WIDTH_AU = 2500.0       # None: 全z領域を積分
CENTER_CODE = (0.0, 0.0, 0.0) # 密度プロファイルと同じ中心
N_RADIAL_BINS = 100           # 密度プロファイルと同じ対数ビン数
ANNULUS_SUBDIVISIONS = 8      # 投影セルを各辺8分割して円環への面積配分を近似
SIGMA_LIMITS = None           # 色の固定範囲 (min, max) [Msun/AU^2]、None: 時刻ごと
OVERWRITE = False
DPI = 200
AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
KYR_CGS = 3.15576e10


def metadata_units(directory):
    preferred = directory / 'Toyouchi.metadata.json'
    paths = [preferred] if preferred.is_file() else sorted(directory.glob('*.metadata.json'))
    if len(paths) != 1:
        raise FileNotFoundError(f'Expected one metadata JSON in {directory}, found {len(paths)}')
    meta = json.loads(paths[0].read_text(encoding='utf-8'))
    units = meta['units']
    length = float(units['length_cgs']) / AU_CGS
    mass = float(units['mass_cgs']) / MSUN_CGS
    time = float(units['time_cgs']) / KYR_CGS
    if not all(np.isfinite(v) and v > 0 for v in (length, mass, time)):
        raise ValueError('Invalid units in metadata')
    return length, mass / length**2, time


def collect_steps(directory):
    groups = defaultdict(list)
    for path in directory.glob(VTK_GLOB):
        match = re.search(r'\.(\d+)\.vtk$', path.name)
        if match:
            groups[int(match.group(1))].append(path)
    if not groups:
        raise FileNotFoundError(f'No VTK files: {directory / VTK_GLOB}')
    return {step: sorted(paths) for step, paths in sorted(groups.items())}


def vtk_time(path):
    with path.open('rb') as handle:
        header = handle.read(4096).decode('ascii', errors='ignore')
    match = re.search(r'\btime\s*=\s*([-+0-9.eEdD]+)', header, re.I)
    if not match:
        raise ValueError(f'VTK time not found: {path}')
    return float(match.group(1).replace('D', 'E').replace('d', 'e'))


def faces(grid):
    """Athena Cartesian VTK (x-fastest ordering), ghost_zones=false を想定。"""
    if isinstance(grid, pv.RectilinearGrid):
        axes = [np.asarray(grid.x), np.asarray(grid.y), np.asarray(grid.z)]
    elif isinstance(grid, pv.ImageData):
        axes = [grid.origin[d] + np.arange(grid.extent[2*d], grid.extent[2*d+1]+1)
                * grid.spacing[d] for d in range(3)]
    else:
        raise ValueError(f'Expected axis-aligned Cartesian VTK, got {type(grid).__name__}')
    axes = [np.asarray(a, dtype=float) for a in axes]
    if any(len(a) < 2 or not np.all(np.isfinite(a)) or np.any(np.diff(a) <= 0) for a in axes):
        raise ValueError('Invalid cell faces')
    if np.prod([len(a)-1 for a in axes]) != grid.n_cells:
        raise ValueError('VTK dimensions and cell count disagree')
    if 'vtkGhostType' in grid.cell_data and np.any(grid.cell_data['vtkGhostType'] != 0):
        raise ValueError('Use VTK with ghost_zones=false')
    return axes


def inspect_domain(paths, center):
    low, high = np.full(3, np.inf), np.full(3, -np.inf)
    rmin, rmax = np.inf, 0.0
    for path in paths:
        grid = pv.read(path)
        axes = faces(grid)
        low = np.minimum(low, [a[0] for a in axes])
        high = np.maximum(high, [a[-1] for a in axes])
        radius = np.linalg.norm(grid.cell_centers().points - center, axis=1)
        valid = radius[np.isfinite(radius) & (radius > 0)]
        if valid.size:
            rmin, rmax = min(rmin, valid.min()), max(rmax, valid.max())
        del grid
    if not 0 < rmin < rmax or np.any(center <= low) or np.any(center >= high):
        raise ValueError('Invalid radial extent or center outside domain')
    return low, high, rmin, rmax


def overlap_matrix(target_edges, source_edges):
    """各画素幅に対するsourceセルとの重なり幅の割合。"""
    overlap = np.maximum(0.0, np.minimum(target_edges[1:, None], source_edges[None, 1:])
                         - np.maximum(target_edges[:-1, None], source_edges[None, :-1]))
    return overlap / np.diff(target_edges)[:, None]


def project_column(target, xedges, yedges, xfaces, yfaces, column):
    # ブロックと交差する画素だけを処理し、メモリと行列演算量を抑える。
    ix = np.flatnonzero((xedges[:-1] < xfaces[-1]) & (xedges[1:] > xfaces[0]))
    iy = np.flatnonzero((yedges[:-1] < yfaces[-1]) & (yedges[1:] > yfaces[0]))
    if not ix.size or not iy.size:
        return
    xs, ys = slice(ix[0], ix[-1]+1), slice(iy[0], iy[-1]+1)
    wx = overlap_matrix(xedges[ix[0]:ix[-1]+2], xfaces)
    wy = overlap_matrix(yedges[iy[0]:iy[-1]+2], yfaces)
    target[ys, xs] += wy @ column @ wx.T


def annular_mass(column, xfaces, yfaces, center, radial_edges, subdivisions):
    mass = np.zeros(len(radial_edges)-1)
    area = np.zeros_like(mass)
    dx, dy = np.diff(xfaces), np.diff(yfaces)
    patch_area = dy[:, None]*dx[None, :] / subdivisions**2
    for a in range(subdivisions):
        x = xfaces[:-1] + (a+0.5)/subdivisions*dx - center[0]
        for b in range(subdivisions):
            y = yfaces[:-1] + (b+0.5)/subdivisions*dy - center[1]
            radius = np.hypot(y[:, None], x[None, :])
            index = np.searchsorted(radial_edges, radius.ravel(), side='right')-1
            valid = (index >= 0) & (index < len(mass))
            mass += np.bincount(index[valid], weights=(column*patch_area).ravel()[valid], minlength=len(mass))
            area += np.bincount(index[valid], weights=patch_area.ravel()[valid], minlength=len(mass))
    return mass, area


def aggregate(paths, center, length_au, modes):
    low, high, rmin, rmax = inspect_domain(paths, center)
    edges = np.geomspace(rmin, np.nextafter(rmax, np.inf), N_RADIAL_BINS+1)
    complete_radius = np.min(np.r_[center[:2]-low[:2], high[:2]-center[:2]])
    maps = {}
    for mode in modes:
        if mode == 'full':
            lo, hi = low[:2], high[:2]
        else:
            half = (high[:2]-low[:2])*ZOOM_FRACTION/2
            lo, hi = np.maximum(low[:2], center[:2]-half), np.minimum(high[:2], center[:2]+half)
        xe = np.linspace(lo[0], hi[0], MAP_RESOLUTION+1)
        ye = np.linspace(lo[1], hi[1], MAP_RESOLUTION+1)
        maps[mode] = (xe, ye, np.zeros((MAP_RESOLUTION, MAP_RESOLUTION)))
    half_z = np.inf if Z_HALF_WIDTH_AU is None else Z_HALF_WIDTH_AU/length_au
    zlo, zhi = max(low[2], center[2]-half_z), min(high[2], center[2]+half_z)
    mass, area = np.zeros(N_RADIAL_BINS), np.zeros(N_RADIAL_BINS)
    for path in paths:
        grid = pv.read(path)
        xf, yf, zf = faces(grid)
        depth = np.maximum(0, np.minimum(zf[1:], zhi)-np.maximum(zf[:-1], zlo))
        if not np.any(depth):
            continue
        name = next((n for n in ('rho', 'dens', 'density', 'prim_dens', 'prim_density') if n in grid.cell_data), None)
        if name is None:
            raise KeyError(f'{path}: density cell data missing')
        rho = np.asarray(grid.cell_data[name], dtype=float).reshape(len(zf)-1, len(yf)-1, len(xf)-1)
        if np.any(~np.isfinite(rho)) or np.any(rho < 0):
            raise ValueError(f'{path}: negative or non-finite density')
        column = np.sum(rho*depth[:, None, None], axis=0)
        for xe, ye, target in maps.values():
            project_column(target, xe, ye, xf, yf, column)
        block_mass, block_area = annular_mass(column, xf, yf, center, edges, ANNULUS_SUBDIVISIONS)
        mass += block_mass
        # Fraction of the selected z interval avoids counting stacked blocks repeatedly.
        area += block_area*depth.sum()/(zhi-zlo)
        del grid, rho
    ann_area = np.pi*np.diff(edges**2)
    sigma = mass/ann_area
    # No fictitious zero points for unsampled bins; no incomplete outer annuli.
    sigma[(area == 0) | (edges[1:] > complete_radius)] = np.nan
    return maps, edges, sigma, (rmin, rmax), (zlo, zhi)


def draw_map(path, step, time, mode, data, length, sigma_unit, zlabel):
    xe, ye, column = data
    sigma = column*sigma_unit
    positive = sigma[np.isfinite(sigma) & (sigma > 0)]
    if not positive.size:
        raise ValueError(f'No positive surface density: step {step}, {mode}')
    limits = SIGMA_LIMITS or (positive.min(), positive.max())
    vmin, vmax = limits
    if vmin == vmax:
        vmin, vmax = vmin/1.01, vmax*1.01
    extent = np.array([xe[0], xe[-1], ye[0], ye[-1]])*length
    exponent = int(np.floor(np.log10(np.max(np.abs(extent)))))
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    im = ax.imshow(np.ma.masked_less_equal(sigma, 0), origin='lower', extent=extent,
                   norm=LogNorm(vmin, vmax), cmap='inferno', interpolation='nearest', aspect='equal')
    formatter = FuncFormatter(lambda value, position: f'{value/10.0**exponent:g}')
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel(rf'$x\ [10^{{{exponent}}}\ {{\rm AU}}]$')
    ax.set_ylabel(rf'$y\ [10^{{{exponent}}}\ {{\rm AU}}]$')
    ax.grid(False, which='both')
    ax.set_title(f'Surface density: {mode}\nstep={step:05d}, t={time:.1f} kyr; {zlabel}')
    fig.colorbar(im, ax=ax, label=r'$\Sigma$ [$M_\odot\,{\rm AU}^{-2}$]')
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def main():
    if MAP_MODE not in ('full', 'zoom', 'both') or not 0 < ZOOM_FRACTION <= 1:
        raise ValueError('Invalid MAP_MODE / ZOOM_FRACTION')
    if min(MAP_RESOLUTION, N_RADIAL_BINS, ANNULUS_SUBDIVISIONS) < 1:
        raise ValueError('Resolution and bin counts must be positive')
    if Z_HALF_WIDTH_AU is not None and Z_HALF_WIDTH_AU <= 0:
        raise ValueError('Z_HALF_WIDTH_AU must be positive or None')
    if SIGMA_LIMITS is not None and not 0 < SIGMA_LIMITS[0] < SIGMA_LIMITS[1]:
        raise ValueError('SIGMA_LIMITS must satisfy 0 < min < max')
    center = np.asarray(CENTER_CODE, dtype=float)
    if center.shape != (3,) or not np.all(np.isfinite(center)):
        raise ValueError('Invalid CENTER_CODE')
    directory = DATADIR.expanduser().resolve()
    length, sigma_unit, time_unit = metadata_units(directory)
    root = directory / '面密度解析'  # VTKファイルと同じディレクトリに置く親フォルダ
    modes = ('full', 'zoom') if MAP_MODE == 'both' else (MAP_MODE,)
    for folder in ('profiles', 'tables') + modes:
        (root/folder).mkdir(parents=True, exist_ok=True)
    print(f'Output: {root}; map mode={MAP_MODE}; z half-width={Z_HALF_WIDTH_AU} AU')
    for step, paths in collect_steps(directory).items():
        profile_path = root/'profiles'/f'Sigma_profile_{step:05d}.png'
        csv_path = root/'tables'/f'Sigma_profile_{step:05d}.csv'
        map_paths = {mode: root/mode/f'Sigma_{mode}_{step:05d}.png' for mode in modes}
        expected = [profile_path, csv_path, *map_paths.values()]
        if not OVERWRITE and all(p.is_file() for p in expected):
            print(f'[SKIP] step={step:05d}')
            continue
        time = vtk_time(paths[0])*time_unit
        maps, edges, sigma, extent, zrange = aggregate(paths, center, length, modes)
        radius = np.sqrt(edges[:-1]*edges[1:])*length
        sigma *= sigma_unit
        zlabel = f'z=[{zrange[0]*length:g}, {zrange[1]*length:g}] AU'
        if OVERWRITE or not profile_path.is_file():
            fig, ax = plt.subplots(figsize=(7.2, 5.2))
            ax.plot(radius, np.where(sigma > 0, sigma, np.nan), color='black', lw=1.6)
            ax.set(xscale='log', yscale='log', xlim=np.asarray(extent)*length,
                   xlabel=r'Radius $R$ [AU]', ylabel=r'Surface density $\Sigma$ [$M_\odot\,{\rm AU}^{-2}$]',
                   title=f'Azimuthally averaged surface density\nstep={step:05d}, t={time:.1f} kyr; {zlabel}')
            ax.grid(True, which='both', linestyle=':', alpha=.4)
            fig.tight_layout()
            fig.savefig(profile_path, dpi=DPI)
            plt.close(fig)
        if OVERWRITE or not csv_path.is_file():
            np.savetxt(csv_path, np.column_stack((edges[:-1]*length, edges[1:]*length, radius, sigma)),
                       delimiter=',', header='R_inner_AU,R_outer_AU,R_AU,Sigma_Msun_AU-2', comments='')
        for mode, path in map_paths.items():
            if OVERWRITE or not path.is_file():
                draw_map(path, step, time, mode, maps[mode], length, sigma_unit, zlabel)
        print(f'[OK] step={step:05d}; {len(paths)} blocks; profile + {", ".join(modes)} maps')


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
"""Toyouchi VTK の球平均加速度・動径速度。

必要環境: numpy, matplotlib, pyvista。datadir を変更して実行する。
MeshBlock ごとに処理し、セル体積（または質量）で重み付けする。
領域は VTK bounds から取得。完全な球殻が入る最大半径を既定にする。
中心は領域中央（中心天体がそこにある計算用）。必要なら CENTER を指定。

gr_nfw は現行 Toyouchi.cpp では softened radius s でのスカラー値。
真の動径成分は gr_nfw*r/s だが、epsilon は VTK に保存されていない。
したがって参考線としてのみ表示し、合計には含めない。
gr_star/gr_nfw は初期出力では未評価のゼロの場合がある。
合計は局所動径運動の右辺の平均であり、平均速度の時間微分ではない。
"""
# %% 設定
from pathlib import Path
import re
import warnings

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
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

AU_CM = 1.495978707e13

DATADIR = Path('~/athena-project/results/〇〇').expanduser()
OUTPUT_DIR = Path('./radial_acceleration_profiles')
VTK_PATTERN = '*.block*.vtk'
STEPS = None                 # None: 全ステップ、例: [0, 25, 50, 100]
N_BINS = 60
CENTER = None                # None: 領域中央、または (x, y, z) [code]
WEIGHT = 'volume'            # 'volume' / 'mass'
FULL_SHELLS_ONLY = True      # False: 隅まで含む（外側は部分球殻）

METADATA = load_toyouchi_metadata(DATADIR)
RUN_UNITS = run_units(METADATA)
LENGTH_UNIT_AU = RUN_UNITS['length_cgs'] / AU_CM
VELOCITY_UNIT_KMS = RUN_UNITS['velocity_cgs'] / 1.0e5
ACCELERATION_UNIT_CGS = RUN_UNITS['acceleration_cgs']

VECTORS = {'self': 'gself', 'pressure': 'apres', 'magnetic': 'amag'}
LABELS = {
    'self': 'Gas self-gravity', 'star': 'Central star',
    'nfw': 'NFW stored scalar (reference)', 'pressure': 'Pressure',
    'rotation': 'Azimuthal motion', 'polar': 'Polar motion',
    'magnetic': 'Magnetic', 'sum': 'Sum excluding NFW',
}


def collect_steps(directory):
    groups = {}
    for path in directory.glob(VTK_PATTERN):
        match = re.search(r'\.(\d+)\.vtk$', path.name)
        if match:
            step = int(match.group(1))
            if STEPS is None or step in STEPS:
                groups.setdefault(step, []).append(path)
    if not groups:
        raise FileNotFoundError(f'No matching VTK files: {directory / VTK_PATTERN}')
    return {step: sorted(paths) for step, paths in sorted(groups.items())}


def load(path):
    grid = pv.read(path)
    required = ['rho', 'vel', 'Bcc', 'gr_star', 'gr_nfw']
    required += [f'{prefix}_{axis}' for prefix in VECTORS.values() for axis in 'xyz']
    missing = [name for name in required if name not in grid.cell_data]
    if missing:
        raise KeyError(f'{path}: missing cell variables {missing}')
    return grid


def geometry(paths):
    bounds = []
    for path in paths:
        bounds.append(np.asarray(pv.read(path).bounds).reshape(3, 2))
    bounds = np.asarray(bounds)
    lo, hi = bounds[:, :, 0].min(axis=0), bounds[:, :, 1].max(axis=0)
    center = (lo + hi) / 2 if CENTER is None else np.asarray(CENTER, dtype=float)
    if center.shape != (3,) or np.any(center <= lo) or np.any(center >= hi):
        raise ValueError('CENTER must lie strictly inside the domain')
    if FULL_SHELLS_ONLY:
        rmax = np.min(np.r_[center - lo, hi - center])
    else:
        rmax = np.linalg.norm(np.maximum(center - lo, hi - center))
    rmin = np.inf
    for path in paths:
        radius = np.linalg.norm(pv.read(path).cell_centers().points - center, axis=1)
        positive = radius[(radius > 0) & (radius <= rmax)]
        if positive.size:
            rmin = min(rmin, positive.min())
    if not 0 < rmin < rmax:
        raise ValueError('Insufficient nonzero radial extent')
    edges = np.geomspace(rmin, rmax, N_BINS + 1)
    edges[0] = np.nextafter(edges[0], 0)
    print(f'Domain: {lo} .. {hi}; center={center}; radius={rmin:g} .. {rmax:g} [code]')
    return center, edges


def cell_quantities(grid, center):
    data = grid.cell_data
    xyz = grid.cell_centers().points - center
    r = np.linalg.norm(xyz, axis=1)
    unit = np.divide(xyz, r[:, None], out=np.zeros_like(xyz), where=r[:, None] > 0)
    velocity = np.asarray(data['vel'])
    vr = np.sum(velocity * unit, axis=1)
    R = np.linalg.norm(xyz[:, :2], axis=1)
    vphi = np.divide(-xyz[:, 1]*velocity[:, 0] + xyz[:, 0]*velocity[:, 1],
                     R, out=np.zeros_like(r), where=R > 0)
    # 軸上では方位方向を定義できないため、横方向運動は polar に入れる。
    vtheta2 = np.maximum(np.sum(velocity**2, axis=1) - vr**2 - vphi**2, 0)
    q = {'vr': vr, 'star': np.asarray(data['gr_star']), 'nfw': np.asarray(data['gr_nfw'])}
    for key, prefix in VECTORS.items():
        vector = np.column_stack([data[f'{prefix}_{axis}'] for axis in 'xyz'])
        q[key] = np.sum(vector * unit, axis=1)
    q['rotation'] = np.divide(vphi**2, r, out=np.zeros_like(r), where=r > 0)
    q['polar'] = np.divide(vtheta2, r, out=np.zeros_like(r), where=r > 0)
    q['sum'] = sum(q[key] for key in ('self', 'star', 'pressure', 'rotation', 'polar', 'magnetic'))
    volume = np.asarray(grid.compute_cell_sizes(length=False, area=False, volume=True).cell_data['Volume'])
    rho = np.asarray(data['rho'])
    if not np.all(np.isfinite(rho) & (rho > 0) & np.isfinite(volume) & (volume > 0)):
        raise ValueError('Invalid density or cell volume')
    if not all(np.all(np.isfinite(value)) for value in q.values()):
        raise ValueError('Non-finite acceleration or velocity')
    magnetic = np.asarray(data['Bcc'])
    if not np.all(np.isfinite(magnetic)):
        raise ValueError('Non-finite Bcc')
    has_magnetic = np.any(magnetic != 0) or any(np.any(data[f'amag_{axis}'] != 0) for axis in 'xyz')
    return r, q, volume if WEIGHT == 'volume' else volume*rho, has_magnetic


def analyze(paths, center, edges):
    totals, denominator = {}, np.zeros(N_BINS)
    counts = np.zeros(N_BINS, dtype=int)
    magnetic = False
    for path in paths:
        r, q, weights, block_magnetic = cell_quantities(load(path), center)
        magnetic |= block_magnetic
        index = np.searchsorted(edges, r, side='right') - 1
        index[r == edges[-1]] = N_BINS - 1
        valid = (r > 0) & (index >= 0) & (index < N_BINS)
        bins, w = index[valid], weights[valid]
        denominator += np.bincount(bins, weights=w, minlength=N_BINS)
        counts += np.bincount(bins, minlength=N_BINS)
        for key, values in q.items():
            totals.setdefault(key, np.zeros(N_BINS))
            totals[key] += np.bincount(bins, weights=w*values[valid], minlength=N_BINS)
    profiles = {key: np.divide(value, denominator, out=np.full(N_BINS, np.nan), where=denominator > 0)
                for key, value in totals.items()}
    return profiles, counts, denominator, magnetic


def save_plot(step, edges, profiles, magnetic):
    radius_au = np.sqrt(edges[:-1]*edges[1:]) * LENGTH_UNIT_AU
    acceleration_cgs = {
        key: value * ACCELERATION_UNIT_CGS
        for key, value in profiles.items()
        if key != 'vr'
    }
    radial_velocity_kms = profiles['vr'] * VELOCITY_UNIT_KMS
    fig, (ax, av) = plt.subplots(2, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    keys = ['self', 'star', 'nfw', 'pressure', 'rotation', 'polar']
    if magnetic:
        keys.append('magnetic')
    keys.append('sum')
    for key in keys:
        kwargs = {'color': 'black', 'lw': 2.0} if key == 'sum' else {'lw': 1.3}
        ax.plot(radius_au, acceleration_cgs[key], '--' if key == 'nfw' else '-', label=LABELS[key], **kwargs)
    values = np.concatenate([np.abs(acceleration_cgs[key]) for key in keys])
    finite = values[np.isfinite(values) & (values > 0)]
    threshold = finite.max()*1e-4 if finite.size else 1.0
    ax.set_yscale('symlog', linthresh=threshold)
    ax.set_ylabel(r'Radial acceleration [cm s$^{-2}$] (+ outward)')
    ax.legend(fontsize=9, ncol=2)
    av.plot(radius_au, radial_velocity_kms, color='black', label='Mean radial velocity')
    av.set_ylabel(r'Radial velocity [km s$^{-1}$] (+ outward)')
    av.set_xlabel('Spherical radius [AU]')
    av.legend()
    for panel in (ax, av):
        panel.set_xscale('log')
        panel.axhline(0, color='0.5', lw=.7)
        panel.grid(True, which='both', alpha=.25)
    region = 'complete shells' if FULL_SHELLS_ONLY else 'outer shells partially sampled'
    fig.suptitle(f'Step {step:05d} | {WEIGHT}-weighted | {region}\nNFW scalar shown separately; excluded from sum')
    fig.savefig(OUTPUT_DIR / f'radial_acceleration_{step:05d}.png', dpi=160)
    plt.close(fig)


def main():
    if WEIGHT not in ('volume', 'mass') or N_BINS < 2:
        raise ValueError('Invalid WEIGHT or N_BINS')
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    warnings.warn(f"Physical-unit conversion loaded from {METADATA['_path']}. "
                  'NFW stored scalar is a reference, excluded from sum. '
                  'Initial gr_star/gr_nfw may be uninitialized zeros.')
    for step, paths in collect_steps(DATADIR).items():
        print(f'Step {step:05d}: {len(paths)} blocks')
        center, edges = geometry(paths)
        profiles, counts, weights, magnetic = analyze(paths, center, edges)
        save_plot(step, edges, profiles, magnetic)
        acceleration_keys = [key for key in profiles if key != 'vr']
        table = np.column_stack([
            edges[:-1] * LENGTH_UNIT_AU,
            edges[1:] * LENGTH_UNIT_AU,
            counts,
            weights,
            profiles['vr'] * VELOCITY_UNIT_KMS,
        ] + [profiles[key] * ACCELERATION_UNIT_CGS for key in acceleration_keys])
        np.savetxt(OUTPUT_DIR / f'radial_acceleration_{step:05d}.csv', table, delimiter=',',
                   header=','.join([
                       'r_left_au', 'r_right_au', 'cell_count', 'weight_sum_code',
                       'vr_km_s',
                   ] + [f'{key}_cm_s-2' for key in acceleration_keys]), comments='')
        print(f'Saved PNG/CSV; empty shells={np.count_nonzero(counts == 0)}; magnetic={magnetic}')


if __name__ == '__main__':
    main()

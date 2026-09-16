# -*- coding: utf-8 -*-
# Converted from checks_dt_cfl.ipynb

# %% [markdown] cell 1
# # Toyouchi：時間刻みとCFL制限セルの確認
#
# 上から順に実行。最初の設定セルの `DATA_DIR` と **今回実際に使用した**音速・CFL・単位を確認してください。
# 必要なパッケージは numpy / pandas / matplotlib のみ（必要ならJupyterで `%pip install numpy pandas matplotlib`）。
# 出力は `DATA_DIR/checks/`。既存の同名解析結果は更新します。計算データは変更しません。
#
# 対象：3次元Cartesian、等温、磁場ゼロ、ghostなし・slice/sumなしのAthena++標準binary VTK（ブロック別）。
# 各ブロックを順次読み、全領域を一様な最細格子に展開しません。
#
# **ログのdtは実測、VTKのCFL刻みは保存状態からの再計算です。** 保存間隔中の一過性イベントの場所は特定できません。
# 実際のdtには前ステップの2倍までの増加制限、終端時刻、拡散・ユーザー条件等も作用します。
# VTKは単精度で時刻表示の桁も少ないため、照合にはcycleの完全一致を使用します。

# %% cell 2
from pathlib import Path
import re, json, heapq, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

DATA_DIR = Path("/work/beta/aida/results/Toyouchi-test55").expanduser().resolve()
# 複数指定時はrestartの古い順。NoneならDATA_DIR直下の*.logを名前順に読む。
# 同じcycle/timeの重複は後のファイルを採用。別の計算を混ぜないこと。
LOG_FILES = None
VTK_GLOB = "*.block*.out2.*.vtk"  # rglobでMPIサブディレクトリも検索
EXPECTED_BLOCKS = None  # SMRの総MeshBlock数が分かれば必ず指定。Noneは最多ファイル数を仮採用。
SNAPSHOT_STRIDE = 1
TOP_K = 10
SMALL_DT = 1e-6  # code time
CFL = 0.3
CS = 0.707      # iso_sound_speed。実行時の入力に合わせる
ROOT_DX = np.array([448/40]*3)
L_UNIT_CM, T_UNIT_S = 7.03e15, 3.61e10
AU_CM = 1.495978707e13
L_AU = L_UNIT_CM/AU_CM
T_YR = T_UNIT_S/(365.25*86400)
SINK_RADIUS_AU = 1000.0
# 2つ目のログ推移を重ねたい場合：例 Path('/work/beta/aida/results/Toyouchi-test54/run.log')
COMPARISON_LOG = None

if not DATA_DIR.is_dir():
    raise FileNotFoundError(DATA_DIR)
assert CFL > 0 and CS > 0 and TOP_K >= 1 and SNAPSHOT_STRIDE >= 1
OUT = DATA_DIR / "checks"
OUT.mkdir(exist_ok=True)
(OUT/'settings.json').write_text(json.dumps({
    'data_dir': str(DATA_DIR), 'logs': None if LOG_FILES is None else list(map(str, LOG_FILES)),
    'vtk_glob': VTK_GLOB, 'expected_blocks': EXPECTED_BLOCKS, 'cfl': CFL, 'cs': CS,
    'root_dx': ROOT_DX.tolist(), 'length_unit_cm': L_UNIT_CM, 'time_unit_s': T_UNIT_S,
    'sink_radius_au': SINK_RADIUS_AU, 'small_dt': SMALL_DT,
    'snapshot_stride': SNAPSHOT_STRIDE, 'top_k': TOP_K,
}, indent=2), encoding='utf-8')
print('保存先:', OUT)

def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(OUT/name, dpi=180, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# %% [markdown] cell 3
# ## 1. ログ：dtの時間変化・ステップ消費
#
# 時間軸とcycle軸の両方を描画します。cycle軸なら、物理時間上では細く見える停滞も分かります。
# 各ログ内でcycle/timeが巻き戻った場合、restart区間を分けて線を引きます。
# `small_dt_events.csv` は閾値未満の記録、`dt_bins.csv` は刻み別の記録数・dt合計です。
# 毎cycleのログがある場合のみ、記録数がステップ数に相当します。CPU時間の配分はこのログだけでは分かりません。

# %% cell 4
NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
LOG_RE = re.compile(r"cycle\s*=\s*(\d+)\s+time\s*=\s*("+NUM+r")\s+dt\s*=\s*("+NUM+r")")
def number(s):
    return float(s.replace('D','E').replace('d','e'))

def read_logs(paths):
    rows = []
    for path in paths:
        previous, segment = None, 0
        with Path(path).open(errors='replace') as f:
            for line_no, line in enumerate(f, 1):
                m = LOG_RE.search(line)
                if not m:
                    continue
                cycle, time, dt = int(m[1]), number(m[2]), number(m[3])
                if previous and (cycle < previous[0] or time < previous[1]):
                    segment += 1
                previous = cycle, time
                rows.append(dict(cycle=cycle, time=time, dt=dt, source=str(path),
                                 line=line_no, segment=f'{path}:{segment}'))
    if not rows:
        raise ValueError('cycle/time/dtを含むログ行がありません。LOG_FILESを指定してください。')
    df = pd.DataFrame(rows).drop_duplicates(['cycle','time'], keep='last')
    if (~np.isfinite(df[['time','dt']].to_numpy())).any() or (df.dt <= 0).any():
        raise ValueError('非有限または非正のdtを検出。ログを確認してください。')
    return df.reset_index(drop=True)

paths = sorted(DATA_DIR.glob('*.log')) if LOG_FILES is None else [Path(p).expanduser() for p in LOG_FILES]
log = read_logs(paths)
log.to_csv(OUT/'dt_history.csv', index=False)
log[log.dt < SMALL_DT].to_csv(OUT/'small_dt_events.csv', index=False)
if log.duplicated('cycle', keep=False).any():
    warnings.warn('同じcycleに異なる時刻があります。restart分岐/別計算の混在を確認。VTK照合時は曖昧なcycleを除外します。')
gaps = log.groupby('segment', sort=False).cycle.diff().dropna()
print('ログ行数:', len(log), '毎cycle以外の間隔数:', int((gaps != 1).sum()))
display(log.nsmallest(10, 'dt')[['cycle','time','dt','source','line']])

fig, axes = plt.subplots(2, 1, figsize=(11, 7))
for _, part in log.groupby('segment', sort=False):
    axes[0].semilogy(part.time*T_YR/1000, part.dt, lw=.65, color='C0')
    axes[1].semilogy(part.cycle, part.dt, lw=.65, color='C0')
if COMPARISON_LOG is not None:
    other = read_logs([COMPARISON_LOG])
    for _, part in other.groupby('segment', sort=False):
        axes[0].semilogy(part.time*T_YR/1000, part.dt, color='C1', lw=.8)
    axes[0].plot([], [], color='C0', label=DATA_DIR.name)
    axes[0].plot([], [], color='C1', label=Path(COMPARISON_LOG).parent.name)
    axes[0].legend()
for ax in axes:
    ax.axhline(SMALL_DT, ls='--', color='gray')
    ax.set_ylabel('dt [code]'); ax.grid(alpha=.25)
axes[0].set_xlabel('time [kyr]'); axes[1].set_xlabel('cycle')
savefig(fig, 'dt_history.png')

bins = np.r_[-np.inf, 10.**np.arange(-10, -1), np.inf]
stats = log.assign(dt_bin=pd.cut(log.dt, bins)).groupby('dt_bin', observed=True).agg(
    records=('dt','size'), sum_dt=('dt','sum'))
stats['record_fraction'] = stats.records/len(log)
stats.to_csv(OUT/'dt_bins.csv')
display(stats)
fig, ax = plt.subplots(figsize=(11, 4))
ax.bar(stats.index.astype(str), stats.records)
ax.set(ylabel='logged records', xlabel='dt bin [code]')
ax.tick_params(axis='x', rotation=45)
savefig(fig, 'dt_bins.png')

# %% [markdown] cell 5
# ## 2. VTKをブロック単位で読み、CFL制限候補を抽出
#
# `dt_CFL = CFL × min(dx/(abs(vx)+cs), dy/(abs(vy)+cs), dz/(abs(vz)+cs))`。
# 保存されたleaf blockの全セルを調べ、各時刻の上位TOP_Kセルを保存します。i,j,kはghostを含まない0始まり。
# 磁場が非ゼロのデータは、等温HDの式を誤適用しないよう処理を中止します。
# 非有限値・非正密度、欠損/重複ブロック、時刻不一致、書き込み途中のファイルはスナップショット全体を除外します。
# SMR専用の整合性検査としてブロックID集合の一致も確認します。AMRの出力にはそのまま使わないでください。

# %% cell 6
def read_vtk(path):
    # Athena++ src/outputs/vtk.cppのlegacy binary float形式に限定。
    with Path(path).open('rb') as f:
        def line():
            while True:
                raw = f.readline()
                if not raw:
                    return ''
                s = raw.decode('ascii').strip()
                if s:
                    return s
        def floats(n):
            raw = f.read(4*n)
            if len(raw) != 4*n:
                raise ValueError('未完了のbinary payload')
            return np.frombuffer(raw, dtype='>f4').astype(float)
        if not line().startswith('# vtk DataFile'):
            raise ValueError('VTKヘッダ不正')
        header = line()
        tm = re.search(r'time=('+NUM+')', header)
        cy = re.search(r'cycle=(\d+)', header)
        if tm is None or cy is None:
            raise ValueError('time/cycleなし')
        if line() != 'BINARY' or line() != 'DATASET RECTILINEAR_GRID':
            raise ValueError('binary rectilinear VTKのみ対応')
        dims_line = line().split()
        if dims_line[0] != 'DIMENSIONS':
            raise ValueError('DIMENSIONSなし')
        dims = tuple(map(int, dims_line[1:]))
        if len(dims) != 3 or min(dims) < 2:
            raise ValueError('3次元volume出力のみ対応')
        coords = []
        for axis, count in zip('XYZ', dims):
            spec = line().split()
            if spec != [axis+'_COORDINATES', str(count), 'float']:
                raise ValueError('座標形式不正')
            coords.append(floats(count))
        n = int(np.prod(np.array(dims)-1))
        if line().split() != ['CELL_DATA', str(n)]:
            raise ValueError('CELL_DATA数不正')
        fields = {}
        while True:
            spec = line().split()
            if not spec:
                break
            if len(spec) < 3 or spec[2] != 'float':
                raise ValueError('float fieldsのみ対応')
            kind, name = spec[:2]
            if kind == 'SCALARS':
                if len(spec) > 3 and spec[3] != '1':
                    raise ValueError('複数成分SCALARSは非対応')
                if not line().startswith('LOOKUP_TABLE '):
                    raise ValueError('LOOKUP_TABLEなし')
                fields[name] = floats(n)
            elif kind == 'VECTORS':
                fields[name] = floats(3*n).reshape(n, 3)
            else:
                raise ValueError('未知のfield形式: '+kind)
    return number(tm[1]), int(cy[1]), coords, fields

def field(fields, candidates):
    for key in candidates:
        if key in fields:
            return fields[key]
    raise ValueError(f'必要な変数{candidates}がありません。実際: {list(fields)}')

def block_candidates(path):
    time, cycle, xyz, fields = read_vtk(path)
    rho = field(fields, ['rho','dens','density'])
    vel = field(fields, ['vel','velocity'])
    if vel.shape != (len(rho), 3):
        raise ValueError('速度配列のshape不正')
    for key, value in fields.items():
        if not np.isfinite(value).all():
            raise ValueError(f'非有限値: {key}')
        if (key.lower().startswith('bcc') or key.lower() in ('b','b1','b2','b3','bfield','magnetic_field')) and np.any(value != 0):
            raise ValueError('非ゼロ磁場：HD式の対象外')
    if np.any(rho <= 0):
        raise ValueError('非正密度')
    widths = [np.diff(a) for a in xyz]
    if any(not np.isfinite(a).all() or np.any(a <= 0) for a in widths):
        raise ValueError('セル幅不正')
    nx, ny, nz = map(len, widths)
    k, j, i = np.indices((nz, ny, nx)).reshape(3, -1)
    width = np.column_stack([widths[0][i], widths[1][j], widths[2][k]])
    centers = [(a[:-1]+a[1:])/2 for a in xyz]
    pos = np.column_stack([centers[0][i], centers[1][j], centers[2][k]])
    direction_dt = CFL*width/(np.abs(vel)+CS)
    dt = direction_dt.min(axis=1)
    if not np.isfinite(dt).all() or np.any(dt <= 0):
        raise ValueError('CFL計算結果不正')
    take = np.argsort(dt, kind='stable')[:TOP_K]
    bid = int(re.search(r'\.block(\d+)\.', Path(path).name)[1])
    rows = []
    for q in take:
        direction = int(np.argmin(direction_dt[q]))
        radius = np.linalg.norm(pos[q])*L_AU
        levels = np.log2(ROOT_DX/width[q])
        if not np.allclose(levels, np.rint(levels), atol=2e-3) or np.ptp(levels) > 2e-3:
            raise ValueError('ROOT_DXとセル幅が整合しません')
        row = dict(time=time, cycle=cycle, block=bid, i=int(i[q]), j=int(j[q]), k=int(k[q]),
                   dt_cfl=dt[q], direction='xyz'[direction], level=int(round(levels[0])),
                   rho=rho[q], speed=np.linalg.norm(vel[q]), radius_au=radius,
                   inside_sink=bool(radius < SINK_RADIUS_AU),
                   block_edge_layers=int(min(i[q], nx-1-i[q], j[q], ny-1-j[q], k[q], nz-1-k[q])),
                   file=str(path))
        for a, label in enumerate('xyz'):
            row[label+'_au'] = pos[q,a]*L_AU
            row['d'+label] = width[q,a]
            row['v'+label] = vel[q,a]
        for key in ('gself_x','gself_y','gself_z','gr_star','gr_nfw'):
            if key in fields and fields[key].ndim == 1:
                row[key] = fields[key][q]
        rows.append(row)
    return rows, dict(rho_max=float(rho.max()), speed_max=float(np.linalg.norm(vel, axis=1).max()),
                      cells=len(rho))

# %% cell 7
groups = {}
for p in sorted(DATA_DIR.rglob(VTK_GLOB)):
    if OUT in p.parents:
        continue
    match = re.search(r'\.block(\d+)\.[^.]+\.(\d+)\.vtk$', p.name)
    if match:
        groups.setdefault(int(match[2]), []).append((int(match[1]), p))
if not groups:
    raise FileNotFoundError(f'VTKが見つかりません: {DATA_DIR}/{VTK_GLOB}')
reference = max(groups.values(), key=len)
expected_ids = {bid for bid, _ in reference}
expected_count = EXPECTED_BLOCKS if EXPECTED_BLOCKS is not None else len(reference)
if EXPECTED_BLOCKS is None:
    warnings.warn(f'総ブロック数を{expected_count}と推定。全出力で同じブロックが欠落している場合は検出不能です。実行ログと照合してください。')

all_rows, summaries, audit = [], [], []
selected = sorted(groups)[::SNAPSHOT_STRIDE]
for count, snap in enumerate(selected, 1):
    members = groups[snap]
    try:
        ids = [bid for bid, _ in members]
        if len(ids) != expected_count or len(set(ids)) != len(ids) or set(ids) != expected_ids:
            raise ValueError('ブロック数/ID集合の不一致または重複')
        best, times, cycles = [], set(), set()
        rho_max, speed_max, ncells = 0., 0., 0
        for bid, path in members:
            rows, block = block_candidates(path)
            times.add(rows[0]['time']); cycles.add(rows[0]['cycle'])
            best = sorted(best+rows, key=lambda r: (r['dt_cfl'],r['block'],r['k'],r['j'],r['i']))[:TOP_K]
            rho_max = max(rho_max, block['rho_max'])
            speed_max = max(speed_max, block['speed_max'])
            ncells += block['cells']
        if len(times) != 1 or len(cycles) != 1:
            raise ValueError('ブロック間でtime/cycle不一致：書き込み途中の可能性')
        for rank, row in enumerate(best, 1):
            row.update(snapshot=snap, rank=rank)
        all_rows.extend(best)
        summaries.append(dict(best[0], rho_max=rho_max, speed_max=speed_max, cells=ncells))
        audit.append(dict(snapshot=snap, status='ok', blocks=len(ids), reason=''))
    except (ValueError, OSError, UnicodeError, IndexError) as exc:
        audit.append(dict(snapshot=snap, status='skipped', blocks=len(members), reason=str(exc)))
        print('除外:', snap, str(exc))
    if count % 10 == 0 or count == len(selected):
        print(f'{count}/{len(selected)} snapshots processed')
pd.DataFrame(audit).to_csv(OUT/'snapshot_audit.csv', index=False)
(OUT/'snapshot_inventory.json').write_text(json.dumps({
    'expected_count_used': expected_count, 'reference_block_ids': sorted(expected_ids),
    'selected_snapshot_ids': selected,
}, indent=2), encoding='utf-8')
pd.DataFrame(all_rows).to_csv(OUT/'cfl_top_cells.csv', index=False)
snapshots = pd.DataFrame(summaries)
snapshots.to_csv(OUT/'cfl_snapshots.csv', index=False)
if snapshots.empty:
    raise ValueError('有効なsnapshotなし。snapshot_audit.csvを確認してください。')
snapshots = snapshots.sort_values(['time','cycle'])
print('有効snapshot:', len(snapshots))
display(snapshots.nsmallest(10, 'dt_cfl')[['cycle','time','dt_cfl','block','i','j','k','level','radius_au','rho','speed']])

# %% [markdown] cell 8
# ## 3. ログとの照合・制限候補の軌跡
#
# ログとVTKのcycleが一意に完全一致したものだけdtを比較します。
# 一致しても、増加制限や保存タイミング・精度などによって一致しない場合があります。
# `small_dt_snapshot_coverage.csv` の距離は「そのイベントの前後に保存状態があるか」の確認用で、場所の同定ではありません。
# 位置図は各保存時刻の最小CFLセルを表します。図の円/横線はシンク半径です。

# %% cell 9
unique_log = log[~log.duplicated('cycle', keep=False)][['cycle','time','dt']].rename(
    columns={'time':'log_time','dt':'dt_logged'})
matched = snapshots.merge(unique_log, on='cycle', how='left', validate='many_to_one')
matched['cfl_over_logged'] = matched.dt_cfl/matched.dt_logged
matched.to_csv(OUT/'cfl_log_comparison.csv', index=False)

events = log[log.dt < SMALL_DT].copy()
st = snapshots.time.to_numpy()
if not events.empty:
    et = events.time.to_numpy()
    right = np.searchsorted(st, et).clip(0, len(st)-1)
    left = (right-1).clip(0, len(st)-1)
    near = np.where(abs(st[left]-et) <= abs(st[right]-et), left, right)
    events['nearest_snapshot_time'] = st[near]
    events['nearest_snapshot_cycle'] = snapshots.cycle.to_numpy()[near]
    events['time_distance_code'] = abs(st[near]-et)
    events['same_cycle'] = events.cycle.to_numpy() == snapshots.cycle.to_numpy()[near]
events.to_csv(OUT/'small_dt_snapshot_coverage.csv', index=False)

fig, ax = plt.subplots(figsize=(11, 4))
for _, part in log.groupby('segment', sort=False):
    ax.semilogy(part.time*T_YR/1000, part.dt, color='gray', lw=.6)
ax.semilogy(snapshots.time*T_YR/1000, snapshots.dt_cfl, 'o', ms=3, label='CFL from saved VTK')
ax.set(xlabel='time [kyr]', ylabel='dt [code]'); ax.legend(); ax.grid(alpha=.25)
savefig(fig, 'dt_log_vs_snapshot.png')

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
t = snapshots.time*T_YR/1000
axes[0,0].plot(t, snapshots.radius_au, '.-')
axes[0,0].axhline(SINK_RADIUS_AU, color='k', ls='--')
axes[0,0].set_ylabel('limiting radius [au]')
axes[0,1].step(t, snapshots.level, where='mid'); axes[0,1].set_ylabel('limiting cell level')
axes[1,0].semilogy(t, snapshots.rho, '.-', label='limiting cell')
axes[1,0].semilogy(t, snapshots.rho_max, '--', label='domain maximum')
axes[1,0].set_ylabel('density [code]'); axes[1,0].legend()
axes[1,1].semilogy(t, np.maximum(snapshots.speed, 1e-30), '.-', label='limiting cell')
axes[1,1].semilogy(t, np.maximum(snapshots.speed_max, 1e-30), '--', label='domain maximum')
axes[1,1].set_ylabel('speed [code]'); axes[1,1].legend()
for ax in axes.flat:
    ax.set_xlabel('time [kyr]'); ax.grid(alpha=.25)
savefig(fig, 'cfl_cell_history.png')

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
color = np.log10(snapshots.dt_cfl)
for ax, horizontal, vertical in [(axes[0],'x','y'), (axes[1],'x','z')]:
    sc = ax.scatter(snapshots[horizontal+'_au'], snapshots[vertical+'_au'], c=color, s=20, cmap='viridis')
    ax.add_patch(plt.Circle((0,0), SINK_RADIUS_AU, fill=False, color='gray', ls='--'))
    ax.set(xlabel=horizontal+' [au]', ylabel=vertical+' [au]', aspect='equal')
    ax.grid(alpha=.25)
    fig.colorbar(sc, ax=ax, label='log10(dt_CFL [code])')
savefig(fig, 'cfl_cell_positions.png')
print('完了:', OUT)
print('極小dtの場所が保存時刻から分からない場合は、実行側でCFL最小セルを記録する追加診断が必要です。')


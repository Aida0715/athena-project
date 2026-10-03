# -*- coding: utf-8 -*-
"""Toyouchi: ログ・HST・3D Cartesian binary VTKの健全性確認。

Jupyter用: 設定セルを書き換え、各セルを上から実行。
IDEでは # %% ごとに実行、.ipynbでは各区切りの内容をコードセルへ貼り付ける。
ファイル全体を1つのJupyterセルへ貼り付けても実行できる。
依存: numpy, matplotlib。PyVista不要。VTKを1ブロックずつ読む。
出力: VTK_DIR/checks/。同名の解析出力は更新。入力データは変更しない。
全出力・全セルを検査するが、保存時刻間の異常は検出できない。
対象: ghostなし、slice/sumなし、静的細分化の標準Athena++ VTK。
表示・CSVは時間kyr、質量Msun、その他CGS。内部計算・入力設定はシミュレーション単位。
閾値は物理的な合否基準ではない。質量流量はCGSのg/s。
"""
# %% 設定
from pathlib import Path
import csv
import json
import re
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import PolyCollection

VTK_DIR = Path("~/athena-project/results/〇〇").expanduser()
HST_FILENAME = "Toyouchi.hst"
LOG_GLOB = "*.log"  # この計算のログのみ置く。各ファイルを独立した実行区間として扱う
VTK_GLOB = "*.block*.out2.*.vtk"  # MPIサブディレクトリも再帰検索
EXPECTED_BLOCKS = None  # 実行時の総MeshBlock数。未指定なら完全性は断定しない
EXPECTED_TLIM = 43.86
DFLOOR = 1.0e-6
SINK_RHO = 6.3e-2
SINK_RADIUS = 1000.0 * 1.495978707e13 / 7.03e15
SINK_CENTER = np.array([0.0, 0.0, 0.0])
CS = 0.707
CFL = 0.3
ETA_OHM = 0.1  # 一定Ohmic拡散、Hall/ADなしの場合の刻み目安
SMALL_DT = 1.0e-6
DT_DROP_RATIO = 0.2  # 前記録の20%未満を参考イベントとして保存
BALANCE_RTOL = 1.0e-5  # HST桁落ちも考慮した整合性検査用
SLICE_STEPS = None  # None: 初期・中間・最終。番号リストで追加指定可
SLICE_HALF_WIDTH = None  # None: 領域全体。code長でsink周辺のズーム指定可
ROOT_CELL_WIDTH = np.array([448/40]*3)  # 実行時の各軸の領域幅/nx。SMR level=0のセル幅
SINK_NO_OUTFLOW_CELLS = 1.0  # 実行時のsink_no_outflow_cellsと一致させる
BUFFER_VR_TOL_CGS = 1.e-6  # 外向き速度の判定許容値 [cm/s]。データ自体は丸めない
PROFILE_BINS = 80  # 初期密度プロファイルの球殻数（対数間隔、中心殻は0から）
SHOW_PLOTS = True  # Jupyterのセル出力にも図を表示（PNG保存は常に行う）

# %% 単位換算・物理量一覧・診断基準（最初の出力）
# 入力と内部計算はシミュレーション単位を保ち、表示・保存時に一度だけ換算する。
# 設定の時間・密度・長さ等も入力ファイルと同じシミュレーション単位。
M_UNIT_G = 4.0e33
L_UNIT_CM = 7.03e15
T_UNIT_S = 3.61e10
MSUN_G = 1.98847e33
KYR_S = 1000 * 365.25 * 86400
TIME_KYR = T_UNIT_S / KYR_S
MASS_MSUN = M_UNIT_G / MSUN_G
VELOCITY_CGS = L_UNIT_CM / T_UNIT_S
DENSITY_CGS = M_UNIT_G / L_UNIT_CM**3
ENERGY_CGS = M_UNIT_G * VELOCITY_CGS**2
B_GAUSS = np.sqrt(4*np.pi*DENSITY_CGS*VELOCITY_CGS**2)

# name -> (conversion factor, unit, description); CSV見出しにも単位を付ける。
QUANTITIES = {}

def register(names, factor, unit, description):
    for name in names.split():
        QUANTITIES[name] = (factor, unit, description)

register('time dt dt_wave_est dt_ohm_est', TIME_KYR, 'kyr', '時刻・時間刻み（dt_*_estは保存状態からの参考推定）')
register('mass', MASS_MSUN, 'Msun', '計算領域のガス質量（sink内を含む）')
register('Mstar', MASS_MSUN, 'Msun', '中心星質量')
register('Msink_gas sink_mass', MASS_MSUN, 'Msun', 'sink内に残るガス質量（HST / VTK積分）')
register('Mflux_cum', MASS_MSUN, 'Msun', '毎cycleで負の正味流入を0にした、中心星への累積加算質量')
register('Mreset_cum', MASS_MSUN, 'Msun', 'sinkリセットで除去した累積質量')
register('Mfloor_cum', MASS_MSUN, 'Msun', 'sinkリセットで追加した累積質量（EOS floor全体ではない）')
register('Mdot_flux', M_UNIT_G/T_UNIT_S, 'g/s', '負値切捨て後のsink正味流入率')
register('Mdot_reset', M_UNIT_G/T_UNIT_S, 'g/s', 'sinkリセットによる質量除去率')
register('Mdot_floor', M_UNIT_G/T_UNIT_S, 'g/s', 'sinkリセットによる質量追加率（M_dot_floorとも表記）')
register('1-mom 2-mom 3-mom', M_UNIT_G*VELOCITY_CGS, 'g cm/s', '各方向の領域積分運動量')
register('1-KE 2-KE 3-KE', ENERGY_CGS, 'erg', '各速度成分の領域積分運動エネルギー')
register('grav-E', ENERGY_CGS, 'erg', 'HSTで記録した重力エネルギー（全重力源を含むとは限らない）')
register('1-ME 2-ME 3-ME magnetic_energy', ENERGY_CGS, 'erg', '磁場エネルギー：HST各成分 / VTK全成分の和')
register('rho rho_min rho_max', DENSITY_CGS, 'g/cm^3', 'ガス密度・最小値・最大値')
register('vel v_max va_max va_max_outside Va_max_sink', VELOCITY_CGS, 'cm/s', '速度・最大速度・Alfven速度（全域 / sink外 / sink内）')
register('Bcc B_max Bmax_sink', B_GAUSS, 'G', 'セル中心磁場・全域最大強度・sink内最大強度')
register('phi', VELOCITY_CGS**2, 'cm^2/s^2', '保存された重力ポテンシャル')
register('gr_star gr_nfw gself_x gself_y gself_z apres_x apres_y apres_z amag_x amag_y amag_z', L_UNIT_CM/T_UNIT_S**2, 'cm/s^2', '中心星・NFW・自己重力・圧力勾配・磁場による加速度（有限値検査）')
register('star_identity', MASS_MSUN, 'Msun', 'ΔMstar−ΔMflux_cum：中心星更新の定義上の整合性')
register('sink_proxy', MASS_MSUN, 'Msun', 'ΔMsink−ΔMflux+ΔMreset−ΔMfloor：符号付きfluxがなく参考残差')
register('gas_unclosed', MASS_MSUN, 'Msun', 'ΔMgas+ΔMreset−ΔMfloor：外部境界flux等を含まない未閉合量')
register('floor_added reset_removed accreted', MASS_MSUN, 'Msun', 'HST区間始点からの累積追加・除去・降着質量の増分')
register('floor_to_accreted floor_fraction beta_min', 1., '1', '追加/降着質量比・floorセル割合・最小plasma beta（無次元）')
register('x y z', L_UNIT_CM, 'cm', 'セル中心の位置、断面図の座標')
register('volume div_proxy_volume', L_UNIT_CM**3, 'cm^3', '全セル体積 / 発散参考値の有効体積')
register('div_proxy_max div_proxy_rms', B_GAUSS/L_UNIT_CM, 'G/cm', 'セル中心Bの発散の最大絶対値・体積加重RMS（CT誤差ではない）')
register('div_proxy_sqvol', B_GAUSS**2*L_UNIT_CM, 'G^2 cm', '発散参考値の二乗体積積分（RMS集計用）')
register('cells floor_cells floor_cells_outside sink_target_cells nonfinite nonpositive_rho count blocks', 1., 'count', 'セル数・非有限成分数・非正密度セル数・ブロック数等')

register('radius radius_inner radius_outer', L_UNIT_CM, 'cm', '球殻の代表半径・内縁・外縁（sink中心から）')
register('rho_mean', DENSITY_CGS, 'g/cm^3', '最初の保存出力における球殻体積加重平均密度（sink内を含む）')
register('mdot_floor_to_flux', 1., '1', '瞬間のMdot_floor / Mdot_flux。分母0・非有限値は未定義')
register('mdot_reset_to_floor', 1., '1', '瞬間のMdot_reset / Mdot_floor。分母0・非有限値は未定義')

register('v_max_sink va_max_sink_vtk buffer_vr_max buffer_vr_min buffer_outward_max', VELOCITY_CGS, 'cm/s', 'sink内最大速度・Alfven速度・バッファー内動径速度（外向き正）。保存値をそのまま使用')
register('buffer_cells buffer_outward_cells', 1., 'count', 'バッファーセル数・許容値を超える外向き速度のセル数')
register('smr_level_min smr_level_max smr_level_mean', 1., '1', 'root格子を0としたSMRレベル。球殻内最小・最大・体積加重平均')

DIAGNOSTIC_CRITERIA = [
    ('sink速度', '保存VTKのsink内最大|v|を表示。0は0のまま描画し、微小値への置換や丸めは行わない'),
    ('バッファー', f'外向きvrの最大値と内向きvr、外向きセル数を表示。vr>{BUFFER_VR_TOL_CGS:g} cm/s: WARN。数値面fluxは未記録でUNKNOWN'),
    ('SMRレベル', '各軸log2(rootセル幅/実セル幅)が同一整数か検査。球殻内min/max/体積加重平均を表示。root幅は実行時入力に合わせる'),
    ('初期密度プロファイル', '最初の保存出力を使用。t=0以外はWARN。球殻の体積加重平均、空の殻はNaN。合否閾値なし'),
    ('瞬間Mdot比', 'floor/fluxとreset/floorを同じ図に表示。分母0・非有限値はNaN、合否閾値なし'),
    ('実行終了', 'exit_codeが全て0: OK、非0あり: FAIL、記録なし: UNKNOWN'),
    ('終了時刻', f'最終記録時刻 >= {(EXPECTED_TLIM or 0)*TIME_KYR:.8g} kyr × (1−1e-8): OK、未到達: WARN' if EXPECTED_TLIM is not None else '目標未指定のため検査しない'),
    ('ログのNaN/Inf・dt', 'NaN/Infトークン、非有限時刻/dt、dt<=0をイベントに記録。イベントあり: WARN（要原文確認）'),
    ('小さい/急落dt', f'dt < {SMALL_DT*TIME_KYR:.8g} kyr または前記録の{DT_DROP_RATIO:g}倍未満: WARN。終端調整も含む'),
    ('時刻重複・巻戻り', 'ログ: イベント記録、HST: 区間分割、VTK: WARN。restart区間をまたいで積分しない'),
    ('HSTの値・形式', '非有限値、mass<=0、dt<=0、ヘッダ/列不整合: FAIL。ファイル/有効行なし: UNKNOWN'),
    ('Mstar更新', f'max|ΔMstar−ΔMflux| <= {BALANCE_RTOL:g} × max(|始点Mstar|, max|ΔMflux|, {1e-30*MASS_MSUN:.3g} Msun): OK、それ以外: WARN'),
    ('全域・sink質量収支', 'UNKNOWN。外部境界数値flux・EOS floor・符号付きsink fluxが不足。参考残差に合否閾値を設定しない'),
    ('VTK形式・ブロック', f'読取エラー/重複ID/ID集合変化/時刻cycle不一致/体積変化: FAIL。期待ブロック数={EXPECTED_BLOCKS}（未指定なら完全性UNKNOWN）'),
    ('VTK体積一致', f'初回とnumpy.isclose(rtol=1e-6, atol={1e-8*L_UNIT_CM**3:.8g} cm^3)。不一致時は積分値を無効化'),
    ('VTK番号', '連番の欠落: WARN（意図的間引きを含む）'),
    ('全セル有限値・密度', '保存された全cell_dataのNaN/Inf、rho<=0: FAIL。該当なし: OK（保存時刻のみ）'),
    ('極値・速度・磁場・beta', '値と位置の参考診断。自動合否閾値なし。急変・sink/格子境界との対応を確認'),
    ('floorセル', f'rho <= {DFLOOR*DENSITY_CGS:.8g} g/cm^3 × (1+1e-5)を計数。補正回数・追加質量とは異なる'),
    ('sink目標密度', f'sink内でrho ≈ {SINK_RHO*DENSITY_CGS:.8g} g/cm^3（rtol=1e-5、atol=0）を計数。合否閾値なし'),
    ('人工追加質量', 'floor_to_accreted等を表示。降着増分<=0なら比はNaN（未定義）。許容比の自動判定なし'),
    ('磁場発散', 'CTの離散発散はUNKNOWN。セル中心差分のブロック端1層を除いた参考値のみ。合否閾値なし'),
    ('dt推定', 'MHD波CFLと一定Ohmic拡散の参考値。実測dtの制限要因を断定しない。eta=0ではOhmic推定NaN'),
    ('エネルギー', '等温・重力・拡散・sinkがあるため一定を合否基準にしない'),
    ('解像度・自己重力精度', 'UNKNOWN。このコードはJeans長・ソルバ収束・解像度収束を検査しない'),
]

def display_table(title, headers, rows):
    """JupyterではHTML表、通常Pythonではテキスト表として表示。"""
    try:
        from IPython import get_ipython
        from IPython.display import HTML, display
        from html import escape
        if get_ipython() is not None:
            head = ''.join(f'<th>{escape(str(x))}</th>' for x in headers)
            body = ''.join('<tr>'+''.join(f'<td>{escape(str(x))}</td>' for x in row)+'</tr>' for row in rows)
            display(HTML(f'<h3>{escape(title)}</h3><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'))
            return
    except ImportError:
        pass
    print('\n'+title)
    print(' | '.join(headers))
    for row in rows:
        print(' | '.join(map(str,row)))

quantity_rows = [(key, unit, desc) for key, (_,unit,desc) in QUANTITIES.items()]
display_table('物理量一覧（時間: kyr / 質量: 太陽質量 / その他: CGS）', ['変数名','出力単位','意味'], quantity_rows)
display_table('診断基準一覧', ['診断項目','判定基準・制限'], DIAGNOSTIC_CRITERIA)
print('MdotはCGSのg/s、dtもkyr。無次元量・個数・IDは換算しません。')
print('入力の設定値は従来どおりシミュレーション単位です。単位基準は実行時のToyouchi.cppと照合してください。')
print('OK=検査範囲内で問題なし / WARN=要確認 / FAIL=異常 / UNKNOWN=判定不可 / INFO=説明')


def quantity(name):
    if name not in QUANTITIES:
        raise KeyError(f'物理量 {name!r} の単位が未定義です。QUANTITIESへ追加してください。')
    return QUANTITIES[name]


def physical_rows(rows):
    """内部辞書を変更せずCSV用に換算。metric/value型の極値表にも対応。"""
    converted = []
    for row in rows:
        result = {}
        for key,value in row.items():
            metric = row.get('metric') if key == 'value' else key
            if metric in QUANTITIES:
                factor,unit,_ = quantity(metric)
                result[f'{key} [{unit}]' if key != 'value' else key] = value*factor
                if key == 'value':
                    result['value_unit'] = unit
            else:
                result[key] = value
        converted.append(result)
    return converted


# %% 共通関数・VTK読取（このセルでは解析しない）
NUM = r"[+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?|inf(?:inity)?|nan)"

def number(s):
    return float(s.replace('D','E').replace('d','e'))

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


def write_csv(path, rows):
    if not rows:
        return
    rows = physical_rows(rows)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def note(report, status, item, detail):
    report.append(dict(status=status, item=item, detail=detail))


def save_plot(out, name, series, ylabel=None, log=False):
    # 異なる次元は同じ縦軸に重ねず、単位ごとにパネルを分ける。
    grouped = defaultdict(list)
    for label, t, y in series:
        metric = 'dt' if name == 'dt.png' else label.split(' [')[0]
        factor, unit, _ = quantity(metric)
        grouped[unit].append((label, np.asarray(t)*TIME_KYR, np.asarray(y)*factor))
    fig, axes = plt.subplots(max(1,len(grouped)), 1, figsize=(9,4*max(1,len(grouped))), squeeze=False)
    for ax, (unit, curves) in zip(axes[:,0], grouped.items()):
        for label,t,y in curves:
            ax.plot(t,y,label=label,lw=1)
        if log and unit != 'cm/s':
            finite = np.concatenate([y[np.isfinite(y) & (y != 0)] for _,_,y in curves])
            if finite.size:
                ax.set_yscale('symlog',linthresh=max(float(np.min(np.abs(finite)))*.1, np.finfo(float).tiny))
        if unit == 'cm/s':
            ax.axhline(0,color='0.5',lw=.6)  # 正確な0をそのまま表示
            ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useOffset=False)
        ax.set(xlabel='Time [kyr]',ylabel=unit)
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    if not grouped:
        axes[0,0].set(xlabel='Time [kyr]',title='No available data')
    fig.tight_layout()
    fig.savefig(out/name,dpi=160)
    if SHOW_PLOTS:
        plt.show()
    plt.close(fig)


# %% ログ検査の関数
def check_logs(root, out, report):
    rows, events = [], []
    pattern = re.compile(r'cycle\s*=\s*(\d+)\s+time\s*=\s*('+NUM+r')\s+dt\s*=\s*('+NUM+r')', re.I)
    paths = sorted(root.glob(LOG_GLOB))
    for path in paths:
        previous, segment, exits, final_time = None, 0, [], None
        with path.open(errors='replace') as f:
            for lineno, line in enumerate(f, 1):
                exit_match = re.search(r'exit_code\s*=\s*(-?\d+)', line)
                if exit_match:
                    exits.append(int(exit_match[1]))
                if re.search(r'(?<![\w])[-+]?(?:nan|inf(?:inity)?)(?![\w])', line, re.I):
                    events.append(dict(file=str(path), line=lineno, event='NaN/Inf token', text='非有限トークンあり。原ログの該当行を参照'))
                m = pattern.search(line)
                if not m:
                    continue
                cycle, t, dt = int(m[1]), number(m[2]), number(m[3])
                final_time = t
                if previous and (cycle <= previous['cycle'] or t <= previous['time']):
                    segment += 1
                    events.append(dict(file=str(path), line=lineno, event='time/cycle duplicate or rollback'))
                    previous = None
                row = dict(file=str(path), line=lineno, segment=segment, cycle=cycle, time=t, dt=dt)
                rows.append(row)
                if not np.isfinite([t, dt]).all() or dt <= 0:
                    events.append(dict(**row, event='invalid time/dt'))
                elif dt < SMALL_DT or (previous and dt < previous['dt'] * DT_DROP_RATIO):
                    events.append(dict(**row, event='small/drop dt (including possible terminal adjustment)'))
                previous = row
        note(report, 'OK' if exits and all(x == 0 for x in exits) else ('FAIL' if any(exits) else 'UNKNOWN'),
             '実行終了 '+path.name, f'exit codes={exits}; last recorded time={final_time*TIME_KYR if final_time is not None else None} kyr')
        if EXPECTED_TLIM is not None:
            reached = final_time is not None and np.isfinite(final_time) and final_time >= EXPECTED_TLIM * (1-1e-8)
            note(report, 'OK' if reached else 'WARN', '終了時刻 '+path.name,
                 f'目標={EXPECTED_TLIM*TIME_KYR:.8g} kyr; 最終cycle記録={final_time*TIME_KYR if final_time is not None else None} kyr。restartの途中区間は未到達でも正常。')
    write_csv(out/'log_steps.csv', rows)
    write_csv(out/'log_events.csv', events)
    groups = defaultdict(list)
    for row in rows:
        groups[(row['file'], row['segment'])].append(row)
    save_plot(out, 'dt.png', [(Path(k[0]).name+f':{k[1]}', [r['time'] for r in v], [r['dt'] for r in v]) for k,v in groups.items()], 'kyr', True)
    note(report, 'WARN' if events else ('OK' if rows else 'UNKNOWN'), 'ログ検査',
         f'{len(rows)}記録、{len(events)}参考/異常イベント。log_events.csv参照。ログなしでは正常終了を判定不可。')


# %% HST検査の関数
def check_hst(root, out, report):
    path = root/HST_FILENAME
    if not path.is_file():
        note(report, 'UNKNOWN', 'HST', f'ファイルなし: {path}')
        return
    rows, header, segment, previous = [], None, 0, None
    with path.open(errors='replace') as f:
        for lineno, line in enumerate(f, 1):
            if line.lstrip().startswith('#'):
                pairs = re.findall(r'\[(\d+)\]\s*=\s*(\S+)', line)
                if pairs:
                    header = {int(i)-1: n.replace('\\_', '_') for i,n in pairs}
                continue
            if not line.strip():
                continue
            try:
                values = [number(x) for x in line.split()]
                if not header or len(values) != max(header)+1:
                    raise ValueError('ヘッダと列数が不一致')
                unknown = set(header.values()) - QUANTITIES.keys()
                if unknown:
                    raise ValueError(f'単位未定義のHST列: {sorted(unknown)}')
                row = {name: values[i] for i,name in header.items()}
                t = row['time']
                if previous is not None and t <= previous:
                    segment += 1
                previous = t
                row.update(segment=segment, line=lineno)
                rows.append(row)
                if not np.isfinite(values).all():
                    note(report, 'FAIL', 'HST非有限値', f'line={lineno}')
                if row.get('mass', 1) <= 0 or row.get('dt', 1) <= 0:
                    note(report, 'FAIL', 'HST非正値', f'line={lineno}: mass/dt')
            except (ValueError, KeyError) as e:
                note(report, 'FAIL', 'HST読取', f'line={lineno}: {e}')
                segment += 1
                previous = None
    if not rows:
        note(report, 'UNKNOWN', 'HST記録', '有効なデータ行なし。')
    else:
        note(report, 'INFO', 'HST記録', f'{len(rows)}行を検査。瞬間Mdotの積分はせず累積量の差を使用。')
    groups = defaultdict(list)
    for row in rows:
        groups[row['segment']].append(row)
    if segment:
        note(report, 'WARN', 'HST区間', '時刻巻戻り/重複/不正行で区間を分割。区間をまたいだ差分・積分は行わない。')
    diagnostics = []
    required = {'mass','Mstar','Mflux_cum','Mreset_cum','Mfloor_cum','Msink_gas'}
    for seg, group in groups.items():
        base = group[0]
        if not all(required <= r.keys() for r in group):
            note(report, 'UNKNOWN', 'HST収支', f'区間{seg}: 必要列が不足')
            continue
        for row in group:
            d = {k: row[k]-base[k] for k in required}
            diagnostics.append(dict(time=row['time'], segment=seg,
                star_identity=d['Mstar']-d['Mflux_cum'],
                sink_proxy=d['Msink_gas']-d['Mflux_cum']+d['Mreset_cum']-d['Mfloor_cum'],
                gas_unclosed=d['mass']+d['Mreset_cum']-d['Mfloor_cum'],
                floor_added=d['Mfloor_cum'], reset_removed=d['Mreset_cum'], accreted=d['Mflux_cum'],
                floor_to_accreted=d['Mfloor_cum']/d['Mflux_cum'] if d['Mflux_cum'] > 0 else np.nan))
        scale = max(abs(base['Mstar']), max(abs(r['Mflux_cum']-base['Mflux_cum']) for r in group), 1e-30)
        residual = [r['star_identity'] for r in diagnostics if r['segment']==seg]
        good = np.isfinite(residual).all() and max(map(abs, residual)) <= BALANCE_RTOL*scale
        note(report, 'OK' if good else 'WARN', f'Mstar更新整合性 区間{seg}',
             'ΔMstar−ΔMflux_cum。定義上の一致を確認するだけで、独立した質量保存の証明ではない。')
    # 比は各HST記録から直接計算し、累積収支の必要列がなくても出力する。
    for row in rows:
        for key, numerator, denominator in [
            ('mdot_floor_to_flux','Mdot_floor','Mdot_flux'),
            ('mdot_reset_to_floor','Mdot_reset','Mdot_floor')]:
            a,b = row.get(numerator,np.nan),row.get(denominator,np.nan)
            row[key] = a/b if np.isfinite(a) and np.isfinite(b) and b!=0 else np.nan
    write_csv(out/'mdot_ratios.csv', [
        {k:r[k] for k in ('time','segment','mdot_floor_to_flux','mdot_reset_to_floor')} for r in rows])
    write_csv(out/'history.csv', rows)
    write_csv(out/'history_balance.csv', diagnostics)
    for filename, names, source in [
        ('history_mass.png', ['mass','Mstar','Msink_gas'], rows),
        ('history_sink.png', ['Bmax_sink','Mdot_floor','Mdot_reset'], rows),
        ('mass_balance_proxies.png', ['star_identity','sink_proxy','gas_unclosed'], diagnostics),
        ('floor_mass.png', ['floor_added','reset_removed','accreted'], diagnostics),
        ('mdot_ratios.png', ['mdot_floor_to_flux','mdot_reset_to_floor'], rows),
        ('floor_ratio.png', ['floor_to_accreted'], diagnostics),
        ('history_energy.png', ['1-KE','2-KE','3-KE','grav-E','1-ME','2-ME','3-ME'], rows)]:
        curves = []
        for seg in sorted({r['segment'] for r in source}):
            group = [r for r in source if r['segment']==seg]
            for name in names:
                if all(name in r for r in group):
                    curves.append((f'{name} [{seg}]', [r['time'] for r in group], [r[name] for r in group]))
        save_plot(out, filename, curves, None, filename in ('history_sink.png','mdot_ratios.png'))
    note(report, 'UNKNOWN', '厳密な全領域質量収支',
         'gas_unclosed=ΔMgas+ΔMreset−ΔMfloor。外部境界の累積数値fluxとEOS floor追加量が未記録なので保存誤差ではない。')
    note(report, 'UNKNOWN', '厳密なsink質量収支',
         'sink_proxy=ΔMsink−ΔMflux+ΔMreset−ΔMfloor。Mfluxは各cycleで負値を切り捨てており、符号付きfluxが必要。')

    return rows

# %% VTK検査・断面図の関数
def block_metrics(path):
    time, cycle, coords, fields = read_vtk(path)
    widths = [np.diff(c) for c in coords]
    if not all(np.isfinite(c).all() for c in coords) or any(np.any(w <= 0) for w in widths):
        raise ValueError('座標が非有限または非単調')
    rho = field(fields, ['rho', 'dens'])
    vel = field(fields, ['vel'])
    b = field(fields, ['Bcc', 'bcc'])
    n = len(rho)
    if rho.shape != (n,) or vel.shape != (n,3) or b.shape != (n,3):
        raise ValueError('rho/vel/Bccの形状不正')
    shape = tuple(len(w) for w in widths[::-1])
    centers = [(c[1:]+c[:-1])/2 for c in coords]
    z,y,x = np.meshgrid(centers[2], centers[1], centers[0], indexing='ij')
    xyz = np.column_stack([x.ravel(),y.ravel(),z.ravel()])
    volume = (widths[2][:,None,None]*widths[1][None,:,None]*widths[0][None,None,:]).ravel()
    sink = np.linalg.norm(xyz-SINK_CENTER, axis=1) < SINK_RADIUS
    nonfinite = {k: int(np.count_nonzero(~np.isfinite(v))) for k,v in fields.items()}
    speed, bmag = np.linalg.norm(vel,axis=1), np.linalg.norm(b,axis=1)
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        va = bmag / np.sqrt(np.where(rho>0,rho,np.nan))
        beta = 2*rho*CS**2 / np.where(bmag>0,bmag*bmag,np.nan)
    result = dict(time=time, cycle=cycle, cells=n, volume=float(volume.sum()),
                  nonfinite=sum(nonfinite.values()), nonpositive_rho=int(np.count_nonzero(rho<=0)),
                  mass=float(np.sum(rho*volume)), magnetic_energy=float(np.sum(.5*bmag**2*volume)),
                  sink_mass=float(np.sum(rho[sink]*volume[sink])),
                  floor_cells=int(np.count_nonzero(rho<=DFLOOR*(1+1e-5))),
                  floor_cells_outside=int(np.count_nonzero((rho<=DFLOOR*(1+1e-5)) & ~sink)),
                  sink_target_cells=int(np.count_nonzero(sink & np.isclose(rho,SINK_RHO,rtol=1e-5,atol=0))))
    positions = []
    for name, values, minimum in [('rho_min',rho,True),('rho_max',rho,False),('v_max',speed,False),('B_max',bmag,False),('va_max',va,False),('beta_min',beta,True),('va_max_outside',np.where(sink,np.nan,va),False),('v_max_sink',np.where(sink,speed,np.nan),False),('va_max_sink_vtk',np.where(sink,va,np.nan),False)]:
        valid = np.flatnonzero(np.isfinite(values))
        idx = valid[np.argmin(values[valid]) if minimum else np.argmax(values[valid])] if valid.size else None
        result[name] = float(values[idx]) if idx is not None else np.nan
        if idx is not None:
            positions.append(dict(metric=name,value=result[name],x=xyz[idx,0],y=xyz[idx,1],z=xyz[idx,2],file=str(path),time=time))
    # 保存時刻の波速からの目安。実際のintegrator/dt成長制限等は再現しない。
    dt_wave = np.full(n,np.inf)
    dxs = [np.broadcast_to(w.reshape((1,1,-1) if a==0 else ((1,-1,1) if a==1 else (-1,1,1))),shape).ravel() for a,w in enumerate(widths)]
    with np.errstate(invalid='ignore',divide='ignore'):
        for a in range(3):
            total = CS**2+va**2
            cf = np.sqrt(.5*(total+np.sqrt(np.maximum(total**2-4*CS**2*b[:,a]**2/rho,0))))
            dt_wave = np.minimum(dt_wave,CFL*dxs[a]/(np.abs(vel[:,a])+cf))
    radius = np.linalg.norm(xyz-SINK_CENTER,axis=1)
    dx_min = np.minimum.reduce(dxs)
    buffer = (radius>=SINK_RADIUS)&(radius<SINK_RADIUS+SINK_NO_OUTFLOW_CELLS*dx_min)&(radius>0)
    vr = np.divide(np.sum(vel*(xyz-SINK_CENTER),axis=1),radius,
                   out=np.full(n,np.nan),where=radius>0)
    valid_buffer = buffer & np.isfinite(vr)
    result['buffer_cells'] = int(buffer.sum())
    result['buffer_outward_cells'] = int(np.count_nonzero(valid_buffer & (vr*VELOCITY_CGS>BUFFER_VR_TOL_CGS)))
    result['buffer_vr_max'] = float(vr[valid_buffer].max()) if valid_buffer.any() else np.nan
    result['buffer_vr_min'] = float(vr[valid_buffer].min()) if valid_buffer.any() else np.nan
    result['buffer_outward_max'] = max(0.,result['buffer_vr_max']) if valid_buffer.any() else np.nan
    result['dt_wave_est'] = float(np.min(dt_wave))
    result['dt_ohm_est'] = float(CFL*min(w.min() for w in widths)**2/(6*ETA_OHM)) if ETA_OHM>0 else np.nan
    # セル中心Bの微分。ブロック端1層を除外し、CTの離散発散とは明確に区別。
    result.update(div_proxy_max=np.nan, div_proxy_sqvol=0., div_proxy_volume=0.)
    if min(shape)>=3:
        div = sum(np.gradient(b[:,a].reshape(shape),centers[a],axis=2-a,edge_order=2) for a in range(3))
        interior = (slice(1,-1),)*3
        values = div[interior].ravel()
        vv = volume.reshape(shape)[interior].ravel()
        good = np.isfinite(values)
        if good.any():
            result.update(div_proxy_max=float(np.max(np.abs(values[good]))),
                          div_proxy_sqvol=float(np.sum(values[good]**2*vv[good])),div_proxy_volume=float(vv[good].sum()))
    return result, positions, nonfinite, (coords, rho.reshape(shape), bmag.reshape(shape), va.reshape(shape))


def save_slice(out, step, blocks):
    # y=0を含むセルのx-z矩形。補間なし、各ブロックの実際のセル幅を使用。
    polygons, values = [], [[],[],[]]
    for coords, arrays in blocks:
        xc,yc,zc = coords
        j = np.searchsorted(yc,SINK_CENTER[1],side='right')-1
        if not 0 <= j < len(yc)-1:
            continue
        for k in range(len(zc)-1):
            for i in range(len(xc)-1):
                if SLICE_HALF_WIDTH is not None and (abs((xc[i]+xc[i+1])/2-SINK_CENTER[0])>SLICE_HALF_WIDTH or abs((zc[k]+zc[k+1])/2-SINK_CENTER[2])>SLICE_HALF_WIDTH):
                    continue
                polygons.append([(xc[i],zc[k]),(xc[i+1],zc[k]),(xc[i+1],zc[k+1]),(xc[i],zc[k+1])])
                for dest, arr in zip(values,arrays):
                    dest.append(arr[k,j,i])
    fig,axes = plt.subplots(1,3,figsize=(15,5))
    polygons = np.asarray(polygons)*L_UNIT_CM
    for ax,metric,data in zip(axes,['rho','Bcc','va_max'],values):
        factor,unit,_ = quantity(metric)
        name = f'{metric} [{unit}]'
        a = np.asarray(data)*factor
        good = np.isfinite(a)&(a>0)
        if good.any():
            lo,hi = a[good].min(),a[good].max()
            coll = PolyCollection(polygons,array=np.ma.masked_where(~good,a),norm=LogNorm(lo,max(hi,lo*(1+1e-6))),cmap='viridis',edgecolors='none')
            ax.add_collection(coll)
            ax.autoscale_view()
            fig.colorbar(coll,ax=ax,label=unit)
        ax.add_patch(plt.Circle((SINK_CENTER[0]*L_UNIT_CM,SINK_CENTER[2]*L_UNIT_CM),SINK_RADIUS*L_UNIT_CM,fill=False,color='red',lw=.7))
        ax.set(title=name,xlabel='x [cm]',ylabel='z [cm]',aspect='equal')
    fig.suptitle(f'Output {step}: y={SINK_CENTER[1]*L_UNIT_CM:.4g} cm cell slice')
    fig.tight_layout()
    fig.savefig(out/f'slice_{step:05d}.png',dpi=160)
    if SHOW_PLOTS:
        plt.show()
    plt.close(fig)


def save_initial_density_profile(paths, out, report):
    """最初の保存出力の球殻体積加重平均。セル中心でbin分類し、sink内も含む。"""
    def geometry(coords):
        centers = [(c[:-1]+c[1:])/2 for c in coords]
        z,y,x = np.meshgrid(centers[2],centers[1],centers[0],indexing='ij')
        r = np.sqrt((x-SINK_CENTER[0])**2+(y-SINK_CENTER[1])**2+(z-SINK_CENTER[2])**2).ravel()
        dx,dy,dz = [np.diff(c) for c in coords]
        v = (dz[:,None,None]*dy[None,:,None]*dx[None,None,:]).ravel()
        return r,v

    rmin,rmax = np.inf,0.
    for path in paths:
        _,_,coords,_ = read_vtk(path)
        r,_ = geometry(coords)
        positive = r[r>0]
        if positive.size:
            rmin = min(rmin,float(positive.min()))
        rmax = max(rmax,float(r.max()))
    if not np.isfinite(rmin) or rmax<=0:
        note(report,'UNKNOWN','初期密度プロファイル','正の半径を持つセルがない。')
        return
    # 中心セルも取りこぼさないよう最初の殻はr=0から。
    edges = np.r_[0.,np.geomspace(rmin*.5,rmax*(1+1e-12),PROFILE_BINS)]
    mass = np.zeros(PROFILE_BINS)
    volumes = np.zeros(PROFILE_BINS)
    counts = np.zeros(PROFILE_BINS,dtype=int)
    level_min = np.full(PROFILE_BINS,np.inf)
    level_max = np.full(PROFILE_BINS,-np.inf)
    level_sum = np.zeros(PROFILE_BINS)
    levels_valid = True
    for path in paths:
        time,_,coords,fields = read_vtk(path)
        r,v = geometry(coords)
        rho = field(fields,['rho','dens'])
        valid = np.isfinite(rho)&(rho>0)
        if not valid.all():
            note(report,'FAIL','初期密度プロファイル','非有限/非正密度があるためプロファイル作成を中止。')
            return
        widths = [np.diff(c) for c in coords]
        levels = [np.log2(ROOT_CELL_WIDTH[a]/w) for a,w in enumerate(widths)]
        level = int(round(float(levels[0][0])))
        if level<0 or not all(np.allclose(v,level,rtol=0,atol=1e-4) for v in levels):
            levels_valid = False
        shell = np.searchsorted(edges,r,side='right')-1
        ok = (shell>=0)&(shell<PROFILE_BINS)
        np.minimum.at(level_min,shell[ok],level)
        np.maximum.at(level_max,shell[ok],level)
        level_sum += np.histogram(r,bins=edges,weights=v*level)[0]
        mass += np.histogram(r,bins=edges,weights=rho*v)[0]
        volumes += np.histogram(r,bins=edges,weights=v)[0]
        counts += np.histogram(r,bins=edges)[0]
    mean = np.divide(mass,volumes,out=np.full(PROFILE_BINS,np.nan),where=volumes>0)
    radius = np.sqrt(edges[:-1]*edges[1:])
    radius[0] = edges[1]/2
    rows = [dict(time=time,radius=radius[i],radius_inner=edges[i],radius_outer=edges[i+1],
                 rho_mean=mean[i],volume=volumes[i],cells=int(counts[i])) for i in range(PROFILE_BINS)]
    write_csv(out/'initial_density_profile.csv',rows)
    if levels_valid:
        level_mean = np.divide(level_sum,volumes,out=np.full(PROFILE_BINS,np.nan),where=volumes>0)
        level_min[volumes==0] = np.nan
        level_max[volumes==0] = np.nan
        write_csv(out/'initial_smr_profile.csv', [dict(time=time,radius=radius[i],
            smr_level_min=level_min[i],smr_level_max=level_max[i],smr_level_mean=level_mean[i]) for i in range(PROFILE_BINS)])
        fig,ax = plt.subplots(figsize=(8,5))
        ax.fill_between(radius*L_UNIT_CM,level_min,level_max,alpha=.2,label='Shell min-max')
        ax.semilogx(radius*L_UNIT_CM,level_mean,marker='.',label='Volume-weighted mean')
        ax.axvline(SINK_RADIUS*L_UNIT_CM,color='red',ls='--',label='Sink radius')
        ax.set(xlabel='Radius [cm]',ylabel='SMR level (root = 0)',title=f'First saved mesh: t={time*TIME_KYR:.6g} kyr')
        ax.grid(alpha=.25); ax.legend(); fig.tight_layout()
        fig.savefig(out/'initial_smr_profile.png',dpi=160)
        if SHOW_PLOTS:
            plt.show()
        plt.close(fig)
    else:
        note(report,'UNKNOWN','SMRレベル','設定rootセル幅から同一整数レベルを復元できません。ROOT_CELL_WIDTHを確認してください。')

    fig,ax = plt.subplots(figsize=(8,5))
    ax.loglog(radius*L_UNIT_CM,mean*DENSITY_CGS,marker='.',label='Volume-weighted shell mean')
    ax.axvline(SINK_RADIUS*L_UNIT_CM,color='red',ls='--',label='Sink radius')
    ax.set(xlabel='Radius [cm]',ylabel='Density [g/cm^3]',title=f'First saved density profile: t={time*TIME_KYR:.6g} kyr')
    ax.grid(alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(out/'initial_density_profile.png',dpi=160)
    if SHOW_PLOTS:
        plt.show()
    plt.close(fig)
    note(report,'INFO' if time==0 else 'WARN','初期密度プロファイル',
         f'最初の出力 t={time*TIME_KYR:.8g} kyr。sink内を含む球殻体積加重平均。t=0でなければ初期条件そのものではない。外側の殻は計算領域内の部分のみ。')


def check_vtk(root, out, report, hst_data=None):
    groups = defaultdict(list)
    for path in sorted(root.rglob(VTK_GLOB)):
        if out in path.parents:
            continue
        match = re.search(r'\.out\d+\.(\d+)\.vtk$',path.name)
        if match:
            groups[int(match[1])].append(path)
    steps = sorted(groups)
    if not steps:
        note(report,'UNKNOWN','VTK','対応するVTKなし。全セル検査・断面図は未実施。')
        return
    selected = set(SLICE_STEPS if SLICE_STEPS is not None else [steps[0],steps[len(steps)//2],steps[-1]])
    rows, positions, bad_fields = [], [], []
    baseline_ids, baseline_volume, previous_time = None, None, None
    sum_keys = ['cells','volume','nonfinite','nonpositive_rho','mass','magnetic_energy','sink_mass','floor_cells','floor_cells_outside','sink_target_cells','div_proxy_sqvol','div_proxy_volume','buffer_cells','buffer_outward_cells']
    min_keys = ['buffer_vr_min','rho_min','beta_min','dt_wave_est','dt_ohm_est']
    max_keys = ['v_max_sink','va_max_sink_vtk','buffer_vr_max','buffer_outward_max','rho_max','v_max','B_max','va_max','va_max_outside','div_proxy_max']
    for step in steps:
        print(f'[VTK] output={step:05d}, blocks={len(groups[step])}',flush=True)
        blocks, slices, ids = [], [], []
        failed = False
        for path in groups[step]:
            match = re.search(r'\.block(\d+)\.',path.name)
            ids.append(int(match[1]) if match else path.name)
            try:
                metrics, extrema, invalid, data = block_metrics(path)
                blocks.append(metrics)
                positions.extend(dict(output=step,**p) for p in extrema)
                bad_fields.extend(dict(output=step,file=str(path),field=k,count=v) for k,v in invalid.items() if v)
                if step in selected:
                    coords, *arrays = data
                    j = np.searchsorted(coords[1],SINK_CENTER[1],side='right')-1
                    if 0 <= j < len(coords[1])-1:
                        # 切断面の1層のみ保持。全領域をメモリに載せない。
                        sliced_coords = [coords[0], coords[1][j:j+2], coords[2]]
                        slices.append((sliced_coords,[a[:,j:j+1,:].copy() for a in arrays]))
            except (ValueError,OSError,UnicodeError,IndexError) as e:
                failed = True
                note(report,'FAIL',f'VTK読取 output={step}',f'{path}: {e}')
        if not blocks:
            continue
        current_ids = set(ids)
        complete = not failed and len(ids)==len(current_ids)
        if baseline_ids is None:
            baseline_ids = current_ids
        elif current_ids != baseline_ids:
            complete = False
        if EXPECTED_BLOCKS is not None and len(ids)!=EXPECTED_BLOCKS:
            complete = False
        time,cycle = blocks[0]['time'],blocks[0]['cycle']
        if not np.isfinite(time) or any(b['time']!=time or b['cycle']!=cycle for b in blocks):
            complete = False
        if previous_time is not None and time<=previous_time:
            note(report,'WARN','VTK時刻',f'output={step}: 番号順で時刻の重複/巻戻り。')
        previous_time = time
        row = dict(output=step,time=time,cycle=cycle,blocks=len(blocks),consistent=complete)
        for key in sum_keys:
            row[key] = sum(b[key] for b in blocks)
        for key in min_keys+max_keys:
            vals = np.array([b[key] for b in blocks])
            finite = vals[np.isfinite(vals)]
            row[key] = float((np.min if key in min_keys else np.max)(finite)) if finite.size else np.nan
        if baseline_volume is None:
            baseline_volume = row['volume']
        if not np.isclose(row['volume'],baseline_volume,rtol=1e-6):
            complete = False
            row['consistent'] = False
        if not complete:
            note(report,'FAIL','VTKブロック整合性',f'output={step}: 欠落・重複・時刻不一致・体積変化の可能性。積分値を無効化。')
            for key in ['mass','magnetic_energy','sink_mass']:
                row[key] = np.nan
        row['div_proxy_rms'] = np.sqrt(row['div_proxy_sqvol']/row['div_proxy_volume']) if row['div_proxy_volume']>0 else np.nan
        row['floor_fraction'] = row['floor_cells']/row['cells']
        rows.append(row)
        if step == steps[0] and complete:
            save_initial_density_profile(groups[step],out,report)
        elif step == steps[0]:
            note(report,'UNKNOWN','初期密度プロファイル','最初の出力が不完全なため作成しない。')
        if slices:
            save_slice(out,step,slices)
    write_csv(out/'vtk_history.csv',rows)
    write_csv(out/'extrema_locations.csv',positions)
    write_csv(out/'nonfinite_fields.csv',bad_fields)
    missing = [(a+1,b-1) for a,b in zip(steps,steps[1:]) if b>a+1]
    note(report,'WARN' if missing else 'OK','VTK番号連続性',f'欠落番号範囲={missing}。意図的な間引きも含む。')
    note(report,'UNKNOWN' if EXPECTED_BLOCKS is None else 'INFO','VTK完全性',
         f'EXPECTED_BLOCKS={EXPECTED_BLOCKS}。ブロックID集合・体積の時間的一致は検査するが、全時刻共通の欠落やghost/空間重複は保証しない。')
    note(report,'FAIL' if bad_fields or any(r['nonpositive_rho'] for r in rows) else ('OK' if rows else 'UNKNOWN'),
         '全セル有限値・密度正値','全cell_dataを走査。FAILした出力の派生極値は有限セルだけの参考値。')
    for filename,keys in [('extrema.png',['rho_min','rho_max','v_max','B_max','va_max','va_max_outside','v_max_sink','va_max_sink_vtk']),
                          ('floor_cells.png',['floor_cells','floor_cells_outside','sink_target_cells','blocks']),
                          ('sink_buffer.png',['buffer_vr_min','buffer_vr_max','buffer_outward_max','buffer_cells','buffer_outward_cells']),
                          ('vtk_mass_energy.png',['mass','sink_mass','magnetic_energy']),
                          ('dt_estimates.png',['dt_wave_est','dt_ohm_est']),
                          ('divB_proxy.png',['div_proxy_max','div_proxy_rms'])]:
        curves = [(k,[r['time'] for r in rows],[r[k] for r in rows]) for k in keys]
        if filename == 'extrema.png':
            for segment in sorted({r['segment'] for r in (hst_data or [])}):
                group = [r for r in hst_data if r['segment']==segment and 'Va_max_sink' in r]
                if group:
                    curves.append((f'Va_max_sink [HST {segment}]',[r['time'] for r in group],[r['Va_max_sink'] for r in group]))
        save_plot(out,filename,curves,None,True)
    note(report,'UNKNOWN' if not any(r['buffer_cells'] for r in rows) else ('WARN' if any(r['buffer_outward_cells'] for r in rows) else 'INFO'),'バッファー速度',
         f'外向き正、許容値={BUFFER_VR_TOL_CGS:g} cm/s。対象セル数0なら未評価。初期出力は処理前の可能性あり。sink_buffer.png参照。')
    note(report,'UNKNOWN','バッファー数値質量flux','VTKは処理後セル中心速度のみ。外向き面fluxを厳密に禁止した証明にはならない。')
    note(report,'UNKNOWN','CTの離散発散',
         '面中心磁場が未保存。divB_proxyはセル中心Bの微分（各ブロック端1層を除外）。CT誤差や細分化境界の評価には使用不可。')
    note(report,'INFO','保存状態の時間刻み推定',
         'dt_wave_est: fast波CFL目安。dt_ohm_est: CFL*dx_min^2/(6*eta)。実測dtの制限要因を断定しない。')
    note(report,'INFO','floor適用数',
         'floor_cellsは保存密度がhydro/dfloor以下のセル数。補正回数や追加質量ではない。sink_target_cellsも目標密度との一致数。')


# %% 出力先の準備（設定変更後はこのセルから下を再実行）
root = VTK_DIR.expanduser().resolve()
if not root.is_dir():
    raise FileNotFoundError(f'設定セルのVTK_DIRを変更してください: {root}')
if np.shape(ROOT_CELL_WIDTH)!=(3,) or not np.isfinite(ROOT_CELL_WIDTH).all() or np.any(ROOT_CELL_WIDTH<=0):
    raise ValueError('ROOT_CELL_WIDTHには各軸の正のrootセル幅を3つ指定してください')
if not np.isfinite([SINK_NO_OUTFLOW_CELLS,BUFFER_VR_TOL_CGS]).all() or SINK_NO_OUTFLOW_CELLS<0 or BUFFER_VR_TOL_CGS<0:
    raise ValueError('バッファー幅・速度許容値は非負の有限値にしてください')
if not isinstance(PROFILE_BINS,int) or PROFILE_BINS < 2:
    raise ValueError('PROFILE_BINSは2以上の整数にしてください')
if CS <= 0 or CFL <= 0 or ETA_OHM < 0 or SINK_RADIUS <= 0 or DFLOOR <= 0:
    raise ValueError('音速・CFL・sink半径・floorは正、etaは非負に設定してください')
out = root/'checks'
out.mkdir(exist_ok=True)
# 設定ファイルにも物理単位と換算基準を明示する。
settings = {
    'VTK_DIR': str(root), 'HST_FILENAME': HST_FILENAME,
    'LOG_GLOB': LOG_GLOB, 'VTK_GLOB': VTK_GLOB, 'EXPECTED_BLOCKS': EXPECTED_BLOCKS,
    'EXPECTED_TLIM [kyr]': EXPECTED_TLIM*TIME_KYR if EXPECTED_TLIM is not None else None,
    'DFLOOR [g/cm^3]': DFLOOR*DENSITY_CGS, 'SINK_RHO [g/cm^3]': SINK_RHO*DENSITY_CGS,
    'SINK_RADIUS [cm]': SINK_RADIUS*L_UNIT_CM, 'SINK_CENTER [cm]': (SINK_CENTER*L_UNIT_CM).tolist(),
    'CS [cm/s]': CS*VELOCITY_CGS, 'CFL': CFL, 'ETA_OHM [cm^2/s]': ETA_OHM*L_UNIT_CM**2/T_UNIT_S,
    'SMALL_DT [kyr]': SMALL_DT*TIME_KYR, 'DT_DROP_RATIO': DT_DROP_RATIO, 'BALANCE_RTOL': BALANCE_RTOL,
    'SLICE_STEPS': SLICE_STEPS, 'PROFILE_BINS': PROFILE_BINS,
    'SLICE_HALF_WIDTH [cm]': SLICE_HALF_WIDTH*L_UNIT_CM if SLICE_HALF_WIDTH is not None else None,
    'SHOW_PLOTS': SHOW_PLOTS,
    'ROOT_CELL_WIDTH [cm]': (ROOT_CELL_WIDTH*L_UNIT_CM).tolist(),
    'SINK_NO_OUTFLOW_CELLS': SINK_NO_OUTFLOW_CELLS, 'BUFFER_VR_TOL_CGS': BUFFER_VR_TOL_CGS,
    'input_unit_scales': {'mass [g]': M_UNIT_G, 'length [cm]': L_UNIT_CM, 'time [s]': T_UNIT_S,
                          'magnetic_field [G]': float(B_GAUSS), 'solar_mass [g]': MSUN_G},
}
(out/'settings.json').write_text(json.dumps(settings,ensure_ascii=False,indent=2),encoding='utf-8')

write_csv(out/'quantity_definitions.csv', [dict(variable=k, unit=u, description=d) for k,u,d in quantity_rows])
write_csv(out/'diagnostic_criteria.csv', [dict(item=k, criterion=v) for k,v in DIAGNOSTIC_CRITERIA])

# 各検査セルの再実行で判定が重複しないよう、結果を分けて保持する。
log_report, hst_report, vtk_report = None, None, None
hst_data = None
print(f'保存先: {out}')

# %% 1. ログ・時間刻みの確認
log_report = []
for name in ['log_events.csv', 'log_steps.csv']:
    (out/name).write_text('', encoding='utf-8')
check_logs(root, out, log_report)

# %% 2. HST・質量収支・floorの確認
hst_report = []
for name in ['history.csv', 'history_balance.csv', 'mdot_ratios.csv']:
    (out/name).write_text('', encoding='utf-8')
hst_data = check_hst(root, out, hst_report)

# %% 3. VTK全セル・極値・磁場・断面図の確認（時間がかかるセル）
vtk_report = []
for name in ['vtk_history.csv', 'extrema_locations.csv', 'nonfinite_fields.csv', 'initial_density_profile.csv', 'initial_smr_profile.csv']:
    (out/name).write_text('', encoding='utf-8')
check_vtk(root, out, vtk_report, hst_data)

# %% 4. 判定結果の表示・保存（検査セルを再実行した後はこのセルも実行）
report = []
for label, results in [('ログ', log_report), ('HST', hst_report), ('VTK', vtk_report)]:
    if results is None:
        note(report, 'UNKNOWN', label, 'この検査セルは未実行。')
    else:
        report.extend(results)
note(report,'INFO','エネルギー', '等温・重力・Ohmic拡散・sink処理のため、総エネルギー一定を合否基準にしない。')
note(report,'UNKNOWN','解像度収束・自己重力精度', 'Jeans長の解像度・重力ソルバ収束・解像度変更比較はこの検査には含まない。')
write_csv(out/'check_report.csv',report)
lines = ['Toyouchi 健全性確認',f'入力: {root}',
         'OK=その検査範囲で問題なし / WARN=要確認 / FAIL=異常 / UNKNOWN=判定不可 / INFO=説明',
         '入力は3D Cartesian・ghostなし・静的細分化を想定。実行時設定・commitを照合すること。',
         'CSVと図は時間kyr・質量Msun・その他CGS（Mdotはg/s）。欠落項目を正常とは判定しない。総合的な物理妥当性を保証しない。',
         'この実行で生成しなかった過去のPNGが残る場合がある。今回のCSV・reportを基準にする。','']
lines += [f"[{r['status']}] {r['item']}: {r['detail']}" for r in report]
(out/'report.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('\n'.join(lines))
print(f'保存先: {out}')
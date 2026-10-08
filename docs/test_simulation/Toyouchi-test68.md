# Toyouchi-test68

## 1. 目的
Bz=10μGとして比較計算

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 57034e1
変更点: Bz=10μG

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlld
  #Magnetic fields:              ON
  #Number of scalars:            0
  #Number of chemical species:   0
  #Special relativity:           OFF
  #General relativity:           OFF
  #Radiative Transfer:           OFF
  #Implicit Radiation:           OFF
  #Cosmic Ray Transport:         OFF
  #Cosmic Ray Diffusion:         OFF
  #Frame transformations:        OFF
  #Self-Gravity:                 Multigrid
  #Super-Time-Stepping:          OFF
  #Chemistry:                    OFF
  #KIDA rates:                   OFF
  #ChemRadiation:                OFF
  #chem_ode_solver:              OFF
  #Debug flags:                  OFF
  #Code coverage flags:          OFF
  #Linker flags:                  
  #Floating-point precision:     double
  #Number of ghost cells:        4
  #MPI parallelism:              ON
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test68.sh

## 6.計算ログ
cycle=183114 time=4.3859861066499775e+01 dt=1.3893350022442519e-04
Sink: Mstar=2.3549914521597711e+02 Mdot_flux=4.5695011923169329e+00 dM_flux=6.3485679492827581e-04 Mdot_reset=4.9198385512999527e+00 Mdot_floor=3.5033735898299506e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=1.2774024259574424e-01 Bmax_sink=3.2062546427887383e-02
cycle=183115 time=4.3859999999999999e+01 dt=2.0840115086638412e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=183115
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 50440125440
cpu time used  = 9.2768708071000001e+04
zone-cycles/cpu_second = 5.4371917523520906e+05
end_time   = 2026-10-01T17:55:40+09:00
exit_code  = 0

## 7.結果・考察等
現在までのテスト計算では、磁場強度が高いほど中心星質量進化が抑えられるようである。
詳細は解析してみないと分からないが、おそらく50kyr程度では磁気張力よりも磁気圧が卓越し、降着を抑えている可能性がある。
初期磁場の向きを、円盤角運動量ベクトルに平行ではなく、角度をつけてみると早い段階から磁気張力が効いて降着率が上がるかもしれないので、試してみる価値はありそう。
2✕10^4AUの範囲で解像度を調べたところ、特にリングやスパイラルなどの構造体を十分な解像度（最低条件Nj>4cell）で観測できていないようである。
Novaでのテスト計算時は計算時間節約のためSMRを用い、CFCAでの本計算では上記の範囲に限定してAMRを用いて計算したほうが良さそう。

## 8.備考
VTKファイルはすべて削除

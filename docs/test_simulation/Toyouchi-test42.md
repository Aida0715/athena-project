# Toyouchi-test42

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: 非回転、sink_va_cap = 0.0、sink_cfl_va_cap = 10.0、sink_cfl_cap_radius_factor = 1.0、sink_cfl_max_relax = 1.2

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlle
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
  #Self-Gravity:                 multigrid
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
Toyouchi-test42.sh

## 6.計算ログ
cycle=1261 time=1.9999529622053316e+01 dt=4.7037794668369770e-04
Sink: Mstar=1.5312991133041516e+01 Mdot_flux=1.7529044144993049e+00 dM_flux=8.2452757922497238e-04 Mdot_reset=1.7529044191253425e+00 Mdot_floor=4.6260376001253385e-09 Msink_gas=4.6224609375000084e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=2.2114079247440806e+01 Bmax_sink=2.2114079247440807e-03
cycle=1262 time=2.0000000000000000e+01 dt=9.4075589336739540e-04
Terminating on time limit
time=2.0000000000000000e+01 cycle=1262
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 288; 224  created, 0 destroyed during this simulation.
zone-cycles = 1488715776
cpu time used  = 1.0594159433000001e+04
zone-cycles/cpu_second = 1.4052231188467520e+05
end_time   = 2026-08-26T02:25:26+09:00
exit_code  = 0

## 7.結果・考察等
Aの結果と一致。CFLタイムステップ緩和1.2倍は問題なし。

## 8.備考
VTKファイルはすべて削除

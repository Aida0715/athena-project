# Toyouchi-test43

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: 非回転、sink_va_cap = 0.0、sink_cfl_va_cap = 10.0、sink_cfl_cap_radius_factor = 1.0、sink_cfl_max_relax = 1.5

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
Toyouchi-test43.sh

## 6.計算ログ
cycle=1152 time=1.9995531416608948e+01 dt=4.4685833910520500e-03
Sink: Mstar=1.5312320330922523e+01 Mdot_flux=1.7415396032258545e+00 dM_flux=7.7822149458344307e-03 Mdot_reset=1.7415396426542857e+00 Mdot_floor=3.9428431021401991e-08 Msink_gas=4.6224609375000084e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=2.2542544311754263e+01 Bmax_sink=2.2542544311754263e-03
cycle=1153 time=2.0000000000000000e+01 dt=6.5785893464324358e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=1153
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 288; 224  created, 0 destroyed during this simulation.
zone-cycles = 1360134144
cpu time used  = 1.0071785806000000e+04
zone-cycles/cpu_second = 1.3504399023157702e+05
end_time   = 2026-08-26T02:12:56+09:00
exit_code  = 0

## 7.結果・考察等
Aと一致。CFLタイムステップ1.5倍も問題なし。

## 8.備考
VTKファイルはすべて削除

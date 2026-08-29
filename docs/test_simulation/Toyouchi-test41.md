# Toyouchi-test41

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: 非回転、sink_va_cap = 10.0、sink_cfl_va_cap = 0.0、sink_cfl_cap_radius_factor = 1.0、sink_cfl_max_relax = 1.0

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
Toyouchi-test41.sh

## 6.計算ログ
cycle=1271 time=1.9996563165687164e+01 dt=3.4368343128363676e-03
Sink: Mstar=1.5312263784784530e+01 Mdot_flux=1.7462061167705227e+00 dM_flux=6.0014210994016812e-03 Mdot_reset=1.7462066326735459e+00 Mdot_floor=4.2488626848842043e-07 Msink_gas=4.9787012728549508e-07 Msink_magfloor=3.5624033535494450e-08 Va_max_sink=1.0000000000000002e+01 Bmax_sink=1.7618911515939822e-03
cycle=1272 time=2.0000000000000000e+01 dt=6.8736686256727353e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=1272
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 288; 224  created, 0 destroyed during this simulation.
zone-cycles = 1500512256
cpu time used  = 1.0704615471999999e+04
zone-cycles/cpu_second = 1.4017432573125875e+05
end_time   = 2026-08-26T02:24:56+09:00
exit_code  = 0

## 7.結果・考察等
他の非回転テストと降着率は同じだが、シンク内磁場強度などの物理量が異なる。
人工的に質量増加させる方法がMHD状態を変更している可能性がある。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test44

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: 回転あり、sink_va_cap = 0.0、sink_cfl_va_cap = 0.0、sink_cfl_cap_radius_factor = 1.0、sink_cfl_max_relax = 1.0

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
Toyouchi-test44.sh

## 6.計算ログ
cycle=3690 time=1.9999933884077414e+01 dt=6.6115922585652243e-05
Sink: Mstar=4.4254023802999754e+01 Mdot_flux=2.3039028954681937e+00 dM_flux=1.5232466548163514e-04 Mdot_reset=2.3039763895795269e+00 Mdot_floor=7.3494111332922534e-05 Msink_gas=4.6224609375000084e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=8.1938645434697222e+01 Bmax_sink=8.1938645434697227e-03
cycle=3691 time=2.0000000000000000e+01 dt=1.3223184517130449e-04
Terminating on time limit
time=2.0000000000000000e+01 cycle=3691
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 512; 448  created, 0 destroyed during this simulation.
zone-cycles = 5958107136
cpu time used  = 1.5104321236000000e+04
zone-cycles/cpu_second = 3.9446374603046081e+05
end_time   = 2026-08-26T23:52:14+09:00
exit_code  = 0

## 7.結果・考察等
回転基準として正常。

## 8.備考
VTKファイルはすべて削除

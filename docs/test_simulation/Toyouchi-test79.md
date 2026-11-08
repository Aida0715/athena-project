# Toyouchi-test79

## 1. 目的
磁場を5μGとし、比較計算用

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a4554e8
変更点: Bz0=5μG

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
Toyouchi-test79.sh

## 6.計算ログ
cycle=322802 time=1.9999097720683782e+01 dt=9.0227931621811308e-04
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=7.3392564307988160e-05 Mdot_out=8.4265512750467147e-05 Mdot_net=-1.0872948442478992e-05 Msun/yr Mstar_diag=2.8837616261907137e+01 Msun (gravity fixed)
cycle=322803 time=2.0000000000000000e+01 dt=1.8045586324362262e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=322803
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 88918023168
cpu time used  = 1.5862480392499999e+05
zone-cycles/cpu_second = 5.6055560648662283e+05
end_time   = 2026-10-27T15:56:42+09:00
exit_code  = 0

## 7.結果・考察等
こちらも分裂は抑えられず、不安定な結果になった。

## 8.備考
VTKファイルはすべて削除

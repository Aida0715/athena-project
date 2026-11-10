# Toyouchi-test84

## 1. 目的
fM=0.8(もとのガス質量の80%)としてより安定な定常円盤ができるかテスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: f9614ca
変更点: fM=0.8

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
Toyouchi-test84.sh

## 6.計算ログ
cycle=505125 time=1.9997976341516612e+01 dt=2.0236584833881466e-03
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=6.1960585927547329e-02 Mdot_out=8.9137747437469518e-04 Mdot_net=6.1069208453172638e-02 Msun/yr Mstar_diag=1.1998492779293540e+03 Msun (gravity fixed)
cycle=505126 time=2.0000000000000000e+01 dt=3.0643211731635760e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=505126
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 139139987456
cpu time used  = 2.2644493691500000e+05
zone-cycles/cpu_second = 6.1445395667304576e+05
end_time   = 2026-11-07T11:23:15+09:00
exit_code  = 0

## 7.結果・考察等
fM=1のモデルよりは比較的安定した円盤が形成された。
なお本モデルとfM=1のモデルのみ、Mdot > Mdot_crit=0.04Msun/yrとなった。

## 8.備考
VTKファイルはすべて削除

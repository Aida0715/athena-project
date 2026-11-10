# Toyouchi-test83

## 1. 目的
fM=0.6(もとのガス質量の60%)としてより安定な定常円盤ができるかテスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: f9614ca
変更点: fM=0.6

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
Toyouchi-test83.sh

## 6.計算ログ
cycle=215609 time=1.9998894198356876e+01 dt=1.1058016431242379e-03
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=3.0684098133417763e-02 Mdot_out=7.6840878908510026e-04 Mdot_net=2.9915689344332662e-02 Msun/yr Mstar_diag=7.4940464900484722e+02 Msun (gravity fixed)
cycle=215610 time=2.0000000000000000e+01 dt=1.2112503829925883e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=215610
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 59391068160
cpu time used  = 8.8577142705999999e+04
zone-cycles/cpu_second = 6.7050106094669714e+05
end_time   = 2026-11-04T17:54:07+09:00
exit_code  = 0

## 7.結果・考察等
安定な円盤は形成されたが、Mdot < Mdot_critであった。

## 8.備考
VTKファイルはすべて削除

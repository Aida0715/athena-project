# Toyouchi-test85

## 1. 目的
fM=0.4(もとのガス質量の40%)、Bz0=10μGとして磁場の影響を確かめる比較テスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: f9614ca
変更点: fM=0.4、Bz0=10μG

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
Toyouchi-test85.sh

## 6.計算ログ
cycle=245156 time=1.9999949755628684e+01 dt=5.0244371315955050e-05
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=3.6447832487826320e-02 Mdot_out=6.1353236657601621e-03 Mdot_net=3.0312508822066157e-02 Msun/yr Mstar_diag=5.1952497646093002e+02 Msun (gravity fixed)
cycle=245157 time=2.0000000000000000e+01 dt=6.0888162901617384e-05
Terminating on time limit
time=2.0000000000000000e+01 cycle=245157
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 67529966592
cpu time used  = 1.1822388505700001e+05
zone-cycles/cpu_second = 5.7120408925354946e+05
end_time   = 2026-11-09T06:26:22+09:00
exit_code  = 0

## 7.結果・考察等
磁場が全く影響していないようである。

## 8.備考
VTKファイルはすべて削除

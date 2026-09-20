# Toyouchi-test50

## 1. 目的
最適な拡散係数を決定するためのテスト。自己重力はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: Bz=3μG、eta_ohm=0.01。

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
  #Self-Gravity:                 OFF
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
Toyouchi-test50.sh

## 6.計算ログ
cycle=25599 time=4.3858609398244909e+01 dt=1.3906017550908700e-03
Sink: Mstar=1.2287292187787983e+02 Mdot_flux=1.9125527942203109e+00 dM_flux=2.6595992723467119e-03 Mdot_reset=1.9125528029302890e+00 Mdot_floor=8.7099780378406075e-09 Msink_gas=4.6626562499999965e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=7.0233102957878302e+01 Bmax_sink=7.5643366868924098e-03
cycle=25600 time=4.3859999999999999e+01 dt=2.7812035101817401e-03
Terminating on time limit
time=4.3859999999999999e+01 cycle=25600
tlim=4.3859999999999999e+01 nlim=-1
Number of MeshBlocks = 512; 896  created, 448 destroyed during this simulation.
zone-cycles = 54691987456
cpu time used  = 1.4052247399000000e+04
zone-cycles/cpu_second = 3.8920455855262019e+06
end_time   = 2026-09-09T23:52:46+09:00
exit_code  = 0

## 7.結果・考察等
磁場の拡散が弱く、40数kyrで磁場が円盤のミッドプレーンに溜まってしまう。

## 8.備考
VTKファイルはすべて削除

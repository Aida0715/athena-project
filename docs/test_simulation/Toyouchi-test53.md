# Toyouchi-test53

## 1. 目的
最適な拡散係数を決定するためのテスト。自己重力はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: Bz=3μG、eta_ohm=0.3。

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
Toyouchi-test53.sh

## 6.計算ログ
cycle=23684 time=4.3859820710183946e+01 dt=1.7928981605308536e-04
Sink: Mstar=1.2290624053098649e+02 Mdot_flux=2.6187415983545108e+00 dM_flux=4.6951369945954298e-04 Mdot_reset=2.6187416024452195e+00 Mdot_floor=4.0907089570419216e-09 Msink_gas=4.6626562499999965e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=6.8782039819478754e+01 Bmax_sink=7.4080524039756221e-03
cycle=23685 time=4.3859999999999999e+01 dt=3.5857963210617072e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=23685
tlim=4.3859999999999999e+01 nlim=-1
Number of MeshBlocks = 512; 896  created, 448 destroyed during this simulation.
zone-cycles = 50564694016
cpu time used  = 1.5167777762000000e+04
zone-cycles/cpu_second = 3.3336916461606054e+06
end_time   = 2026-09-11T02:00:06+09:00
exit_code  = 0

## 7.結果・考察等
なぜかeta_ohm=0.1より磁場の拡散が弱く、test50、52と同様に40数kyrで磁場が円盤のミッドプレーンに溜まってしまう。

## 8.備考
VTKファイルはすべて削除

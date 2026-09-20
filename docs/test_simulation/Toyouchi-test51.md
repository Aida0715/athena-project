# Toyouchi-test51

## 1. 目的
最適な拡散係数を決定するためのテスト。自己重力はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: Bz=3μG、eta_ohm=0.1。

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
Toyouchi-test51.sh

## 6.計算ログ
cycle=20552 time=4.3857354873101670e+01 dt=2.6451268983294085e-03
Sink: Mstar=1.2309626200777750e+02 Mdot_flux=1.4823577873355083e+00 dM_flux=3.9210244562292180e-03 Mdot_reset=1.4823577993184374e+00 Mdot_floor=1.1982929052197717e-08 Msink_gas=4.6626562499999965e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=6.0596669395669103e+01 Bmax_sink=6.5264610291824341e-03
cycle=20553 time=4.3859999999999999e+01 dt=3.8915112897754938e-03
Terminating on time limit
time=4.3859999999999999e+01 cycle=20553
tlim=4.3859999999999999e+01 nlim=-1
Number of MeshBlocks = 512; 896  created, 448 destroyed during this simulation.
zone-cycles = 43701436416
cpu time used  = 1.1199397704999999e+04
zone-cycles/cpu_second = 3.9021238076481009e+06
end_time   = 2026-09-10T05:16:53+09:00
exit_code  = 0

## 7.結果・考察等
eta_ohmを決めるtest50~53の４モデルのうち、計算時間含め最も効率化できた。
今後の計算では、特に問題がない限りeta_ohm=0.1を使用することにする。

## 8.備考
VTKファイルはすべて削除

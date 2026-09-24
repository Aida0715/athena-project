# Toyouchi-test59

## 1. 目的
磁場なし、自己重力OFF。test57との比較用


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: Bz=0、自己重力OFF

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
Toyouchi-test59.sh

## 6.計算ログ
cycle=1507 time=4.3856350619596249e+01 dt=3.6493804037505129e-03
Sink: Mstar=3.9170357573270294e+01 Mdot_flux=6.7812061374553223e-01 dM_flux=2.4747200791822160e-03 Mdot_reset=9.6918799569857572e-01 Mdot_floor=2.9106738195304832e-01 Msink_gas=2.6146890000000056e+00 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=1508 time=4.3859999999999999e+01 dt=7.2987608075010257e-03
Terminating on time limit
time=4.3859999999999999e+01 cycle=1508
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 415387648
cpu time used  = 1.6641996300000000e+02
zone-cycles/cpu_second = 2.4960205525343134e+06
end_time   = 2026-09-22T16:04:38+09:00
exit_code  = 0

## 7.結果・考察等
自己重力OFFにすると計算はかなり早くなるが、test57との比較により、中心星質量進化・降着率・円盤形成・スパイラル形成に早い段階で差が出ることがわかった。
特に中心星質量進化と構造形成には明確な差が見られた。
計算時間節約のため円盤質量〜中心星質量となるまでは自己重力OFFの設定にすることも考えたが、以上の結果を考慮すると、はじめから自己重力はONにすべきと思われる。

## 8.備考
VTKファイルはすべて削除

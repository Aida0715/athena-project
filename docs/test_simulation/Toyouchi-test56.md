# Toyouchi-test56

## 1. 目的
磁場なし、回転のみ、自己重力をON。密度フロア、シンクフロアともに1e-6に上げ、test55で見られた速度爆発からのタイムステップ持っていかれる問題が解決するかテスト。


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: 密度フロア、シンクフロアともに1e-6に上げた

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
Toyouchi-test56.sh

## 6.計算ログ
cycle=151127 time=4.3859969174964135e+01 dt=3.0825035864268102e-05
Sink: Mstar=2.7293430050170366e+02 Mdot_flux=1.6785051180453510e+00 dM_flux=5.1739980462105507e-05 Mdot_reset=1.6785051219665907e+00 Mdot_floor=0.0000000000000000e+00 Msink_gas=4.1502999999999999e-05 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=151128 time=4.3859999999999999e+01 dt=6.1650071728536204e-05
Terminating on time limit
time=4.3859999999999999e+01 cycle=151128
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 41629114368
cpu time used  = 7.2479694459999999e+04
zone-cycles/cpu_second = 5.7435554437904293e+05
end_time   = 2026-09-19T17:21:32+09:00
exit_code  = 0

## 7.結果・考察等
計算時間短縮はかなり改善した。
降着率はtest55とほとんど変わらず、~0.01M☀/yr であった。

## 8.備考
VTKファイルはすべて削除

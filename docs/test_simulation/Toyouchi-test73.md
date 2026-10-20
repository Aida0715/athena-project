# Toyouchi-test73

## 1. 目的
これまでの計算で円盤が初期で不安定な挙動であったので、Toyouchi+23のfig.2を根拠に音速をT=530K→10^3Kの値に変更し、安定性を比較。
新単位系に変更し、シンク半径を10^3AUに戻した。
（そもそもToyouchi+23の初期条件がr<10^3AUの範囲は対応していないため）

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: d8209c3
変更点: 温度を530K→10^3K、新単位系に変更、シンク半径をもとに戻した

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
Toyouchi-test73.sh

## 6.計算ログ
cycle=40547 time=5.0173162796914788e+01 dt=4.4291076214619807e-05
Sink: Mstar=2.8051121842458156e+02 Mdot_flux=2.2725626102571490e+00 dM_flux=1.0065424377339472e-04 Mdot_reset=2.6446153242440968e+00 Mdot_floor=3.7205271398693635e-01 Msink_gas=2.5342574765707417e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=40548 time=5.0173207087991003e+01 dt=8.8582152429239613e-05
Terminating on time limit
time=5.0173207087991003e+01 cycle=40548
tlim=5.0173207087991003e+01 nlim=-1
zone-cycles = 11169189888
cpu time used  = 1.9069121771999999e+04
zone-cycles/cpu_second = 5.8572125248055183e+05
end_time   = 2026-10-14T00:21:43+09:00
exit_code  = 0

## 7.結果・考察等
安定性はほぼ変わらない。
温度はもっと上げたほうが良いか？

## 8.備考
VTKファイルはすべて削除

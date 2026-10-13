# Toyouchi-test72

## 1. 目的
これまでの計算で、シンク境界でリングのような構造物ができ、それが分裂して不安定を起こしているようであった。
それが人工的なものなのか確かめるため、シンク半径を小さくして実験。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 57034e1
変更点: シンク半径を1000AU→500AU

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
Toyouchi-test72.sh

## 6.計算ログ
cycle=421765 time=4.3859740202658983e+01 dt=2.5979734101611029e-04
Sink: Mstar=2.7027410942167290e+02 Mdot_flux=1.2036922664256080e+00 dM_flux=3.1271605021902836e-04 Mdot_reset=1.3410444371956074e+00 Mdot_floor=1.3735217076999753e-01 Msink_gas=3.6735300000000010e-01 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=421766 time=4.3859999999999999e+01 dt=5.1959468203222059e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=421766
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 116177975296
cpu time used  = 2.0316607725199999e+05
zone-cycles/cpu_second = 5.7183746847608301e+05
end_time   = 2026-10-10T23:42:55+09:00
exit_code  = 0

## 7.結果・考察等
シンクを小さくしたが、やはりシンク境界でリング構造が生じる。
おそらくシンク内部での、密度フロア一定・速度除去などの不連続処理が原因でshockが生じることが原因ではないかと考えられる。
指導教員が解決法を調べてくれるとのことなので、お願いしてみることにする。
なお初期条件としているToyouchi+23のfig.2では、10^3AUより内側のプロファイルは存在しない（Toyouchi+23のシンク半径が10^3AUであるため）。
なのでシンク半径を500AUに縮めた場合、500<r<10^3AUの範囲では密度はr=10^3AUでの値がコピーされ、ガスに速度は与えられていない。
従って、どっちみちもとの10^3AUに戻す必要がある。

## 8.備考
VTKファイルはすべて削除

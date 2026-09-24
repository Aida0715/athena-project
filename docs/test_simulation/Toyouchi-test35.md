# Toyouchi-test35

## 1. 目的
Toyouchi+23の設定に磁場Bz=3μGを入れてテスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: cd2fe2c
変更点: Toyouchi+23の設定に磁場Bz=3μGを入れてテスト

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlle
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
Toyouchi-test35.sh

## 6.計算ログ
cycle=95646 time=2.9272281969151447e+01 dt=4.2635030243348048e-04
Sink: Mstar = 5.5155190122668905e+01, Mdot = 6.5809127718160720e-01, dM = 2.8057741505521368e-04
cycle=95647 time=2.9272708319453880e+01 dt=3.7111638341423244e-04
Sink: Mstar = 5.5155437342926902e+01, Mdot = 6.6615290794652793e-01, dM = 2.4722025799798954e-04
cycle=95648 time=2.9273079435837293e+01 dt=7.4223276682846488e-04
^Cend_time   = 2026-08-23T21:10:26+09:00
exit_code  = 1
ERROR: Athena++/MPI terminated abnormally. See /work/beta/aida/results/Toyouchi-test35/Toyouchi-test35.20260806-141534.log

## 7.結果・考察等
途中で止めた。
本データは解析対象外。記録だけ残す。

## 8.備考
VTKファイルはすべて削除

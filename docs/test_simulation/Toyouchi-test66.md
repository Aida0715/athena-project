# Toyouchi-test66

## 1. 目的
クランプがメッシュブロック由来なのかどうかを確かめるため、SMRをOFF、磁場なし、シンクなし、50kyr


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 8331c2e
変更点: シンクなし、50kyr

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
Toyouchi-test66.sh

## 6.計算ログ
cycle=74 time=4.3650516231453253e+01 dt=2.0948376854674677e-01
cycle=75 time=4.3859999999999999e+01 dt=4.0685804493485106e-01
Terminating on time limit
time=4.3859999999999999e+01 cycle=75
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 4800000
cpu time used  = 4.4763385000000000e+01
zone-cycles/cpu_second = 1.0723049653193117e+05
end_time   = 2026-09-28T20:48:34+09:00
exit_code  = 0

## 7.結果・考察等
リング形成すら見られなかった。
意味のないテストだったかもしれない。

## 8.備考
VTKファイルはすべて削除

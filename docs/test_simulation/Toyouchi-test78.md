# Toyouchi-test78

## 1. 目的
磁場を0.1μGとし、比較計算用

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: b3b903d
変更点: Bz0=0.1μG

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
Toyouchi-test78.sh

## 6.計算ログ
cycle=254152 time=1.9994993190014185e+01 dt=1.9407603835963740e-03
cycle=254153 time=1.9996933950397782e+01 dt=3.0660496022179018e-03
cycle=254154 time=2.0000000000000000e+01 dt=4.9778831428276072e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=254154
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 70008244224
cpu time used  = 1.2487056534000000e+05
zone-cycles/cpu_second = 5.6064649049501936e+05
end_time   = 2026-10-25T04:52:52+09:00
exit_code  = 0

## 7.結果・考察等
磁場を入れたが、クランプの分裂は抑えられず結果は同じ。

## 8.備考
VTKファイルはすべて削除

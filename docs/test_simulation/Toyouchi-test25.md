# Toyouchi-test25

## 1. 目的
NFWを追加して、実装状況を確認する。
自己重力・中心星重力はOFFにした。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 5e17b2d
変更点: NFWを加速度に直接追加

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlle
  #Magnetic fields:              OFF
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
  #MPI parallelism:              OFF
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test25.sh

## 6.計算ログ
Terminating on time limit
time=1.0000000000000000e+00 cycle=1
tlim=1.0000000000000000e+00 nlim=50000
Number of MeshBlocks = 288; 280  created, 0 destroyed during this simulation.
zone-cycles = 9437184
cpu time used  = 4.7181959999999998e+00
zone-cycles/cpu_second = 2.0001678607671238e+06

## 7.結果・考察等
自己重力と中心星重力はOFFにして、NFWが問題なく実装できているか確認した。
NFW加速度を直接VTKに出力するよう改造し出力データをプロファイルにしたところ、理論曲線とデータの曲線が全半径で
一致しており、問題なく実装できた。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test19

## 1. 目的
領域をToyouchi+23に合わせて、10^7auに拡張
AMR、回転(vr,vphi)、中心星重力、自己重力は完全OFF

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a2f95b5
変更点: 領域を10^7auに拡張

## 4. ビルド設定（configure）
　#Problem generator:          Toyouchi
  #Coordinate system:          cartesian
  #Equation of state:          isothermal
  #Riemann solver:             hlle
  #Magnetic fields:            OFF
  #Relativistic dynamics:      OFF 
  #General relativity:         OFF 
  #Radiative Transfer:         OFF
  #Implicit Radiation:         OFF
  #Cosmic Ray Transport:       OFF
  #Cosmic Ray Diffusion:       OFF
  #Self-Gravity:               OFF
  #Super-Time-Stepping:        OFF
  #Floating-point precision:   double
  #Number of ghost cells:      4
  #MPI parallelism:            OFF
  #OpenMP parallelism:         OFF
  #FFT:                        OFF
  #HDF5 output:                OFF
  #Compiler:                   g++
  #Compilation command:        g++ -O3 -std=c++11

## 5. 対応run
Toyouchi-test19.sh

## 6.計算ログ


## 7.結果・考察等

## 8.備考
VTKファイルはすべて削除

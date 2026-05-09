# Toyouchi-test18

## 1. 目的
ラディアル方向の速度はtest17でVr=0.85(弱い外向き流)で導入していた。
今回はより物理的と考えられるvr<0(弱い降着流)に変更して挙動を確認。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 7cef4da
変更点: Vr=0.85を導入

## 4. ビルド設定（configure）
　#Problem generator:          jeans
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
  #Self-Gravity:               Multigrid
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
Toyouchi-test18.sh

## 6.計算ログ
Terminating on time limit
time=2.0000000000000000e+01 cycle=235
tlim=2.0000000000000000e+01 nlim=50000
Number of MeshBlocks = 568; 56  created, 0 destroyed during this simulation.
zone-cycles = 4373872640
cpu time used  = 1.0754128728000000e+04
zone-cycles/cpu_second = 4.0671566712903127e+05

## 7.結果・考察等

## 8.備考
VTKファイルはすべて削除

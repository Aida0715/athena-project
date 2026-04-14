# Toyouchi-test17

## 1. 目的
ラディアル方向の速度をToyouchi+23 fig.2に合わせて、Vr=0.85で導入し、挙動をみる。

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
Toyouchi-test17.sh

## 6.ログ記録
Terminating on time limit
time=2.0000000000000000e+01 cycle=183
tlim=2.0000000000000000e+01 nlim=50000
Number of MeshBlocks = 568; 56  created, 0 destroyed during this simulation.
zone-cycles = 3406036992
cpu time used  = 8.1328504690000000e+03
zone-cycles/cpu_second = 4.1879990354953619e+05

## 7.結果・考察等
Vr>0として外向き流を入れている。VφのみのToyouchi-test16よりは強めの外向き流が生じ、期待通り。
結果を総合すると、本テストで物理的な破綻などの問題は見られなかった。
ただし、初期条件としてのVrの符号の正・負は要確認。

## 8.備考
VTKファイルはすべて削除

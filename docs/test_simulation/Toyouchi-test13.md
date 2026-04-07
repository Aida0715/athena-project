# Toyouchi-test13

## 1. 目的
Toyouchi-test12でコアダンプの原因がAMRにあることが判明。
本計算ではさらにその原因がjeansリファインと密度リファインのいずれにあるのか、jeansリファインのみ有効にしてテスト。
（Cells per Jeans length = 8）

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 9038d69
変更点: jeansリファインのみ有効

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
Toyouchi-test13.sh

## 6.ログ記録
Terminating on time limit
time=2.0000000000000000e+01 cycle=109
tlim=2.0000000000000000e+01 nlim=50000
Number of MeshBlocks = 512; 0  created, 0 destroyed during this simulation.
zone-cycles = 1828716544
cpu time used  = 3.8693345220000001e+03
zone-cycles/cpu_second = 4.7261784516236768e+05

## 7.結果・考察等
ログにある通り、AMRは発火なし。
これは中心よりmassが外側へ拡散（密度プロファイルを見ればわかる）し密度が低下。
するとjeans長が伸び、リファインする必要がなくなるためと思われる。
今回AMRの暴走等でコアダンプしなかったのは、そもそも発火していないからである。

## 8.備考
VTKファイルはすべて削除

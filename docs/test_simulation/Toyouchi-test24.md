# Toyouchi-test24

## 1. 目的
計算時間1Myrで回して、重力崩壊の時間発展を確認する。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a2f95b5
変更点: ファイル出力間隔dt=10,tlim=1000

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
  #MPI parallelism:              OFF
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test24.sh

## 6.計算ログ
Terminating on Interrupt signal
time=5.7146855982527802e+02 cycle=7271
tlim=1.0000000000000000e+03 nlim=50000
Number of MeshBlocks = 288; 280  created, 0 destroyed during this simulation.
zone-cycles = 68617764864
cpu time used  = 1.6633584816299999e+05
zone-cycles/cpu_second = 4.1252541542793805e+05

## 7.結果・考察等
numlevel=6に変更したところ、レベル5までAMRは生成された。
ただ、dt=10、tlim=1000では計算終了までに3日くらい掛かりそうだったため、57ファイルめで計算ストップさせた。
密度はr~18あたりで約1000倍に潰れることが確認できた。
密度mapの解析が途中。

## 8.備考
VTKファイルはすべて削除

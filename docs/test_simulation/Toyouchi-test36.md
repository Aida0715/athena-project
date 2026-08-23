# Toyouchi-test36

## 1. 目的
シンク内部に磁場処理を入れた

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 909df56
変更点: シンク内部に磁場処理入れた

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
Toyouchi-test36.sh

## 6.計算ログ
Terminating on time limit
time=2.0000000000000000e+01 cycle=5244
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 512; 504  created, 0 destroyed during this simulation.
zone-cycles = 60022456320
cpu time used  = 1.1411429600200000e+05
zone-cycles/cpu_second = 5.2598542358748836e+05
end_time   = 2026-08-22T05:15:27+09:00
exit_code  = 0

## 7.結果・考察等
シンク内部には磁場処理を入れておらず、内部でアルフベン速度が上昇しタイムステップが持ってかれていた。
シンク内部のアルフベン速度に上限をつけることで計算時間の短縮を試みた。
解析した限り、問題なく実装はできていると思われる。

## 8.備考
VTKファイルはすべて削除

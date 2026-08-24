# Toyouchi-test37

## 1. 目的
タイムステップと計算時間節約のため、計算領域縮小＆sink_va_cap下げ＆AMRレベル下げ
問題ないかテスト計算

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 002ad4a
変更点: 計算領域をx,y各方向半分(224)、sink_vA_cap=5、numlevel=5

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
Toyouchi-test37.sh

## 6.計算ログ
Terminating on time limit
time=2.0000000000000000e+01 cycle=775
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 1472; 1596  created, 140 destroyed during this simulation.
zone-cycles = 1956331520
cpu time used  = 3.5616210200000000e+02
zone-cycles/cpu_second = 5.4928121465320867e+06
end_time   = 2026-08-23T21:29:33+09:00
exit_code  = 0

## 7.結果・考察等
計算時間は大幅に短縮できた。
ただ、中心星質量成長や降着率が数倍差ができた。
これが領域サイズの問題なのか解像度の問題なのか切り分ける必要がありそう。
test38では、領域サイズのみをもとに戻しテストしてみる。

## 8.備考
VTKファイルはすべて削除

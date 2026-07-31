# Toyouchi-test29

## 1. 目的
Bz=3μGに強めて安定性テスト。AMRと回転はOFFにした。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 3cc79ae
変更点: MHD計算できるようにした

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
  #MPI parallelism:              OFF
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test29.sh

## 6.計算ログ
cycle=16 time=1.0000000000000000e+03 dt=1.8155913621538048e+01
Terminating on time limit
time=1.0000000000000000e+03 cycle=16
tlim=1.0000000000000000e+03 nlim=50000
Number of MeshBlocks = 8; 0  created, 0 destroyed during this simulation.
zone-cycles = 4194304
cpu time used  = 1.2102976999999999e+01
zone-cycles/cpu_second = 3.4655143110657821e+05

## 7.結果・考察等
磁場をBz＝１μG→３μGに強めてテストした。
Bmaxは3μG→33.56μG（約11倍）となった。Test28でも1μG→11.19μGと約11倍になった。
つまり磁場を強めても計算がおかしくなることはなく、物理が保たれていることがわかる。
またプラズマβはβ∝B^-2であるため、Bz＝1μGのときと3μGで比較すると、βmin_3μG＝βmin_1μG / 9 となるはず。
実際Test28の計算値βmin_1μG~14を代入して計算すると、βmin_3μG~1.56となりTest29の計算値におおよそ一致する。
したがって、本テストは問題なくできたと考えられる。

## 8.備考
VTKファイルはすべて削除

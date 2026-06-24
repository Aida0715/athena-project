# Toyouchi-test23

## 1. 目的
test22で潰れて中心付近密度が約100倍になったがAMRはレベル1までしか作動していない。
そこでAMRの動作条件を、jeans_cells=16.0、refine_thr=0.2に厳し目に設定し、動作状況の条件依存性を確認。 
回転(vr,vphi)はOFF

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a2f95b5
変更点: jeans_cells=16.0、refine_thr=0.2

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
Toyouchi-test23.sh

## 6.計算ログ
Terminating on time limit
time=1.0000000000000000e+03 cycle=67
tlim=1.0000000000000000e+03 nlim=50000
Number of MeshBlocks = 64; 56  created, 0 destroyed during this simulation.
zone-cycles = 140509184
cpu time used  = 3.1171634899999998e+02
zone-cycles/cpu_second = 4.5075975145596231e+05

## 7.結果・考察等
AMRの動作条件を変え、test22と同様のテストを行った。
中心付近密度が潰れて約100倍になっているにも関わらず、AMRはレベル1までしか作動しなかった。
これは入力ファイルで、numlevel=2（レベル0+1くらいまでしか生成されない？）となっているのが原因と考えて、numlevel=6に変更し再度同じテストをやってみる。

## 8.備考
VTKファイルはすべて削除

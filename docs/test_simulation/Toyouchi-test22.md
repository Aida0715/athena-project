# Toyouchi-test22

## 1. 目的
AMRはjeansリファイン・密度勾配両方ONにして、AMRの動作及び解像度の確認
中心星重力・自己重力ともにONにして潰し、AMRの動作状況を確認する
回転(vr,vphi)はOFF

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a2f95b5
変更点: 中心星重力・自己重力ON

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
Toyouchi-test22.sh

## 6.計算ログ
Terminating on time limit
time=1.0000000000000000e+03 cycle=67
tlim=1.0000000000000000e+03 nlim=50000
Number of MeshBlocks = 64; 56  created, 0 destroyed during this simulation.
zone-cycles = 140509184
cpu time used  = 3.1174804899999998e+02
zone-cycles/cpu_second = 4.5071391609575081e+05

## 7.結果・考察等
中心星重力・自己重力両方をONにして潰し、AMRの動作状況を確認した。
結果として重力完全OFFにしたtest21と変わらず、AMRはレベル1までしか作動していなかった。
おそらく現状の動作条件（jeans_cells = 8,refine_thr = 0.3）では緩すぎるのではないかと考えられる。
一度jeans_cells = 16、refine_thr = 0.2のように条件を少しきつくして、同様のテストを行う。

## 8.備考
VTKファイルはすべて削除

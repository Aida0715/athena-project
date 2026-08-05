# Toyouchi-test29

## 1. 目的
Bz=3μG,Takasao+22に合わせて計算時間を26.3Myrにして長期安定性テスト。AMRと回転はOFFにした。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: c88f238
変更点: 計算時間を26.3Myr

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
Toyouchi-test30.sh

## 6.計算ログ
end_time   = 2026-08-05T22:32:38+09:00
exit_code  = 1
ERROR: Athena++/MPI terminated abnormally. See /work/beta/aida/results/Toyouchi-test30/Toyouchi-test30.20260801-133206.log

## 7.結果・考察等
計算時間だけTakasao+22に合わせても意味ないとのことなので、計算途中で止めた。
役に立つデータではないと思うが、念の為残しておく。

## 8.備考
VTKファイルはすべて削除

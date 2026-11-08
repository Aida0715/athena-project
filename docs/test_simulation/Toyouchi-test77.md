# Toyouchi-test77

## 1. 目的
test76で温度を10^4Kとしたところ、ガス圧が強くなりガスが拡散してしまった。
そのため温度を10^3Kに戻して安定性を確認。
磁場はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 4c342db
変更点: 温度を10^3Kに戻した

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlld
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
  #MPI parallelism:              ON
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test77.sh

## 6.計算ログ
cycle=258545 time=1.9991976472931483e+01 dt=6.2592514744896589e-03
cycle=258546 time=1.9998235724405973e+01 dt=1.7642755940272536e-03
cycle=258547 time=2.0000000000000000e+01 dt=3.5285511880545073e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=258547
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 71218322432
cpu time used  = 1.2043045701500001e+05
zone-cycles/cpu_second = 5.9136471119701490e+05
end_time   = 2026-10-22T03:42:14+09:00
exit_code  = 0

## 7.結果・考察等
ガスの拡散は抑えられたが、まだ不安定で早い段階でクランプに分裂してしまう。

## 8.備考
VTKファイルはすべて削除

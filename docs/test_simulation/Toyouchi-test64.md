# Toyouchi-test64

## 1. 目的
クランプがメッシュブロック由来なのかどうかを確かめるため、SMRをOFF、磁場なし、tlim=17.54の短時間テスト


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 8331c2e
変更点: SMRはOFF,磁場なし、tlim=17.54

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
Toyouchi-test64.sh

## 6.計算ログ
cycle=17 time=1.7191719672739850e+01 dt=3.4828032726014868e-01
Sink: Mstar=1.0000000000000000e+00 Mdot_flux=0.0000000000000000e+00 dM_flux=0.0000000000000000e+00 Mdot_reset=0.0000000000000000e+00 Mdot_floor=0.0000000000000000e+00 Msink_gas=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=18 time=1.7539999999999999e+01 dt=6.2197740545405544e-01
Terminating on time limit
time=1.7539999999999999e+01 cycle=18
tlim=1.7539999999999999e+01 nlim=-1
zone-cycles = 1152000
cpu time used  = 3.5297429999999999e+00
zone-cycles/cpu_second = 3.2636937023460347e+05
end_time   = 2026-09-28T18:56:27+09:00
exit_code  = 0

## 7.結果・考察等
SMRをOFFにしたユニフォームグリッドだと解像度が足りなすぎてよくわからなかった。レゾリューション上げて再度実行。

## 8.備考
VTKファイルはすべて削除

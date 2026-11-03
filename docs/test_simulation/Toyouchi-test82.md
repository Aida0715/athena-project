# Toyouchi-test82

## 1. 目的
fM=0.4(もとのガス質量の40%)としてより安定な定常円盤ができるかテスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: f9614ca
変更点: fM=0.4

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
Toyouchi-test82.sh

## 6.計算ログ
cycle=126300 time=1.9999148207931682e+01 dt=8.5179206831753618e-04
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=2.5020375540968073e-02 Mdot_out=1.1317190946014933e-02 Mdot_net=1.3703184594953139e-02 Msun/yr Mstar_diag=5.2443236876513356e+02 Msun (gravity fixed)
cycle=126301 time=2.0000000000000000e+01 dt=1.4251136963379164e-03
Terminating on time limit
time=2.0000000000000000e+01 cycle=126301
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 34790368256
cpu time used  = 6.0152017593999997e+04
zone-cycles/cpu_second = 5.7837408698108653e+05
end_time   = 2026-11-03T13:15:54+09:00
exit_code  = 0

## 7.結果・考察等
今回からWise+19のガス分布をfMをファクターとして分布を保ったままその量を変更できるようにした。
安定した円盤になった。
20kyrではクランプに分裂することもなく、またtest81(fM=1相当)よりも円盤は小さく、スパイラル構造もきれいにできた。

## 8.備考
VTKファイルはすべて削除

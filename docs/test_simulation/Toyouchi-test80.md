# Toyouchi-test80

## 1. 目的
音速を10^3K→5e3Kに上げ、BH質量を２M☀→100M☀に変更して定常円盤ができるかテスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a4554e8
変更点: 温度を5e3、BH質量を100M☀

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
Toyouchi-test80.sh

## 6.計算ログ
cycle=1241 time=1.9988754231458699e+01 dt=1.1245768541300549e-02
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=5.6392555581647336e-05 Mdot_out=1.7266469897581913e-07 Mdot_net=5.6219890882671519e-05 Msun/yr Mstar_diag=1.7800632696131197e+02 Msun (gravity fixed)
cycle=1242 time=2.0000000000000000e+01 dt=1.5771551500635431e-02
Terminating on time limit
time=2.0000000000000000e+01 cycle=1242
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 342116352
cpu time used  = 7.0061476300000004e+02
zone-cycles/cpu_second = 4.8830879688442987e+05
end_time   = 2026-10-27T18:24:36+09:00
exit_code  = 0

## 7.結果・考察等
圧力勾配力が強すぎたためか、初期からほとんどaccretionしないうえに、途中からガスが外側に拡散した。
温度をもとの10^3Kに下げて再計算してみる必要がある。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test74

## 1. 目的
シンク設計を改良して、境界でリング状構造ができないか確認する。
磁場はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: c79c93d
変更点: シンク設計を改良

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
Toyouchi-test74.sh

## 6.計算ログ
cycle=73419 time=5.0173059446192724e+01 dt=1.4764179827864155e-04
Sink(drain): Mstar=2.3469225002004740e+02 dM_star=2.9623023447927619e-04 dM_flux=3.6847379637350373e-04 Msink_gas=4.3157181625136420e-01 tau=2.3233428976347503e-02
cycle=73420 time=5.0173207087991003e+01 dt=2.9528359655728309e-04
Terminating on time limit
time=5.0173207087991003e+01 cycle=73420
tlim=5.0173207087991003e+01 nlim=-1
zone-cycles = 20223979520
cpu time used  = 3.4834481282000001e+04
zone-cycles/cpu_second = 5.8057358042102738e+05
end_time   = 2026-10-15T07:02:40+09:00
exit_code  = 0

## 7.結果・考察等
実装自体はできていたが、やはり改良前と比較しても、ほぼ同時刻にリングができて分裂しクランプができてしまう。
シンク自体廃止する方向で進めた方が良さそうなので、その方針で。

## 8.備考
VTKファイルはすべて削除

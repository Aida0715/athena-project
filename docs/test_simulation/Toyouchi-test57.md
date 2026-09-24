# Toyouchi-test57

## 1. 目的
磁場なし、回転のみ、自己重力をON。シンクフロアの値をWise+19のfig.4のr<1e3AUでフィッティングして平均を取り、6.3e-2(code unit)とした。
密度フロアは1e-6で変更なし。


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: シンクフロアの値を6.3e-2とした

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
Toyouchi-test57.sh

## 6.計算ログ
cycle=120373 time=4.3859832656755934e+01 dt=1.6734324406542100e-04
Sink: Mstar=2.3907928418936703e+02 Mdot_flux=1.6392755223900952e+00 dM_flux=2.7432168383379623e-04 Mdot_reset=2.0548823387597932e+00 Mdot_floor=4.1560681636961272e-01 Msink_gas=2.6146890000000056e+00 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=120374 time=4.3859999999999999e+01 dt=2.6082309440997830e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=120374
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 33157740544
cpu time used  = 5.7777949762999997e+04
zone-cycles/cpu_second = 5.7388226269727631e+05
end_time   = 2026-09-20T13:24:48+09:00
exit_code  = 0

## 7.結果・考察等
シンクフロアの値はこれまでの1e-8、1e-6のように適当に決めていたが、今回Wise+19のfig.4をフィッティングして6.3e-2とした。これでシンクフロアの値に設定根拠をもたせた。
test56よりも計算時間は20%ほど短縮できた。降着率はtest55、test56とほぼ変わらず~0.04M☀/yrであった。
今後はシンク密度フロア値はこの値で計算する。

## 8.備考
VTKファイルはすべて削除

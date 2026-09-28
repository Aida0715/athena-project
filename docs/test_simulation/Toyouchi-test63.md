# Toyouchi-test63

## 1. 目的
比較のためBz=3μGで入れた


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 8331c2e
変更点: Bz=3μG

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
Toyouchi-test63.sh

## 6.計算ログ
cycle=104000 time=4.3859170575562125e+01 dt=8.2942443787459297e-04
Sink: Mstar=2.5142264626147926e+02 Mdot_flux=5.0085543556648195e+00 dM_flux=4.1542173810116372e-03 Mdot_reset=5.3959530905754018e+00 Mdot_floor=3.8739873491058074e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=2.2374406836396965e-02 Bmax_sink=5.6159315452272256e-03
cycle=104001 time=4.3859999999999999e+01 dt=1.2692307823841426e-03
Terminating on time limit
time=4.3859999999999999e+01 cycle=104001
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 28647699456
cpu time used  = 5.2025236152999998e+04
zone-cycles/cpu_second = 5.5065006090026279e+05
end_time   = 2026-09-28T06:37:08+09:00
exit_code  = 0

## 7.結果・考察等
test57,58,60,61,62,63のうち最も初期磁場強度が強いはずだが、降着率や中心星質量進化はBz＝0(test57)とほとんど変わらないという結果になった。
またこれまで上記のテストで与えていた初期磁場強度は天体形成前の種磁場強度であり、初期に種BHが形成されている本研究の設定としては弱すぎるという指摘があった。そのため、次回以降天体形成中に増幅があった前提で、10μG付近の値でテスト・比較してみる。

## 8.備考
VTKファイルはすべて削除

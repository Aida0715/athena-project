# Toyouchi-test69

## 1. 目的
Bz=5μGとして比較計算

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 57034e1
変更点: Bz=5μG

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
Toyouchi-test69.sh

## 6.計算ログ
cycle=134750 time=4.3858730059838699e+01 dt=1.2699401613005534e-03
Sink: Mstar=2.6899295633935520e+02 Mdot_flux=8.4493891308661553e+00 dM_flux=1.0730218595743309e-02 Mdot_reset=8.8680532047028908e+00 Mdot_floor=4.1866407383673560e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=3.3406985187840754e-01 Bmax_sink=8.3850867341046936e-02
cycle=134751 time=4.3859999999999999e+01 dt=1.5043070793782044e-03
Terminating on time limit
time=4.3859999999999999e+01 cycle=134751
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 37117971456
cpu time used  = 6.8493618688999995e+04
zone-cycles/cpu_second = 5.4191868040345062e+05
end_time   = 2026-10-02T16:28:42+09:00
exit_code  = 0

## 7.結果・考察等
現在までのテスト計算では、磁場強度が高いほど中心星質量進化が抑えられるようである。
詳細は解析してみないと分からないが、おそらく50kyr程度では磁気張力よりも磁気圧が卓越し、降着を抑えている可能性がある。

## 8.備考
VTKファイルはすべて削除

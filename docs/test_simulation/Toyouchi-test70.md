# Toyouchi-test70

## 1. 目的
Bz=0.1μGとして比較計算

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 57034e1
変更点: Bz=0.1μG

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
Toyouchi-test70.sh

## 6.計算ログ
cycle=178148 time=4.3859946280831267e+01 dt=5.3719168732868638e-05
Sink: Mstar=2.8021008517616116e+02 Mdot_flux=1.2614691333510349e+01 dM_flux=6.7765073225789811e-04 Mdot_reset=1.3111337172612998e+01 Mdot_floor=4.9664583910272453e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=1.4526510361846411e-03 Bmax_sink=3.6461251634369821e-04
cycle=178149 time=4.3859999999999999e+01 dt=1.0743833746573728e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=178149
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 49072210944
cpu time used  = 8.9612392745000005e+04
zone-cycles/cpu_second = 5.4760518540822051e+05
end_time   = 2026-10-03T17:59:07+09:00
exit_code  = 0

## 7.結果・考察等
現在までのテスト計算では、磁場強度が高いほど中心星質量進化が抑えられるようである。
詳細は解析してみないと分からないが、おそらく50kyr程度では磁気張力よりも磁気圧が卓越し、降着を抑えている可能性がある。

## 8.備考
VTKファイルはすべて削除

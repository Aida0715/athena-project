# Toyouchi-test67

## 1. 目的
カーテシアン格子由来と考えられるクランプ形成を抑制するため、White noiseとして初期密度に１％のゆらぎを与え、20kyrでテスト計算
SMRはON、磁場はOFF

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 57034e1
変更点: White noiseとして初期密度に１％のゆらぎ与えた

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
Toyouchi-test67.sh

## 6.計算ログ
cycle=57980 time=1.7538730572127594e+01 dt=1.2694278724048047e-03
Sink: Mstar=1.3558469791425426e+02 Mdot_flux=3.0000439237661012e+00 dM_flux=3.8083393752673641e-03 Mdot_reset=3.3768416215141692e+00 Mdot_floor=3.7679769774806976e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=57981 time=1.7539999999999999e+01 dt=1.4418489173574566e-03
Terminating on time limit
time=1.7539999999999999e+01 cycle=57981
tlim=1.7539999999999999e+01 nlim=-1
zone-cycles = 15971214336
cpu time used  = 2.7246251981000001e+04
zone-cycles/cpu_second = 5.8618023305140913e+05
end_time   = 2026-09-30T02:37:34+09:00
exit_code  = 0

## 7.結果・考察等
リングが不自然に四角く変形する挙動は抑えられた。
その後形成されていた４つのクランプに関しても、クランプ形成自体は見られるが、White noiseなしと比較して大分抑えられたようである。

## 8.備考
VTKファイルはすべて削除

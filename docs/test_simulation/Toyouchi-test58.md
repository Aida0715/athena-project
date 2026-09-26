# Toyouchi-test58

## 1. 目的
Latif+13を参考にBz=3e-20Gで入れて、磁場なしのtest57と比較。


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: Bz=3e-20とした

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
Toyouchi-test58.sh

## 6.計算ログ
cycle=169369 time=4.3859676681793907e+01 dt=3.2331820609243778e-04
Sink: Mstar=2.5464951743066209e+02 Mdot_flux=1.1221566567241142e+01 dM_flux=3.6281367720672813e-03 Mdot_reset=1.1498548737126725e+01 Mdot_floor=2.7698216988558721e-01 Msink_gas=2.6146890000000056e+00 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=1.5593822690800887e-15 Bmax_sink=3.9140184318759416e-16
cycle=169370 time=4.3859999999999999e+01 dt=3.5137707884468045e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=169370
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 46653982720
cpu time used  = 8.2994748391999994e+04
zone-cycles/cpu_second = 5.6213174476587796e+05
end_time   = 2026-09-24T20:21:07+09:00
exit_code  = 0

## 7.結果・考察等
降着率、中心星質量進化はBz=0の場合とほぼ一致した。

## 8.備考
VTKファイルはすべて削除

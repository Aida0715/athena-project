# Toyouchi-test61

## 1. 目的
比較のためBz=3e-20μGで入れた


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 8331c2e
変更点: Bz=3e-20μG

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
Toyouchi-test61.sh

## 6.計算ログ
cycle=103713 time=4.3859670853108824e+01 dt=3.2914689117546914e-04
Sink: Mstar=2.7774352323295079e+02 Mdot_flux=1.0844914261595481e+01 dM_flux=3.5695698142686606e-03 Mdot_reset=1.1236204216608886e+01 Mdot_floor=3.9128995501341440e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=6.9338487568374751e-22 Bmax_sink=1.7403822254636726e-22
cycle=103714 time=4.3859999999999999e+01 dt=6.5829378235093827e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=103714
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 28568643584
cpu time used  = 5.1139897482000000e+04
zone-cycles/cpu_second = 5.5863709140315477e+05
end_time   = 2026-09-26T10:56:37+09:00
exit_code  = 0

## 7.結果・考察等
test62(3e-30μG)より10桁大きいが、降着率・中心星質量進化ともにほとんど一致している。
一定以上磁場が弱すぎると、ほとんど変わらない結果になるのか？
１Myr回して進化を調べたい。

## 8.備考
VTKファイルはすべて削除

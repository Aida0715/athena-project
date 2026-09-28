# Toyouchi-test60

## 1. 目的
Latif+14を根拠にBz=3e-5μGで入れた


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 8331c2e
変更点: Bz=3e-5μG

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
Toyouchi-test60.sh

## 6.計算ログ
cycle=118423 time=4.3859675957904308e+01 dt=3.2404209569136810e-04
Sink: Mstar=2.9378752468934920e+02 Mdot_flux=6.3578562394119995e+00 dM_flux=2.0602130599235049e-03 Mdot_reset=6.7271926135388940e+00 Mdot_floor=3.6933637412690296e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=1.5214250553053623e-06 Bmax_sink=3.8187465814241750e-07
cycle=118424 time=4.3859999999999999e+01 dt=6.4808419138273621e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=118424
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 32620601344
cpu time used  = 5.8820604973000001e+04
zone-cycles/cpu_second = 5.5457779393757682e+05
end_time   = 2026-09-25T15:21:40+09:00
exit_code  = 0

## 7.結果・考察等
計算時間＜50kyrで、現時点では最も降着率が高く、中心星質量進化も大きい。
１Myr回して進化を見たい。

## 8.備考
VTKファイルはすべて削除

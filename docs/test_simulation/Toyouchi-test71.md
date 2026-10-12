# Toyouchi-test71

## 1. 目的
Bz=20μGとして比較計算

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 57034e1
変更点: Bz=20μG

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
Toyouchi-test71.sh

## 6.計算ログ
cycle=359946 time=4.3859277083510008e+01 dt=7.2291648999112113e-04
Sink: Mstar=2.3672599497293291e+02 Mdot_flux=3.0121351695760650e+00 dM_flux=2.1775221841687392e-03 Mdot_reset=3.3852810536047215e+00 Mdot_floor=3.7314588402866272e-01 Msink_gas=2.6146890000000056e+00 Va_max_sink=4.0965922685585426e-01 Bmax_sink=1.0282364988334436e-01
cycle=359947 time=4.3859999999999999e+01 dt=8.1530648101261279e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=359947
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 99149560832
cpu time used  = 1.7965396381399999e+05
zone-cycles/cpu_second = 5.5189186326360097e+05
end_time   = 2026-10-08T01:11:35+09:00
exit_code  = 0

## 7.結果・考察等
磁場強度を変更した他の計算結果とも比較して、おそらくまだ磁場が効いていない可能性がある。
中心星質量進化に50kyrで数十太陽質量の違いは出るが、これは磁場由来というよりは磁場を加えたことによる初期条件の違いに由来しているのでは？という意見もあった。とりあえず、例えばB0=0、10μGの２つで結果の違いはどこに由来するのか検証する。

## 8.備考
VTKファイルはすべて削除

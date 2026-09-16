# Toyouchi-test54

## 1. 目的
AMRをSMRに変更。磁場なし、回転のみ、自己重力はOFF。
SMR化されたかグリッド構造を調べる。


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: AMRをSMRに変更。

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
  #Self-Gravity:                 OFF
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
Toyouchi-test54.sh

## 6.計算ログ
cycle=3036 time=4.3853838827045927e+01 dt=6.1611729540729243e-03
Sink: Mstar=1.2788404204323818e+02 Mdot_flux=3.0390635515569118e+00 dM_flux=1.8724196159561252e-02 Mdot_reset=3.0390635543357782e+00 Mdot_floor=2.7788652612715258e-09 Msink_gas=4.8143480000000020e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=3037 time=4.3859999999999999e+01 dt=9.2179036964235364e-03
Terminating on time limit
time=4.3859999999999999e+01 cycle=3037
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 836559872
cpu time used  = 3.0812263500000000e+02
zone-cycles/cpu_second = 2.7150224520181711e+06
end_time   = 2026-09-15T15:08:58+09:00
exit_code  = 0

## 7.結果・考察等
グリッド構造を、メッシュ:40^3、メッシュブロック:8^3にしたところ、レベル５で回してかなり計算軽くなった。
（dx~5✕10^3AU、dxmin~156AU(rsinkを約6セルで解像)）
とりあえずこの構造で問題ないと思う。

## 8.備考
VTKファイルはすべて削除

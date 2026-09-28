# Toyouchi-test65

## 1. 目的
クランプがメッシュブロック由来なのかどうかを確かめるため、SMRをOFF、磁場なし、rsink=25000AU、tlim=17.54の短時間テスト


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 8331c2e
変更点: SMRはOFF,磁場なし、rsink=25000AU、tlim=17.54

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
Toyouchi-test65.sh

## 6.計算ログ
cycle=204 time=1.7496655484500486e+01 dt=4.3344515499512681e-02
Sink: Mstar=1.2458414741001820e+03 Mdot_flux=0.0000000000000000e+00 dM_flux=0.0000000000000000e+00 Mdot_reset=3.6730274032132792e+02 Mdot_floor=9.6025933874020950e+02 Msink_gas=3.8236520447999879e+04 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=205 time=1.7539999999999999e+01 dt=8.6689030999025363e-02
Terminating on time limit
time=1.7539999999999999e+01 cycle=205
tlim=1.7539999999999999e+01 nlim=-1
zone-cycles = 13120000
cpu time used  = 6.5445065000000000e+01
zone-cycles/cpu_second = 2.0047348107913102e+05
end_time   = 2026-09-28T19:41:15+09:00
exit_code  = 0

## 7.結果・考察等
シンク半径が大きすぎて、ガスが希薄な領域しか残らず、今回の計算時間程度ではガスリング形成まで至らなかった。
test57の密度mapにセルグリッドを重ね書きしたところ、クランプに分裂する前の不自然な四角い構造はおそらくメッシュブロック由来ではないと考える。
引き続き調べる。

## 8.備考
VTKファイルはすべて削除

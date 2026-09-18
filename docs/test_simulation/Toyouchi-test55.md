# Toyouchi-test55

## 1. 目的
磁場なし、回転のみ、自己重力をONにしてOFFの場合と計算時間、円盤形成、降着率等を比較。


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: 自己重力ON

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
Toyouchi-test55.sh

## 6.計算ログ
cycle=191948 time=4.3859183444318916e+01 dt=8.1655568108374155e-04
Sink: Mstar=2.4653527798047821e+02 Mdot_flux=4.2975152616072618e+00 dM_flux=3.5091605014094915e-03 Mdot_reset=4.2975152656956483e+00 Mdot_floor=4.0883889604802194e-09 Msink_gas=4.8143480000000020e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=191949 time=4.3859999999999999e+01 dt=7.6466832473533382e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=191949
tlim=4.3859999999999999e+01 nlim=-1
zone-cycles = 52873503744
cpu time used  = 9.0704992310000001e+04
zone-cycles/cpu_second = 5.8291723969608697e+05
end_time   = 2026-09-15T21:55:20+09:00
exit_code  = 0

## 7.結果・考察等
計算は問題なく動き、物理的にも正しそうな結果が得られた。
しかし、あるセルで密度を減少させる更新が続いてフロアに達したため速度が爆発し、タイムステップが無駄に持っていかれる問題が発生。
タイムステップは大小を変動する挙動を見せている。
ただそのセルは局所的であり、計算が壊れたり全体として非物理的な結果が出ているようには見えないので、密度フロアを大きく取ってしまえばひとまず問題は解決しそう。
密度フロアを1e-6で設定し、降着率・中心星質量成長など今回の結果と比較して影響を調べてみる。
また自己重力の有無で早い段階から降着率や中心星質量の成長に影響が出ているので、はじめから自己重力入れる方針のほうが良いかもしれない。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test49

## 1. 目的
ガス質量分布とvrの初期条件をToyouchi+23のfig.2を７次多項式フィッティングすることで見直し、実装。
Bz=0で回転のみを与えた。
なお計算時間短縮のため、自己重力はOFFにしてある（ガスの自己重力が効いてくるのは、円盤質量と中心星質量がコンパラくらいになってから）。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 33e9e2f
変更点: 初期条件見直し。AMR許可領域を2✕10^4AUに制限。密度勾配リファインをy,z方向にも拡張（以前はx方向のみ）。

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
Toyouchi-test49.sh

## 6.計算ログ
cycle=2394 time=4.3859763943441941e+01 dt=2.3605655805880588e-04
Sink: Mstar=1.2623777813778368e+02 Mdot_flux=1.5370499340780930e+00 dM_flux=3.6283071700298912e-04 Mdot_reset=1.5370499342194328e+00 Mdot_floor=1.4133933039739937e-10 Msink_gas=4.6626562499999965e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=0.0000000000000000e+00 Bmax_sink=0.0000000000000000e+00
cycle=2395 time=4.3859999999999999e+01 dt=4.7211311611761175e-04
Terminating on time limit
time=4.3859999999999999e+01 cycle=2395
tlim=4.3859999999999999e+01 nlim=-1
Number of MeshBlocks = 512; 896  created, 448 destroyed during this simulation.
zone-cycles = 4728160256
cpu time used  = 1.2015380709999999e+03
zone-cycles/cpu_second = 3.9350898403617875e+06
end_time   = 2026-09-09T18:18:32+09:00
exit_code  = 0

## 7.結果・考察等
密度mapを見る限り、初期条件はr=10^5AUでうまく与えられているようである。自己重力をOFFにしたためと思われるが、test27で見られたようなスパイラル構造は
見られない。
AMRが親ブロックをリファインしないと子ブロックをリファインできないという性質から、今回の設定では原点に接している幅~5✕10^4AUのメッシュブロック８つが
まずリファインされるので、無駄な領域までリファインすることになった。
メッシュサイズを変更して、中心付近だけSMR化するように実装し直した方が良さそうである。
円盤構造もでき、見直した初期条件(fig.2を７次多項式フィッティング)も問題なく実装できていると思われる。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test46

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: 回転あり、sink_va_cap = 0.0、sink_cfl_va_cap = 10.0、sink_cfl_cap_radius_factor = 1.0、sink_cfl_max_relax = 1.2

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlle
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
Toyouchi-test43.sh

## 6.計算ログ
cycle=3190 time=1.9999965165879591e+01 dt=3.4834120409499292e-05
Sink: Mstar=4.4222221430256845e+01 Mdot_flux=2.3433074124492914e+00 dM_flux=8.1627052561730829e-05 Mdot_reset=2.3433631081929076e+00 Mdot_floor=5.5695743616495584e-05 Msink_gas=4.6224609375000084e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=8.4368638450973037e+01 Bmax_sink=8.4368638450973039e-03
cycle=3191 time=2.0000000000000000e+01 dt=6.9668240818998584e-05
Terminating on time limit
time=2.0000000000000000e+01 cycle=3191
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 512; 448  created, 0 destroyed during this simulation.
zone-cycles = 5115166720
cpu time used  = 1.5147480277000001e+04
zone-cycles/cpu_second = 3.3769093119513028e+05
end_time   = 2026-08-26T23:47:58+09:00
exit_code  = 0

## 7.結果・考察等
物理量の時間変化はEと一致しており、dtは意図通り1.2倍に緩和。
ただし、CFL緩和は意図通りcycle数は減らしたが、総計算時間短縮にはほぼつながっていない（テストEとほぼ同じ）。
E~Gを各16コアで並走させていたので、CPU競合の影響か？
ただし、sink_Va_cap制限の用に物理を変えないので、物理面ではかなり良好と考えられる。
次は同じ設定でEとGをtlim=50まで走らせて比較。

## 8.備考
VTKファイルはすべて削除

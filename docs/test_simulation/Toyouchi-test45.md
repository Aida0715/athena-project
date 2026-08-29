# Toyouchi-test45

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: 回転あり、sink_va_cap = 10.0、sink_cfl_va_cap = 0.0、sink_cfl_cap_radius_factor = 1.0、sink_cfl_max_relax = 1.0

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
Toyouchi-test45.sh

## 6.計算ログ
cycle=2196 time=1.9999665837306633e+01 dt=3.3416269336683513e-04
Sink: Mstar=4.4135970419687652e+01 Mdot_flux=2.2518950945218426e+00 dM_flux=7.5249932996498275e-04 Mdot_reset=2.2519434993291103e+00 Mdot_floor=4.3715618087829983e-06 Msink_gas=3.1286927739731000e-06 Msink_magfloor=2.6664466802230995e-06 Va_max_sink=1.0000000000000002e+01 Bmax_sink=1.0753537264413015e-02
cycle=2197 time=2.0000000000000000e+01 dt=6.6832538673367026e-04
Terminating on time limit
time=2.0000000000000000e+01 cycle=2197
tlim=2.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 512; 448  created, 0 destroyed during this simulation.
zone-cycles = 3341975552
cpu time used  = 1.3159563654000000e+04
zone-cycles/cpu_second = 2.5395793051118139e+05
end_time   = 2026-08-26T23:19:29+09:00
exit_code  = 0

## 7.結果・考察等
他の回転系のテストと降着率は同じだが、非回転でsink_Va_capで制限したテストBと同様磁場強度などの物理量に違いが見られた。
テストBと同様、人工質量がMHD計算に影響を与えていると考えられる。

## 8.備考
VTKファイルはすべて削除

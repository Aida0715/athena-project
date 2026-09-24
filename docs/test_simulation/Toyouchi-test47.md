# Toyouchi-test47

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: test44と同様の設定で、50kyrで計算

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
Toyouchi-test47.sh

## 6.計算ログ
cycle=967784 time=4.9999972178840920e+01 dt=2.7821159079621793e-05
Sink: Mstar=1.0474770629243608e+02 Mdot_flux=5.2829891575589394e-01 dM_flux=1.4697888176836436e-05 Mdot_reset=5.2833731402465234e-01 Mdot_floor=3.7245604693982475e-05 Msink_gas=4.6224609375000084e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=2.3321563216330082e+02 Bmax_sink=2.3321563216330082e-02
cycle=967785 time=5.0000000000000000e+01 dt=5.5642318159243587e-05
Terminating on time limit
time=5.0000000000000000e+01 cycle=967785
tlim=5.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 736; 896  created, 224 destroyed during this simulation.
zone-cycles = 1856427761664
cpu time used  = 1.2465089281980000e+06
zone-cycles/cpu_second = 1.4893016164334430e+06
end_time   = 2026-09-12T14:24:44+09:00
exit_code  = 0

## 7.結果・考察等
Gは破綻しなかったが、速度向上が小さい一方、終盤の非線形イベントで積分量にも3%程度の差が生じたため、安全側として無制限のEを採用
つまり、シンク内のvAに関する人工処理は一切用いないこととする。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test48

## 1. 目的
シンク内部のアルフベン速度上限を、これまで通り人工的に密度をいじって調節すべきか、密度いじらず直接タイムスリップをいじって上限をつけるか決定するため
比較テスト計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 814a606
変更点: test46と同様の設定で、50kyrで計算

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
Toyouchi-test48.sh

## 6.計算ログ
cycle=937052 time=4.9999969365785603e+01 dt=3.0634214397196047e-05
Sink: Mstar=1.0127280028044380e+02 Mdot_flux=6.3172700454461861e-01 dM_flux=1.9352460497718286e-05 Mdot_reset=6.3175238154423297e-01 Mdot_floor=2.5376999614480246e-05 Msink_gas=4.6224609375000084e-07 Msink_magfloor=0.0000000000000000e+00 Va_max_sink=3.2085192552113512e+02 Bmax_sink=3.2085192552113513e-02
cycle=937053 time=5.0000000000000000e+01 dt=6.1268428794392094e-05
Terminating on time limit
time=5.0000000000000000e+01 cycle=937053
tlim=5.0000000000000000e+01 nlim=-1
Number of MeshBlocks = 736; 896  created, 224 destroyed during this simulation.
zone-cycles = 1779385974784
cpu time used  = 1.2288628745110000e+06
zone-cycles/cpu_second = 1.4479939232374232e+06
end_time   = 2026-09-12T09:23:54+09:00
exit_code  = 0

## 7.結果・考察等
Gは破綻しなかったが、速度向上が小さい一方、終盤の非線形イベントで積分量にも3%程度の差が生じたため、安全側として無制限のEを採用。
つまり、シンク内のvAに関する人工処理は一切用いないこととする。

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test86

## 1. 目的
fM=0.4(もとのガス質量の40%)、Bz0=100μGとして磁場の影響を確かめる比較テスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: f9614ca
変更点: fM=0.4、Bz0=100μG

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
Toyouchi-test86.sh

## 6.計算ログ
cycle=140696 time=1.9999907144150995e+01 dt=9.2855849004536140e-05
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=2.2521640559484808e-02 Mdot_out=1.5386472609517184e-04 Mdot_net=2.2367775833389638e-02 Msun/yr Mstar_diag=5.3606263849992115e+02 Msun (gravity fixed)
cycle=140697 time=2.0000000000000000e+01 dt=1.8571169800907228e-04
Terminating on time limit
time=2.0000000000000000e+01 cycle=140697
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 38755832832
cpu time used  = 6.8051460275000005e+04
zone-cycles/cpu_second = 5.6950773246283573e+05
end_time   = 2026-11-10T10:29:21+09:00
exit_code  = 0

## 7.結果・考察等
磁場はtest85の10μGよりはわずかに影響していると考えられるが、やはりほとんど影響していないようである。
現時点での考察として、初期磁場強度の大幅な増幅よりも初期ガス分布の20％程度の違いのほうがより大きな影響を与えているようである。
その考察の裏付けのためにfM=1のモデルに同じく10、100μGの磁場を入れ、磁場なしの場合とどれだけ差が出るかテストする予定。

## 8.備考
VTKファイルはすべて削除

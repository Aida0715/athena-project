# Toyouchi-test34

## 1. 目的
Toyouchi+23の設定で回転させ、１Myr計算。磁場BzはOFF

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: cd2fe2c
変更点: 計算時間を1Myr

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
Toyouchi-test34.sh

## 6.計算ログ
cycle=94714 time=6.0284774499848119e+01 dt=8.8097934021867432e-11
Sink: Mstar = 1.4778390399632258e+02, Mdot = 1.3575730391964713e+01, dM = 1.1959938003699678e-09
cycle=94715 time=6.0284774499936219e+01 dt=8.8259618443803518e-11
Sink: Mstar = 1.4778390399752075e+02, Mdot = 1.3575617984695606e+01, dM = 1.1981788634680711e-09
cycle=94716 time=6.0284774500024476e+01 dt=8.8421273517384808e-11
^Cend_time   = 2026-08-23T21:10:41+09:00
exit_code  = 1
ERROR: Athena++/MPI terminated abnormally. See /work/beta/aida/results/Toyouchi-test34/Toyouchi-test34.20260811-154747.log

## 7.結果・考察等
タイムステップがかなり小さくなって、計算が全く終わる見込みがなさそうなので、途中で止めた。
タイムステップが持っていかれる原因は、AMRだと考えられるので、その点明らかにして今後改善が必要と思われる。
本データは解析対象外。記録だけ残す。

## 8.備考
VTKファイルはすべて削除

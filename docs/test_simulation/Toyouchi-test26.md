# Toyouchi-test26

## 1. 目的
NFW、回転(vr,vΦ)、自己重力、中心星重力すべてONにして計算。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 9dd8f2a
変更点: 回転・重力等すべてONに戻した

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlle
  #Magnetic fields:              OFF
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
  #MPI parallelism:              OFF
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test26.sh

## 6.計算ログ
cycle=487 time=1.6710835295068821e+02 dt=0.0000000000000000e+00
cycle=488 time=1.6710835295068821e+02 dt=0.0000000000000000e+00
client_loop: send disconnect: Broken pipe

## 7.結果・考察等
同じ計算をローカルとnovaで３回回したが、毎回同じところで計算が落ちる。
物理量そのものが壊れたのではなく、どうやらAMR最細ブロックの境界の１セルで、out13→ out14間で密度が~1/5になる。
モーメンタムρv自体は保存していて、そのセルのvz=ρvz/ρ が密度低下が原因で暴走する(vx,vyは正常)。
原因解明のためにソースコードにデバッグコードを仕込み、密度がどのタイミングで落ちるか追跡する。
また変数補間も原因かもしれないとのことなので入力ファイルで、xorder=3→ 3cに変更してみる。

## 8.備考
VTKファイルはすべて削除

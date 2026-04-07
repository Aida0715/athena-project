# Toyouchi-test15

## 1. 目的
回転をVφ=9.47=constで入れていたが、これはToyouchi+23 fig.2の幾何平均で求めた代表値。
Toyouchi+23に倣い回転を、Vφ(r)=0.5(GM_enc/r)で入れ直し、挙動を確認。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 44fc89d
変更点:Vφ=9.47での計算を取り消しToyouchi+23に合わせた形で入れた

## 4. ビルド設定（configure）
　#Problem generator:          jeans
  #Coordinate system:          cartesian
  #Equation of state:          isothermal
  #Riemann solver:             hlle
  #Magnetic fields:            OFF
  #Relativistic dynamics:      OFF 
  #General relativity:         OFF 
  #Radiative Transfer:         OFF
  #Implicit Radiation:         OFF
  #Cosmic Ray Transport:       OFF
  #Cosmic Ray Diffusion:       OFF
  #Self-Gravity:               Multigrid
  #Super-Time-Stepping:        OFF
  #Floating-point precision:   double
  #Number of ghost cells:      4
  #MPI parallelism:            OFF
  #OpenMP parallelism:         OFF
  #FFT:                        OFF
  #HDF5 output:                OFF
  #Compiler:                   g++
  #Compilation command:        g++ -O3 -std=c++11

## 5. 対応run
Toyouchi-test15.sh

## 6.ログ記録
Terminating on time limit
time=2.0000000000000000e+01 cycle=78
tlim=2.0000000000000000e+01 nlim=50000
zone-cycles = 1308622848
cpu time used  = 2.4830464609999999e+03
zone-cycles/cpu_second = 5.2702310188468115e+05

## 7.結果・考察等
密度プロファイルによると、massが内側から外側に流れるような挙動が見られた。
これは、重力(中心星＋ガス自重) < 遠心力＋圧力勾配力＋数値粘性の影響 によって流れができたと思われる。
AMRを完全OFFにしたことによりグリッドが粗くなり、数値粘性の影響が大きく出ている可能性がある。
AMRをONにして再計算してみる。

## 8.備考
VTKファイルはすべて削除

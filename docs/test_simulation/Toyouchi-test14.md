# Toyouchi-test14

## 1. 目的
Toyouchi-test12でコアダンプの原因がAMRにあることが判明。
本計算ではさらにその原因がjeansリファインと密度リファインのいずれにあるのか、密度リファインのみ有効にしてテスト。
（threshold = 0.3）

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 9038d69
変更点: 密度リファインのみ有効

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
Toyouchi-test14.sh

## 6.ログ記録
cycle=20 time=2.0441869827877754e+00 dt=1.0222737151037083e-01
cycle=21 time=2.1464143542981464e+00 dt=1.0222976329041975e-01
cycle=22 time=2.2486441175885661e+00 dt=1.0223164030361968e-01
cycle=23 time=2.3508757578921857e+00 dt=1.0223379065444879e-01
cycle=24 time=2.4531095485466343e+00 dt=1.0223627386110067e-01
cycle=25 time=2.5553458224077348e+00 dt=1.0223909922742208e-01
cycle=26 time=2.6575849216351570e+00 dt=1.0224227327978090e-01
cycle=27 time=2.7598271949149380e+00 dt=1.0220866899035275e-01
run/Toyouchi-test14.sh: 14 行: 996297 Segmentation fault      (コアダンプ) ../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_6

## 7.結果・考察等
● 挙動
　・初期から中心領域が level 1 にリファイン
　・その後、レベル分布はほぼ固定（derefineなし）
　・dt は安定（CFL破綻なし）
● 異常
  ・最小密度：初期は一定（~3.2×10⁻⁶）、途中（~step 9）から急落
  ・最終的にコアダンプ
● 考えられる原因
  ・初期のrefine自体は正常
　・時間発展により中心密度が数値拡散で減少
　・本来はderefineすべき状況になるが、メッシュは維持される
　　→ 不要なcoarse–fine境界が残存
● 考えられる破綻メカニズム
  ・境界でのフラックス不整合により質量が徐々に流出
  ・誤差が蓄積し、ある時点で非線形的に崩壊
● 考えられる本質
問題はrefineではなく、不適切な位置に(今回は初期にレベル１にrefineされた内側)残り続けるAMR境界
　・refine/derefine条件の非対称性（ヒステリシスなし）が原因
　・物理構造に追従しないメッシュが不安定性を生む
● 改善案
　・refine/derefineにヒステリシスを導入し、AMR境界を物理構造に追従させる
## 8.備考
VTKファイルはすべて削除

# Toyouchi-test11

## 1. 目的
Toyouchi-test10で初期条件の密度フロアを削除して計算したが、VΦ爆発とは別の原因でコアダンプした。
恐らくメモリアクセスの問題で、本テストではAMRの境界条件を見直して再テスト。
＜考えられる原因＞
・RefinementCondition 内で中央差分
  w(i+1) - w(i-1)を使用
・ループが境界まで回っていたため
  → i-1, i+1 が配列外アクセス
・AMRによりブロック分割が変化し、
  これまで表面化していなかったバグが顕在化

＜本計算で講じた対策＞
・ループ範囲を内側に制限
for (k = ks+1; k <= ke-1)
for (j = js+1; j <= je-1)
for (i = is+1; i <= ie-1)

＜今後の対策＞
・AMR使用時は境界セルの扱いが非常に重要
・中央差分を使う場合は必ず
  境界除外 or ガードセル前提を確認
・AMRが「今まで動いていた」は安全性の保証にならない
　→ 今回のようにダイナミクス変化によるAMR分割変化によって計算破綻する可能性ある

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 9038d69
変更点: AMR境界扱い見直し

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
Toyouchi-test11.sh

## 6.ログ記録
cycle=0 time=0.0000000000000000e+00 dt=1.0220249351489546e-01
cycle=1 time=1.0220249351489546e-01 dt=1.0220178102960427e-01
cycle=2 time=2.0440427454449972e-01 dt=1.0220155216486063e-01
cycle=3 time=3.0660582670936032e-01 dt=1.0220179805085781e-01
cycle=4 time=4.0880762476021815e-01 dt=1.0220249665130686e-01
cycle=5 time=5.1101012141152502e-01 dt=1.0220348433510519e-01
cycle=6 time=6.1321360574663020e-01 dt=1.0220478458648924e-01
cycle=7 time=7.1541839033311949e-01 dt=1.0220485992841485e-01
cycle=8 time=8.1762325026153437e-01 dt=1.0220512688011167e-01
cycle=9 time=9.1982837714164600e-01 dt=1.0220574268916251e-01
cycle=10 time=1.0220341198308085e+00 dt=1.0220673550587264e-01
cycle=11 time=1.1242408553366812e+00 dt=1.0220809567505916e-01
cycle=12 time=1.2264489510117402e+00 dt=1.0221009268676511e-01
cycle=13 time=1.3286590436985053e+00 dt=1.0221237693192135e-01
cycle=14 time=1.4308714206304267e+00 dt=1.0221450179396260e-01
cycle=15 time=1.5330859224243893e+00 dt=1.0221604423088686e-01
cycle=16 time=1.6353019666552762e+00 dt=1.0221809345953246e-01
cycle=17 time=1.7375200601148086e+00 dt=1.0222015638347763e-01
cycle=18 time=1.8397402164982863e+00 dt=1.0222224312309557e-01
cycle=19 time=1.9419624596213818e+00 dt=1.0222452316639366e-01
cycle=20 time=2.0441869827877754e+00 dt=1.0222737151037083e-01
cycle=21 time=2.1464143542981464e+00 dt=1.0222976329041975e-01
cycle=22 time=2.2486441175885661e+00 dt=1.0223164030361968e-01
cycle=23 time=2.3508757578921857e+00 dt=1.0223379065444879e-01
cycle=24 time=2.4531095485466343e+00 dt=1.0223627386110067e-01
cycle=25 time=2.5553458224077348e+00 dt=1.0223909922742208e-01
cycle=26 time=2.6575849216351570e+00 dt=1.0224227327978090e-01
cycle=27 time=2.7598271949149380e+00 dt=1.0220866899035275e-01
run/Toyouchi-test11.sh: 14 行: 957397 Segmentation fault      (コアダンプ) ../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_3

## 7.結果・考察等
AMR境界付近の不適切な参照により、非物理的な状態が生成されて計算がSegmentation faultで破綻したと考えられる
● 症状
・一定時間後にSegmentation fault
・dtは安定（CFL起因ではない）
・全質量がステップ的に増加
・密度プロファイルは一見まとも（物理破綻っぽくない）
● 考えられる原因
・RefinementCondition 内で、AMR境界（ゴースト/隣接レベル）を跨いだ参照または 未定義セルの値参照
・特にこれ：w(i+1), w(i-1)
・AMRでは隣が「同じレベル」とは限らない
　→ 値が未定義 or 不整合
・その結果：refine/derefine判定が壊れる
　→ メッシュ構造崩壊
　→ segfault
● 修正
今は原因が物理にあるのか、AMR実装にあるのかを切り分ける必要がある。
まずAMRを完全にOFFにして、同じ計算を回してみる。
落ちる →　物理or初期条件に原因あり
落ちない → AMRに原因あり

## 8.備考
VTKファイルはすべて削除

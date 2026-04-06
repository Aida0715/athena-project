# Toyouchi-test12

## 1. 目的
Toyouchi-test11でAMRの境界条件を見直して再テストしたが、コアダンプ。
原因がAMRなのかどうかを見極めるため、本テストでは一旦AMRを完全OFFにして、VΦ回転のみでテスト。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 9038d69
変更点: AMR完全OFF

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
Toyouchi-test12.sh

## 6.ログ記録
Terminating on time limit
time=2.0000000000000000e+01 cycle=109
tlim=2.0000000000000000e+01 nlim=50000
zone-cycles = 1828716544
cpu time used  = 3.4041335859999999e+03
zone-cycles/cpu_second = 5.3720469476311561e+05

## 7.結果・考察等
AMR完全OFFにしたところ、問題なく計算は回った。
つまりこれまでの計算破綻は、AMRが原因であることが判明。
密度プロファイル、質量保存、密度map(グリッド付き)、流線、速度プロファイルを作成。
●　密度プロファイル
遠心力の影響で内側から外側に向かってmassが流れる様子が見られ、物理的に正常な結果と思われる。
●　質量保存
1.195→ 1.201で、~0.5%upしている。
これまで境界で(恐らく)数値拡散による質量の流入が確認されており、今回も同様に流入が見られた。
ただ回転の遠心力の影響か、0.5%の増加に抑えられた。
（回転なしのときは、~1.3%の増加）
●　密度map
グリッドはレベル０で描画されており、AMRはちゃんとOFFになっていたことがわかった。
こちらでも回転によってmassが外側に流れている様子が見られた。

●　速度プロファイル
時間発展させると内側でVφが落ちる挙動が見られた。
本計算の.cppは、Vφ~9.47の一定値で入れてしまっていた。
（これはToyouchi+23のfig.2を幾何平均したときのある半径での「代表値」）
Toyouchi+23の通り、9.47を初期値としてVφ(R)=0.5(GM_enc/R)で入れ直すべき。
→　入れ直して再計算

## 8.備考
VTKファイルはすべて削除

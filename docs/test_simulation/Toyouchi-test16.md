# Toyouchi-test16

## 1. 目的
Toyouchi-test14でAMR作動によるコアダンプの原因が、jeans refine=trueによることは判明。
原因として、refine_thr=0.3とderefine_thr=0.15(derefine_thrは.cppで指定していた)の差が小さすぎたことがコアダンプの原因だった。
→　一度refineされるとderefineされず、不要なAMR境界が残り続ける
そのため入力ファイルにパラメータとして新たにderefine_thrを追加し、derefine_thr=0.1と差を少し広げた。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 7cef4da
変更点:refine/derefineにヒステリシスを導入

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
Toyouchi-test16.sh

## 6.ログ記録
Terminating on time limit
time=2.0000000000000000e+01 cycle=144
tlim=2.0000000000000000e+01 nlim=50000
Number of MeshBlocks = 568; 56  created, 0 destroyed during this simulation.
zone-cycles = 2680160256
cpu time used  = 6.4258129650000001e+03
zone-cycles/cpu_second = 4.1709278975255700e+05

## 7.結果・考察等
計算は問題なく回った。
密度プロファイルによると、ガスは内側から外側に拡散しているようである。
またVφプロファイルを見ると、r<30で減速、r>30でわずかに加速している。
これはガスの拡散によって包含質量が減少することが主要因と考えられる。
また、r>30でのわずかな加速は（恐らく数値粘性による）角運動量輸送が生じていることが予測される。
合力（中心星重力、圧力勾配力、遠心力）プロファイルを見ると、圧力勾配が支配的で初期から常に外向きの力が働いているようである。
● test14までのAMR暴走の原因
　・refine / derefine が頻繁に往復（チャタリング）
  ・境界付近で数値ノイズに反応して不安定化
● 今回の解決策
　・ヒステリシスを導入し、refine条件 と derefine条件 を分離
　　→ refine条件(厳し目) と derefine条件(緩め)
● 結果
　・不要な再分割が止まる
　・AMR構造が安定
　・計算がクラッシュせず最後まで回る

## 8.備考
VTKファイルはすべて削除

# Toyouchi-test7

## 1. 目的
Toyouchi+23のVΦ=9.47を導入して挙動テスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 23d5653
変更点: VΦ=9.25を導入

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
Toyouchi-test7.sh

## 6. ログ記録
cycle=115 time=1.6861358481112261e+00 dt=1.3912953733732096e-02
cycle=116 time=1.7000488018449582e+00 dt=4.4434289842591793e-10
cycle=117 time=1.7000488022893012e+00 dt=8.3285157904892402e-13
cycle=118 time=1.7000488022901341e+00 dt=2.7311307370374017e-20
cycle=119 time=1.7000488022901341e+00 dt=2.8565092915040307e-26
cycle=120 time=1.7000488022901341e+00 dt=9.4020715739567691e-34
cycle=121 time=1.7000488022901341e+00 dt=9.6739801369226748e-40
cycle=122 time=1.7000488022901341e+00 dt=3.1961388494340570e-47
cycle=123 time=1.7000488022901341e+00 dt=3.2381764214447533e-53
cycle=124 time=1.7000488022901341e+00 dt=1.0739199842816158e-60
cycle=125 time=1.7000488022901341e+00 dt=1.0723511636906763e-66
cycle=126 time=1.7000488022901341e+00 dt=3.5700703269825657e-74
cycle=127 time=1.7000488022901341e+00 dt=3.5165887686678250e-80
cycle=128 time=1.7000488022901341e+00 dt=1.1752923564013332e-87
cycle=129 time=1.7000488022901341e+00 dt=1.1430166861785049e-93
cycle=130 time=1.7000488022901341e+00 dt=3.8350994112122134e-101
cycle=131 time=1.7000488022901341e+00 dt=3.6857010216003210e-107
cycle=132 time=1.7000488022901341e+00 dt=1.2415319779321408e-114
cycle=133 time=1.7000488022901341e+00 dt=1.1800656811648164e-120
cycle=134 time=1.7000488022901341e+00 dt=3.9908980249425581e-128
cycle=135 time=1.7000488022901341e+00 dt=3.7547828388032215e-134
cycle=136 time=1.7000488022901341e+00 dt=1.2749339118647286e-141
cycle=137 time=1.7000488022901341e+00 dt=1.1882893845738378e-147
cycle=138 time=1.7000488022901341e+00 dt=4.0511103986012028e-155
cycle=139 time=1.7000488022901341e+00 dt=3.7434913450493074e-161
cycle=140 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=141 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=142 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=143 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=144 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=145 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=146 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=147 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=148 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=149 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=150 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=151 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=152 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=153 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=154 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=155 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
cycle=156 time=1.7000488022901341e+00 dt=0.0000000000000000e+00
run/Toyouchi-test7.sh: 14 行: 892514 Segmentation fault      (コアダンプ) ../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_2

## 7.結果・考察等
回転導入時のCFL崩壊
● 症状
  ・ 計算初期は正常に進行
  ・後半で dt が指数的に減少
  ・最終的に dt → 0 → segmentation fault
　・密度プロファイルと密度mapを作成したら、やはり早い段階で中心から急激に低密度化していた
● 原因
円筒回転vϕ=constを、
vx = -vphi * y / r_cyl_safe;
vy =  vphi * x / r_cyl_safe;　として実装した際、
r_cyl_safe = max(r_cyl, 1e-12);　としていたため、中心付近（r → 0）でV~1/rが発散し速度が極端に大きくなった。
その結果、dtが極端に小さくなり最終的に数値アンダーフローを起こした。
● 修正
r_cyl_safe = max(r_cyl, dx);
● 教訓
  ・1/r 型の速度場は必ず中心でカットオフが必要
  ・カットオフスケールは物理量ではなく セルサイズ(dx) に合わせる
  ・AMR環境では dx が最小セルにより決まるため特に注意


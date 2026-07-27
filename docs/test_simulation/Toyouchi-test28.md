# Toyouchi-test28

## 1. 目的
MHD計算できるようにした。Bz=1μGを追加して安定性テスト。AMRと回転はOFFにした。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 3cc79ae
変更点: MHD計算できるようにした

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
  #MPI parallelism:              OFF
  #OpenMP parallelism:           OFF
  #FFT:                          OFF
  #HDF5 output:                  OFF
  #Compiler:                     g++
  #Compilation command:          g++  -O3 -std=c++11


## 5. 対応run
Toyouchi-test28.sh

## 6.計算ログ
cycle=15 time=1.0000000000000000e+03 dt=3.3844637485867999e+01
Terminating on time limit
time=1.0000000000000000e+03 cycle=15
tlim=1.0000000000000000e+03 nlim=50000
Number of MeshBlocks = 8; 0  created, 0 destroyed during this simulation.
zone-cycles = 3932160
cpu time used  = 5.4328709999999996e+00
zone-cycles/cpu_second = 7.2377201667405688e+05

## 7.結果・考察等
・３磁場成分の可視化
Bz~1μGでほぼ一定、Bx~Byで発展しており球対称分布が保たれたまま収縮していると考えられる
・プラズマβ分布の可視化
β~32がピーク。分布全体がβ=1より大きく、ガス優勢であることがわかる。
ただしこれはβのセル数分布であり、AMRなどをONにすると変わる可能性がある。物理的には質量加重分布が良いらしい。
・アルフベン・マッハ数
分布のピークが~6なので、流体運動は超アルフベン的であり、磁場が流れを制御している状態ではなく磁場に引きずられて変形して増幅していると思われる。

弱い磁場の理想MHDのテスト計算としてはとりあえず問題ないと思う。
Test29でBz=3μGに磁場を強めて同様の計算をしてみる。

## 8.備考
VTKファイルはすべて削除

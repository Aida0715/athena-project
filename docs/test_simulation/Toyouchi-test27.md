# Toyouchi-test27

## 1. 目的
cycle=71あたりで毎回密度が負になり、Vzが爆発→dt~0になり計算が落ちる。
デバッグコードを複数仕込み、原因を追求するためテスト計算を回す。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 48a42b7
変更点: 密度がいつ落ちるか追跡するデバッグコードを追記。入力ファイルでxorderを３から3cへ変更

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
Toyouchi-test27.sh

## 6.計算ログ
Terminating on cycle limit
time=8.0952256367946541e+02 cycle=50000
tlim=1.0000000000000000e+03 nlim=50000
Number of MeshBlocks = 764; 1008  created, 252 destroyed during this simulation.
zone-cycles = 1323501322240
cpu time used  = 9.4496118533899996e+05
zone-cycles/cpu_second = 1.4005880270788062e+06
run/Toyouchi-test27.sh: 16 行: 1119610 Segmentation fault      (コアダンプ) "$ATHENA/bin/athena" -i "$ATHENA/inputs/hydro/Toyouchi_test/athinput.Toyouchi_19"
aida@nova:~/athena-project$ client_loop: send disconnect: Broken pipe

## 7.結果・考察等
コアダンプになっているが、物理的な破綻や計算破綻ではなさそう。
cycle limitに達して設定通り計算は終了し、少なくともそこまででコアダンプになったわけではない。
終了後に終了処理・デストラクタ固有の不具合or計算中の配列外アクセスなどでメモリが壊れ、終了時のdeleteで初めて表面化したが原因として考えられるらしい。
ｖΦは正しく与えられていたが、ｖrが誤っていた。本来中心に向かう弱い降着流のはずが円筒半径を使ってしまったため、vrはｚ軸に向かう垂直な流れになってしまった。ここは修正が必要である（test32以降で反映）。
解析の結果、期待通り円盤構造ができていた。

## 8.備考
VTKファイルはすべて削除

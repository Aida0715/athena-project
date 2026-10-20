# Toyouchi-test76

## 1. 目的
シンクは廃止し、20kyrの軽いテスト。
10^3auより内側のMenc、VrプロファイルはWise+19のfig.4を７次多項式フィッティングして取得。
音速は10^3K→10^4Kの場合に上げて変更。
ただしVrは570au、Mencは140auより内側のプロファイルは得られなかったため、それぞれより内側は中心で０になるよう滑らかに接続するよう実装。
SMRレベルは５に戻した。
磁場はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a04dffd
変更点: 上記の通り

## 4. ビルド設定（configure）
　#Problem generator:            Toyouchi
  #Coordinate system:            cartesian
  #Equation of state:            isothermal
  #Riemann solver:               hlld
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
  #Self-Gravity:                 Multigrid
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
Toyouchi-test76.sh

## 6.計算ログ
cycle=874 time=1.9959660459214025e+01 dt=2.5359842175130941e-02
cycle=875 time=1.9985020301389156e+01 dt=1.4979698610844139e-02
cycle=876 time=2.0000000000000000e+01 dt=2.5362579315566919e-02
Terminating on time limit
time=2.0000000000000000e+01 cycle=876
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 241299456
cpu time used  = 4.7852287000000001e+02
zone-cycles/cpu_second = 5.0425898348390328e+05
end_time   = 2026-10-18T19:43:44+09:00
exit_code  = 0

## 7.結果・考察等
密度mapと密度プロファイルを可視化したところ、ガスがアクリージョンせず外側に広がった。
温度を10^3K→10^4Kに上げると、圧力が19倍(μの変更も考慮すると)になり、強い圧力勾配力でガスが外側に押しのけられたものと思われる。
また赤道面でのガス自己重力と遠心力を比較すると、遠心力は自己重力のおよそ1/4。
このことからも強すぎる圧力購買力によって今回のような結果になったと考えられる。
温度のみ10^3Kに戻してテストし、比較してみる。

## 8.備考
VTKファイルはすべて削除

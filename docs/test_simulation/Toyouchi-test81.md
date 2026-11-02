# Toyouchi-test81

## 1. 目的
音速を10^3Kに戻して定常円盤ができるかテスト

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: a4554e8
変更点: 温度を10^3K

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
Toyouchi-test81.sh

## 6.計算ログ
cycle=598534 time=1.9999508735663106e+01 dt=4.9126433689394844e-04
VirtualSphere: time=2.0000000000000000e+04 yr Mdot_in=2.0936430776534581e-01 Mdot_out=1.3183214236160926e-01 Mdot_net=7.7532165403736572e-02 Msun/yr Mstar_diag=2.0449640214488034e+03 Msun (gravity fixed)
cycle=598535 time=2.0000000000000000e+01 dt=9.8252867378789688e-04
Terminating on time limit
time=2.0000000000000000e+01 cycle=598535
tlim=2.0000000000000000e+01 nlim=-1
zone-cycles = 164870056960
cpu time used  = 3.2733051896500000e+05
zone-cycles/cpu_second = 5.0368067567090748e+05
end_time   = 2026-11-01T17:11:32+09:00
exit_code  = 0

## 7.結果・考察等
面密度map、面密度プロファイル等を見ると前回までのようにクランプに分裂したりせず、おおよそ定常円盤が形成されたと考えられる。
ただし分裂はしないがまだ少し不安定そうな要素があるので（スパイラル構造が安定しないなど）、ガスの量を減らしてより安定なモデルを考える。
ただし単に安定させるのが重要ではなく、初期条件による結果の違いと原因となる物理を解析することが必要。
Wise+19の分布を保ったままファクターを変更してガスの量を調節し、比較する予定。
また、中心星質量を上げたせいかタイムステップが持っていかれ、前回のモデルと比較して計算にかなりの時間を要した。
特にr<10^3auで最低密度が実際の平均密度より桁で小さい1e-7g/cm^3付近のセルが存在し、そこで局所的に速度が上がりタイムステップを持っていっているようである。
r<10^3auより内側の密度フロアを底上げするなどして対策を考えたかったが、その範囲でρ<1e-7g/cm^3（底上げ候補）であるガスの割合が十数kyrあたりから30％以上
になっていたため、局所的な密度のオーバーシュートではないと判断し、底上げはやめる方針で。

## 8.備考
VTKファイルはすべて削除

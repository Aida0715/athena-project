# Toyouchi-test75

## 1. 目的
SMRレベルを５から４へ下げて、クランプの出現する時刻及び分裂が解像度に依存するかテスト。
test74などと比較すれば良い。
磁場はOFF。

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: c79c93d
変更点: SMRレベルを５から４へ

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
Toyouchi-test75.sh

## 6.計算ログ
cycle=69193 time=5.0172893149408907e+01 dt=3.1052831964927088e-04
Sink(drain): Mstar=2.6888659600575136e+02 dM_star=5.0209342166087332e-04 dM_flux=4.9172023790634237e-04 Msink_gas=4.1993804113553096e-01 tau=2.1705922113391542e-02
cycle=69194 time=5.0173203677728559e+01 dt=3.4102624439924512e-06
Sink(drain): Mstar=2.6888660150909470e+02 dM_star=5.5033433235090753e-06 dM_flux=5.4013651400693235e-06 Msink_gas=4.1993793915734756e-01 tau=2.1705901847588124e-02
cycle=69195 time=5.0173207087991003e+01 dt=6.8205248879849023e-06
Terminating on time limit
time=5.0173207087991003e+01 cycle=69195
tlim=5.0173207087991003e+01 nlim=-1
zone-cycles = 17076218880
cpu time used  = 2.7913531084999999e+04
zone-cycles/cpu_second = 6.1175416424388939e+05
end_time   = 2026-10-16T00:00:52+09:00
exit_code  = 0

## 7.結果・考察等
test74と比較したところ、リングが形成される時刻や生成するクランプ数は変わらなかったが、その後の不安定の進化が解像度によって異なっていた。

## 8.備考
VTKファイルはすべて削除

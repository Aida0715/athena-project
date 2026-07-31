#!/bin/bash

# commit: bb069ed //本シミュレーションに対応するToyouchi.cppをgitの履歴から追跡可
# Takasao+22の計算時間tlim=26.3Myrに合わせた。(AMRと回転をOFF)Bz=3μGの磁場はTest29と同じ。

ATHENA=$HOME/athena-project
WORK=/work/beta/aida

OUTDIR=$WORK/results/Toyouchi-test30
mkdir -p $OUTDIR
cd $OUTDIR

# 実行（ファイル名はinputに任せる）
# 環境作成中のテスト計算用入力ファイルはToyouchi_testディレクトリに格納
# 並列計算のコア数は-np32のように指定する
mpirun -np 32 "$ATHENA/bin/athena" \
    -i "$ATHENA/inputs/hydro/Toyouchi_test/athinput.Toyouchi_22"

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



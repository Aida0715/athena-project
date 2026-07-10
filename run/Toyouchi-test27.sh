#!/bin/bash

# commit: 48a42b7  //本シミュレーションに対応するToyouchi.cppをgitの履歴から追跡可
# 密度がいつ落ちるか追跡するデバッグコードを追記。入力ファイルでxorderを３から3cへ変更

ATHENA=$HOME/athena-project
WORK=/work/beta/aida

OUTDIR=$WORK/results/Toyouchi-test27
mkdir -p $OUTDIR
cd $OUTDIR

# 実行（ファイル名はinputに任せる）
# 環境作成中のテスト計算用入力ファイルはToyouchi_testディレクトリに格納
"$ATHENA/bin/athena" \
    -i "$ATHENA/inputs/hydro/Toyouchi_test/athinput.Toyouchi_19"

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



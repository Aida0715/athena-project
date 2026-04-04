#!/bin/bash

# commit: 9038d69  //本シミュレーションに対応するjeans.cppをgitの履歴から追跡可
# test11でもコアダンプしたので、原因切り分けのためAMRを完全OFFで再計算

cd "$(dirname "$0")/.."

OUTDIR=results/Toyouchi-test12
mkdir -p $OUTDIR
cd $OUTDIR

# 実行（ファイル名はinputに任せる）
# 環境作成中のテスト計算用入力ファイルはToyouchi_testディレクトリに格納
../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_4

# configure
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

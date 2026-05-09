#!/bin/bash

# commit: a2f95b5  //本シミュレーションに対応するjeans.cppをgitの履歴から追跡可
# 領域をToyouchi+23に合わせて、10^7auに拡張
# AMR、回転(vr,vphi)、中心星重力、自己重力は完全OFF

cd "$(dirname "$0")/.."

OUTDIR=results/Toyouchi-test19
mkdir -p $OUTDIR
cd $OUTDIR

# 実行（ファイル名はinputに任せる）
# 環境作成中のテスト計算用入力ファイルはToyouchi_testディレクトリに格納
../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_11

# configure
 # Problem generator:            Toyouchi
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
 #Self-Gravity:                 OFF
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


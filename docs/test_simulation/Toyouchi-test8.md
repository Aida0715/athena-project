# Toyouchi-test8

## 1. 目的
Toyouchi-test7でVΦの中心カットオフがまずく、早い段階で中心付近で回転速度の爆発が見られ、コアダンプした。
本テストでは中心カットオフを以下のように修正し、再テストを行った。
// --- cutoff radius ---
	Real r_cut = std::max(racc, 1e-12);

	// --- smooth cylindrical rotation ---
	if (vphi0 != 0.0) {

  	  Real vphi_eff = vphi0;

  	  // 中心で滑らかにゼロへ（剛体回転化させて安定化を図った）
  	  if (r_cyl < r_cut) {
    	    vphi_eff *= r_cyl / r_cut;
  	  }

  	  Real inv_r = 1.0 / std::max(r_cyl, 1e-12);

  	  vx += -vphi_eff * y * inv_r;
  	  vy +=  vphi_eff * x * inv_r;
	  }


## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 0294366
変更点: VΦの中心カットオフを修正

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
Toyouchi-test8.sh

## 6. ログ記録
cycle=113 time=1.6612575468192174e+00 dt=1.3431775747503106e-02
cycle=114 time=1.6746893225667205e+00 dt=8.3637062155877614e-03
cycle=115 time=1.6830530287823082e+00 dt=3.5402936508037066e-10
cycle=116 time=1.6830530291363377e+00 dt=4.5328957857247227e-16
cycle=117 time=1.6830530291363381e+00 dt=1.8318887390420209e-22
cycle=118 time=1.6830530291363381e+00 dt=3.7660723326276558e-30
cycle=119 time=1.6830530291363381e+00 dt=1.5296681889845481e-36
cycle=120 time=1.6830530291363381e+00 dt=3.1449521023388872e-44
cycle=121 time=1.6830530291363381e+00 dt=1.2838313100394540e-50
cycle=122 time=1.6830530291363381e+00 dt=2.6396835146081526e-58
cycle=123 time=1.6830530291363381e+00 dt=1.0830070056639907e-64
cycle=124 time=1.6830530291363381e+00 dt=2.2269072084813023e-72
cycle=125 time=1.6830530291363381e+00 dt=9.1826310147653489e-79
cycle=126 time=1.6830530291363381e+00 dt=1.8882731095052993e-86
cycle=127 time=1.6830530291363381e+00 dt=7.8255623313212806e-93
cycle=128 time=1.6830530291363381e+00 dt=1.6093107511736226e-100
cycle=129 time=1.6830530291363381e+00 dt=6.7031126835412753e-107
cycle=130 time=1.6830530291363381e+00 dt=1.3785655579110231e-114
cycle=131 time=1.6830530291363381e+00 dt=5.7709859738749330e-121
cycle=132 time=1.6830530291363381e+00 dt=1.1869360827666492e-128
cycle=133 time=1.6830530291363381e+00 dt=4.9938562902353573e-135
cycle=134 time=1.6830530291363381e+00 dt=1.0271636198321955e-142
cycle=135 time=1.6830530291363381e+00 dt=4.3434474826339645e-149
cycle=136 time=1.6830530291363381e+00 dt=8.9343776124619149e-157
cycle=137 time=1.6830530291363381e+00 dt=3.7970438425563723e-163
cycle=138 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=139 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=140 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=141 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=142 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=143 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=144 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=145 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=146 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=147 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=148 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=149 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=150 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=151 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=152 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=153 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
cycle=154 time=1.6830530291363381e+00 dt=0.0000000000000000e+00
run/Toyouchi-test8.sh: 14 行: 909743 Segmentation fault      (コアダンプ) ../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_3

## 7.結果・考察等
回転導入時のCFL崩壊
● 症状
前回のToyouchi-test7の時と同様に、初期段階で中心が崩壊
● 原因
前回同様、VΦの中心カットオフがまずかったと思われる
● 修正
　・r<r_cutでのカットオフ
　　線形減衰→ ２次減衰でより平滑化
　・VΦの初期条件にケプラー速度で上限をつけ、初速の暴走を防止
　・密度フロアを設定し、真空領域を人工的に補完（今回の問題とは無関係だが安全装置として実装）
## 8.備考
VTKファイルはすべて削除

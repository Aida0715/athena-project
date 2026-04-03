# Toyouchi-test9

## 1. 目的
Toyouchi-test8でVΦの中心カットオフがまずく、早い段階で中心付近で回転速度の爆発が見られ、Test7同様コアダンプした。
本テストでは中心カットオフを以下のように修正し、再テストを行った。
 // --- 回転速度（解析的M_encを使用）---
        if (vphi0 != 0.0) {
            Real vphi_eff = vphi0;

            // 中心で滑らかにゼロへ（二次減衰で安定化）
            if (r_cyl < r_cut) {
                vphi_eff *= (r_cyl / r_cut) * (r_cyl / r_cut);
            }

            // ケプラー上限の計算（M_enc = 星質量 + ガス質量）
            // ρ = 0.11 * r^{-1.75} の解析的積分
            // M_gas(r) = 4π * 0.11 * r^{1.25} / 1.25
            Real M_gas = 4.0 * M_PI * 0.11 * pow(r_cyl, 1.25) / 1.25;
            Real M_enc = Mstar + M_gas;
            
            Real v_max_rot = std::sqrt(gconst * M_enc / std::max(r_cyl, r_cut));
            if (vphi_eff > v_max_rot && v_max_rot > 0) {
                vphi_eff = 0.9 * v_max_rot;
            }

            Real inv_r = 1.0 / std::max(r_cyl, 1e-12);
            vx += -vphi_eff * y * inv_r;
            vy +=  vphi_eff * x * inv_r;
        }

つまり、
・r<r_cutでのカットオフを２次減衰でより平滑化（test8では線形減衰だった）
・初期条件で速度にケプラー速度の上限をつけ、速度が暴走しないようにした
・今回の問題とは無関係だが、密度フロアを設定し真空発生防止の安全策とした

## 2. 参照
論文：
Toyouchi+2023:https: //arxiv.org/pdf/2206.14459

## 3. 使用コード状態
commit ID: 9571459
変更点: VΦの中心カットオフを再修正

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

## 6.ログ記録
cycle=0 time=0.0000000000000000e+00 dt=1.0220249351489546e-01
cycle=1 time=1.0220249351489546e-01 dt=1.0220178102960427e-01
cycle=2 time=2.0440427454449972e-01 dt=1.0220155216486063e-01
cycle=3 time=3.0660582670936032e-01 dt=1.0220179805085781e-01
cycle=4 time=4.0880762476021815e-01 dt=1.0220249665130686e-01
cycle=5 time=5.1101012141152502e-01 dt=1.0220348433510519e-01
cycle=6 time=6.1321360574663020e-01 dt=1.0220478458648924e-01
cycle=7 time=7.1541839033311949e-01 dt=1.0220485992841485e-01
cycle=8 time=8.1762325026153437e-01 dt=1.0220512688011167e-01
cycle=9 time=9.1982837714164600e-01 dt=1.0220574268916251e-01
cycle=10 time=1.0220341198308085e+00 dt=1.0220673550587264e-01
cycle=11 time=1.1242408553366812e+00 dt=1.0220809567505916e-01
cycle=12 time=1.2264489510117402e+00 dt=1.0221009268676511e-01
cycle=13 time=1.3286590436985053e+00 dt=1.0220555026880536e-01
run/Toyouchi-test9.sh: 14 行: 953507 Segmentation fault      (コアダンプ) ../../bin/athena -i ../../inputs/hydro/Toyouchi_test/athinput.Toyouchi_3

## 7.結果・考察等
回転導入時のCFL崩壊とは別の問題。
● 症状
前回のtest8よりも早い段階で計算破綻。
計算ログ、作成した密度プロファイル、密度mapなどを見る限り前回のような中心での回転速度爆発が原因ではないと考えられる。
● 原因
VΦの中心カットオフは問題なく実装できたと思われる。
密度フロアを初期条件にかけており、これがAMRと干渉し、途中で計算が破綻したと思われる。
● 修正
密度フロアの設定を削除し、計算を回してみる。
これで問題なく回れば、密度フロアは時間発展のコードに組み込む予定。

// 密度フロア（最小密度を設定）
Real rho_min = 1e-6 * rho_profile;
rho_final = std::max(rho_profile, rho_min);
P_final = rho_final * cs * cs;
→　これを削除

## 8.備考
VTKファイルはすべて削除

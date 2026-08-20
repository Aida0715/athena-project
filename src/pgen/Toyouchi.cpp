//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file Toyouchi.cpp
//! \brief Linear wave problem generator for 1D/2D/3D problems including self-gravity
//!
//! In 1D, the problem is setup along one of the three coordinate axes (specified by
//! setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
//! automatically sets the wavevector along the domain diagonal.
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <fstream>    // 追加
#include <vector>     // 追加
#include <utility>    // 追加

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

void CentralGravity(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc,
                    AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar);

int RefinementCondition(MeshBlock *pmb); // AMRのリファインメント宣言

Real SinkHistory(MeshBlock *pmb, int iout); // hstファイルへ中心星質量・降着率を出力

namespace {
// with functions A1,2,3 which compute vector potentials
Real cs2, gam, gm1, gconst;

// AMR関連のグローバル変数（RefinementCondition関数で使用）
bool use_jeans_refine = true;   // Jeans長ベースのリファインを使用するか
bool use_grad_refine = false;   // 密度勾配ベースのリファインを使用するか
Real jeans_cells = 8.0;         // Jeans長を何セルで解像するか
Real refine_thr = 0.3;          // 密度勾配の閾値（use_grad_refine=true時)
Real derefine_thr = 0.1;

// rotation parameters
Real vr0   = 0.0;   // radial velocity
bool use_rotation = true;  // 回転のON,OFFは入力ファイルで指定

// 中心星重力
bool use_central_gravity = true;  // 中心星重力のON,OFFは入力ファイルで指定

// NFW halo gravity
bool use_nfw_gravity = false;

// NFW halo parameters (CGS)
const Real Phi0_phys = 1.02e13;   // cm^2 s^-2
const Real Rs_phys   = 5.71e20;   // cm

//===== Toyouchi+23 parameters =====
Real epsilon_soft = 0.5;   // gravitational softening

// ===== unit system (code unit <-> cgs) =====
const Real Munit = 4.0e33;   // g (central star mass)
const Real Lunit = 6.7e15;   // cm (bondi radius)
const Real Tunit = 3.34e10;  // s (free fall time)

const Real Vunit = Lunit/Tunit;                 // velocity unit
const Real Rhounit = Munit/(Lunit*Lunit*Lunit); // density unit
const Real Punit = Rhounit*Vunit*Vunit;         // pressure unit

// unit conversion
const Real AU = 1.496e13;

// 中心カットオフ用（速度減衰）
Real racc = 0.0;   // accretion radius (cutoff scale)

// ===== シンク用変数の宣言 =====
bool use_sink = true;

// シンク半径
Real r_sink = 0.0;

// シンク内部に残す密度フロア
Real sink_rho_floor = 1.0e-8;

// sink内部で許容するAlfven速度の上限 [code velocity]
// <= 0 の場合はAlfven速度制限を使用しない
Real sink_va_cap = 10.0;

// シンク表面の外側で逆流を抑制する幅
// 各セルで sink_no_outflow_cells * dx として使用
Real sink_no_outflow_cells = 1.0;

// AMRで強制的に細分化するシンク外側の幅
Real sink_refine_buffer = 1.0;

// ===== sink diagnostic data =====
//
// ruser_mesh_data[0]に保存する量。
// restartファイルにも保存される。
//
enum SinkDataIndex {
  // 中心星の累積質量 [code mass]
  SINK_MSTAR = 0,

  // sink境界を内向きに通過した正味質量流入率
  // [code mass / code time]
  SINK_MDOT_FLUX = 1,

  // 現在の1 cycleでsink内部セルから除去した質量率
  // 中心星質量には加えない診断量
  // [code mass / code time]
  SINK_MDOT_RESET = 2,

  // Alfven速度floorにより人工的に追加した質量率
  // 中心星質量には加えない診断量
  // [code mass / code time]
  SINK_MDOT_FLOOR = 3,

  // sink内部に現在存在するガス質量
  // 人工的なatmosphere質量を含む [code mass]
  SINK_MGAS = 4,

  // 配列要素数
  NSINK_DATA = 5
};

// ===== magnetic field =====
Real bz_microgauss = 0.0;  // input value [microgauss]
Real bz_code = 0.0;        // Athena++へ実際に渡す無次元磁場(z方向)
Real Bunit = 0.0;          // magnetic-field unit [Gauss]

} // namespace


void Mesh::InitUserMeshData(ParameterInput *pin) {

  epsilon_soft = pin->GetOrAddReal("problem","epsilon",0.5);
  gconst = pin->GetOrAddReal("problem","grav_const",1.0);

  // ===== sink parameters =====
  use_sink =
      pin->GetOrAddBoolean("problem", "use_sink", true);

  // シンク半径の入力はAU単位、内部ではコード単位を使用
  Real r_sink_au =
      pin->GetOrAddReal("problem", "r_sink_au", 1000.0);
  r_sink = r_sink_au * AU / Lunit;
   
  // シンク領域内の密度フロア
  sink_rho_floor =
      pin->GetOrAddReal("problem", "sink_rho_floor", 1.0e-8);

  // sink内部のAlfven速度上限 [code velocity]
  // vA = |B|/sqrt(rho) <= sink_va_cap となるように数値的密度floorを調整する。
  // 0以下なら磁場依存のdensity floorを使用しない。
  sink_va_cap =
      pin->GetOrAddReal("problem", "sink_va_cap", 10.0);

  // シンク表面で外向き流を禁止するセル数
  sink_no_outflow_cells =
      pin->GetOrAddReal("problem", "sink_no_outflow_cells", 1.0);

  // シンク半径に対する外側バッファの比率
  Real sink_refine_buffer_factor =
      pin->GetOrAddReal("problem", "sink_refine_buffer_factor", 1.0);
  sink_refine_buffer = sink_refine_buffer_factor * r_sink;

  // 中心星の初期質量（２M_solar）
  Real mstar_init =
      pin->GetOrAddReal("problem", "mstar_init", 1.0);
  // ===== sink parameters読み込み終わり =====
      
  // restartファイルへ保存されるMeshデータ
  // Athena++のrestart読み込みでは、この初期値はrestartファイル中の値で上書きされる。したがって再スタート時にも成長後の中心星質量を引き継げる
  AllocateRealUserMeshDataField(1);
  ruser_mesh_data[0].NewAthenaArray(NSINK_DATA);

  // 中心星質量 [code mass]
  ruser_mesh_data[0](SINK_MSTAR) = mstar_init;

  // sink境界の正味質量流入率 [code mass / code time]
  ruser_mesh_data[0](SINK_MDOT_FLUX) = 0.0;

  // sink内部セルのリセットで除去した質量率
  // [code mass / code time]
  ruser_mesh_data[0](SINK_MDOT_RESET) = 0.0;

  // Alfven速度floorで人工的に追加した質量率
  // [code mass / code time]
  ruser_mesh_data[0](SINK_MDOT_FLOOR) = 0.0;

  // sink内部に残るガス質量 [code mass]
  ruser_mesh_data[0](SINK_MGAS) = 0.0;

  // 履歴出力を登録する
  // sink診断量をhstへ出力する
  AllocateUserHistoryOutput(5);

  // 中心星の累積質量 [code mass]
  EnrollUserHistoryOutput(
      0,
      SinkHistory,
      "Mstar"
  );

  // sink境界を通過した正味質量流入率
  // [code mass / code time]
  EnrollUserHistoryOutput(
      1,
      SinkHistory,
      "Mdot_flux"
  );

  // sink内部セルをrho_targetへ戻す際に除去した質量率
  // 中心星質量には加えない
  // [code mass / code time]
  EnrollUserHistoryOutput(
      2,
      SinkHistory,
      "Mdot_reset"
  );

  // Alfven速度floorを満たすため人工的に追加した質量率
  // [code mass / code time]
  EnrollUserHistoryOutput(
      3,
      SinkHistory,
      "Mdot_floor"
  );

  // sink内部に残っているガス質量
  // [code mass]
  EnrollUserHistoryOutput(
      4,
      SinkHistory,
      "Msink_gas"
  );

  // 初期設定をログに表示
  if (Globals::my_rank == 0) {
    std::cout
        << "Sink region:" << std::endl
        << "  enabled          = " << use_sink << std::endl
        << "  r_sink           = " << r_sink << " code units" << std::endl
        << "  r_sink           = " << r_sink_au << " AU" << std::endl
        << "  rho floor        = " << sink_rho_floor << std::endl
        << "  initial Mstar    = " << mstar_init << " code units" << std::endl
        << "  Alfven speed cap = " << sink_va_cap << " code velocity"
        << std::endl;
  }

  // 中心星重力のON,OFF読み込み
  use_central_gravity =
    pin->GetOrAddBoolean("problem","use_central_gravity",true);
  use_nfw_gravity =
    pin->GetOrAddBoolean("problem","use_nfw_gravity",false);
  if (use_central_gravity || use_nfw_gravity) {
      EnrollUserExplicitSourceFunction(CentralGravity);
  }

  // EOS関連の設定
  if (NON_BAROTROPIC_EOS) {
    gam = pin->GetReal("hydro","gamma");
    gm1 = gam-1.0;
    cs2 = gam * pin->GetReal("hydro", "gamma") / 1.0;  // p0=1.0/gamを使用
  } else {
    Real iso_cs = pin->GetReal("hydro","iso_sound_speed");
    cs2 = SQR(iso_cs);
  }

  //追加 回転パラメータ読み込み
  vr0   = pin->GetOrAddReal("problem","vr0", 0.0);
  use_rotation = pin->GetOrAddBoolean("problem","use_rotation",true);

  // raccを入力ファイルから
  racc = pin->GetOrAddReal("problem", "racc", 0.44);

  //磁場の強さを入力ファイルから読み込む(μG単位)
  bz_microgauss = pin->GetOrAddReal("problem", "bz_microgauss", 0.0);

  // CGS→コード単位の変換基準
  Bunit = std::sqrt(4.0*M_PI*Punit);

  // 磁場をCGS→コード単位へ変換
  Real bz_phys = bz_microgauss * 1.0e-6;  // μG→Gaussに変換
  bz_code = bz_phys / Bunit;

  // 磁場関連のメッセージ表示
  if (Globals::my_rank == 0) {
    std::cout
        << "Initial vertical magnetic field:" << std::endl
        << "  Bz physical = " << bz_microgauss << " microgauss" << std::endl
        << "  B unit      = " << Bunit << " Gauss" << std::endl
        << "  Bz code     = " << bz_code << std::endl;
  }

  //重力定数（入力ファイルから読み込む）
  if (SELF_GRAVITY_ENABLED) {
    SetGravitationalConstant(gconst);
  }

  // AMR関連の設定
  if (adaptive) {
    // AMRリファイン方式の選択（入力ファイルから読み込み）
    // デフォルト：Jeans長ベースのリファインを使用（自己重力系では必須）
    bool jeans_refine_input =
        pin->GetOrAddBoolean("problem", "use_jeans_refine", true);
    bool grad_refine_input =
        pin->GetOrAddBoolean("problem", "use_grad_refine", false);
    
    // 入力パラメータのチェックと警告
    // 両方OFFの場合はjeansリファインが作動する
    if (!jeans_refine_input && !grad_refine_input) {
      std::cout << "WARNING: AMR enabled but no refinement condition specified!" << std::endl;
      std::cout << "         Enabling Jeans length refinement as default." << std::endl;
      jeans_refine_input = true;
    }

    // namespace内の変数へ必ず反映する（falseの場合も古い初期値を残さない）
    use_jeans_refine = jeans_refine_input;
    use_grad_refine = grad_refine_input;
    
    // Jeans長ベースのリファインが有効な場合のパラメータ読み込み
    if (jeans_refine_input) {
      // Jeans長を何セルで解像するか（デフォルト: 8.0）
      Real jeans_cells = pin->GetOrAddReal("problem", "jeans_cells", 8.0);
      
      // 妥当性チェック（Jeans長は最低4セル以上で解像することが推奨）
      if (jeans_cells < 4.0) {
        std::cout << "WARNING: jeans_cells = " << jeans_cells << " is too small!" << std::endl;
        std::cout << "         Setting to 4.0 (minimum recommended value)" << std::endl;
        jeans_cells = 4.0;
      }
      
      // グローバル変数に保存（RefinementCondition関数で使用）
      // 注：これらの変数はnamespace内で定義されている必要があります
      ::jeans_cells = jeans_cells;
      
      std::cout << "AMR: Jeans length refinement enabled" << std::endl;
      std::cout << "     Cells per Jeans length = " << jeans_cells << std::endl;
    }
    
    // 密度勾配ベースのリファインが有効な場合のパラメータ読み込み
    if (grad_refine_input) {
      // 密度勾配の閾値（必須パラメータ）
      refine_thr = pin->GetReal("problem", "refine_thr");
      derefine_thr = pin->GetOrAddReal("problem", "derefine_thr", 0.1);
      
      // グローバル変数に保存
      ::refine_thr = refine_thr;
      ::derefine_thr = derefine_thr;
      
      std::cout << "AMR: Density gradient refinement enabled" << std::endl;
      std::cout << "     refine_thr   = " << refine_thr << std::endl;
      std::cout << "     derefine_thr = " << derefine_thr << std::endl;
    }
    
    // 両方のリファイン方式が有効な場合のメッセージ
    if (jeans_refine_input && grad_refine_input) {
      std::cout << "AMR: Using BOTH Jeans length AND density gradient refinement" << std::endl;
      std::cout << "     (refine if EITHER condition is met)" << std::endl;
    }
    
    // リファイン条件関数の登録
    EnrollUserRefinementCondition(RefinementCondition);

  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "gr_star");
  SetUserOutputVariableName(1, "gr_nfw");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Determine mesh center (default sphere center)
  // 計算領域が原点対称なので、これらは原点になる
  Real x0 = 0.5*(pmy_mesh->mesh_size.x1min + pmy_mesh->mesh_size.x1max);
  Real y0 = 0.5*(pmy_mesh->mesh_size.x2min + pmy_mesh->mesh_size.x2max);
  Real z0 = 0.5*(pmy_mesh->mesh_size.x3min + pmy_mesh->mesh_size.x3max);
  
  // isothermal sound speed
  Real cs = pin->GetReal("hydro","iso_sound_speed");

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
	Real x = pcoord->x1v(i) - x0;
	Real y = pcoord->x2v(j) - y0;
	Real z = pcoord->x3v(k) - z0;
	Real r_sph = std::sqrt(x*x + y*y +z*z);

        Real dx = std::min({
  	  pcoord->dx1v(i),
  	  pcoord->dx2v(j),
  	  pcoord->dx3v(k)
	});

	Real r_sph_safe = std::max(r_sph, 0.5*dx);

        Real rphys = r_sph_safe * Lunit;

	// Toyouchi+23 density
        Real rho_phys = 1.1e-19 * pow(rphys/(1.0e5*AU), -1.75);
	Real rho_profile = rho_phys / Rhounit;

	// isothermal pressure
	Real P_profile = rho_profile * cs * cs;

        Real rho_final = rho_profile;
	Real P_final   = P_profile;

        // 初期状態ではシンク内部を低密度にしておく。
        // この初期除去分は中心星への降着量には数えない。
        if (use_sink && r_sph < r_sink) {
          rho_final = sink_rho_floor;
          P_final   = sink_rho_floor * cs * cs;
        }

        phydro->u(IDN,k,j,i) = rho_final;

	user_out_var(0,k,j,i) = 0.0;   // star
	user_out_var(1,k,j,i) = 0.0;   // nfw

        // 追加 rotation velocity
        Real vx = 0.0, vy = 0.0, vz = 0.0;

        Real r_cyl = std::sqrt(x*x + y*y);

        // --- cutoff radius ---
        Real r_cut = std::max(racc, 1e-12); // 入力ファイルでracc=0.44としているが念の為下限値を設定

        // --- 回転速度（解析的M_encを使用）---
        // --- gas enclosed mass (Toyouchi+23) ---
	Real M_gas = 4.0 * M_PI * 0.11 * pow(r_cyl, 1.25) / 1.25;

	// --- Toyouchi rotation: 0.5 vkep (gas only) ---
	Real vphi_eff = 0.0;
	if (use_rotation && r_cyl > 0.0) {
    	    vphi_eff = 0.5 * std::sqrt(gconst * M_gas / r_cyl);
	}

        // 中心で滑らかにゼロへ（二次減衰で安定化）
        if (r_cyl < r_cut) {
            vphi_eff *= (r_cyl / r_cut) * (r_cyl / r_cut);
        }

        Real inv_r = 1.0 / std::max(r_cyl, 1e-12);
        vx += -vphi_eff * y * inv_r;
        vy +=  vphi_eff * x * inv_r;

        // --- spherical radial inflow toward the origin ---
        if (vr0 != 0.0) {
            Real vr_eff = vr0;
            // 原点近傍で速度を滑らかに0へ落とす
            if (r_sph < r_cut) {
                vr_eff *= r_sph / r_cut;
            }

            const Real inv_r_sph = 1.0 / r_sph_safe;
            vx += vr_eff * x * inv_r_sph;
            vy += vr_eff * y * inv_r_sph;
            vz += vr_eff * z * inv_r_sph;
        }

        // シンク領域内部の速度をゼロにする
        if (use_sink && r_sph < r_sink) {
          vx = 0.0;
          vy = 0.0;
          vz = 0.0;
        }

        //追加 set momentum
	phydro->u(IM1,k,j,i) = rho_final * vx;
	phydro->u(IM2,k,j,i) = rho_final * vy;
	phydro->u(IM3,k,j,i) = rho_final * vz;

        // set internal energy consistently (use P_final and kinetic energy)
        if (NON_BAROTROPIC_EOS) {
	  Real ke = 0.5 * rho_final *
         	    (vx*vx + vy*vy + vz*vz);
          Real me = 0.0;
          if (MAGNETIC_FIELDS_ENABLED) {
            // Athena++ではPmag,code = Bcode * Bcode / 2
            me = 0.5 * bz_code * bz_code;
          }
          phydro->u(IEN,k,j,i) = P_final/gm1 + ke + me;
	}
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    // Initialize face-centered magnetic field: B = (0, 0, bz_code)

    // Bx on x1 faces
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = 0.0;
        }
      }
    }

    // By on x2 faces
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }

    // Bz on x3 faces
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = bz_code;
        }
      }
    }
  }
}

void CentralGravity(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc,
                    AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar) {

  Coordinates *pcoord = pmb->pcoord;
  const Real x0 = 0.5 * (pmb->pmy_mesh->mesh_size.x1min
                       + pmb->pmy_mesh->mesh_size.x1max);
  const Real y0 = 0.5 * (pmb->pmy_mesh->mesh_size.x2min
                       + pmb->pmy_mesh->mesh_size.x2max);
  const Real z0 = 0.5 * (pmb->pmy_mesh->mesh_size.x3min
                       + pmb->pmy_mesh->mesh_size.x3max);

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        Real x = pcoord->x1v(i) - x0;
        Real y = pcoord->x2v(j) - y0;
        Real z = pcoord->x3v(k) - z0;

        Real r2 = x*x + y*y + z*z + epsilon_soft*epsilon_soft;
        Real r  = std::sqrt(r2);

	Real gx = 0.0;
	Real gy = 0.0;
	Real gz = 0.0;

        if (use_central_gravity) {
        // 降着により増加した中心星質量が次のタイムステップで中心星重力に反映
        Real mstar =
            pmb->pmy_mesh->ruser_mesh_data[0](SINK_MSTAR);
   	    Real fac = -gconst*mstar/(r2*r);

	    Real gr_star = fac * std::sqrt(x*x + y*y + z*z);

	    pmb->user_out_var(0,k,j,i) = gr_star;

   	    gx += fac*x;
   	    gy += fac*y;
   	    gz += fac*z;
	}


	if (use_nfw_gravity) {
   	    // radius in CGS
	    Real rphys = r * Lunit;

	    // x = r/Rs
	    Real xdm = rphys / Rs_phys;

	    if (xdm > 1.0e-12) {

    	        Real bracket =
        	    xdm/(1.0 + xdm)
        	    - std::log(1.0 + xdm);

    		// NFW acceleration (CGS)
    		Real gr_phys =
        	    (Phi0_phys / Rs_phys)
        	    * bracket
        	    / (xdm * xdm);

    		// CGS -> code unit
    		Real gr =
        	    gr_phys * (Tunit * Tunit / Lunit);

    		pmb->user_out_var(1,k,j,i) = gr;

    		gx += gr * (x/r);
    		gy += gr * (y/r);
    		gz += gr * (z/r);
	     }

	 }

	 // ===== Debug print (before source update) =====
        if (Globals::my_rank == 0 &&
            pmb->gid == 435 &&
            i == 0 && j == 31 && k == 31) {

          std::cout
            << "\n===== DEBUG CentralGravity =====\n"
            << "time = " << time
            << "  dt = " << dt << "\n"
            << "gid = " << pmb->gid << "\n"
            << "r = " << r << "\n"
            << "gx = " << gx
            << " gy = " << gy
            << " gz = " << gz << "\n"
            << "rho = " << prim(IDN,k,j,i) << "\n"
            << "before M = ("
            << cons(IM1,k,j,i) << ", "
            << cons(IM2,k,j,i) << ", "
            << cons(IM3,k,j,i) << ")"
            << std::endl;
        }

        cons(IM1,k,j,i) += dt * prim(IDN,k,j,i) * gx;
        cons(IM2,k,j,i) += dt * prim(IDN,k,j,i) * gy;
        cons(IM3,k,j,i) += dt * prim(IDN,k,j,i) * gz;

	// ===== Debug print (after source update) =====
        if (Globals::my_rank == 0 &&
            pmb->gid == 435 &&
            i == 0 && j == 31 && k == 31) {

          std::cout
            << "after  M = ("
            << cons(IM1,k,j,i) << ", "
            << cons(IM2,k,j,i) << ", "
            << cons(IM3,k,j,i) << ")"
            << std::endl;
        }
      }
    }
  }
}
  //  pmy_mesh->tlim=pin->SetReal("time","tlim",TWO_PI/omega*2.0);

void Mesh::UserWorkInLoop() {
  if (!use_sink) {
    ruser_mesh_data[0](SINK_MDOT_FLUX) = 0.0;
    ruser_mesh_data[0](SINK_MDOT_RESET) = 0.0;
    ruser_mesh_data[0](SINK_MDOT_FLOOR) = 0.0;
    ruser_mesh_data[0](SINK_MGAS) = 0.0;
    return;
  }

  // NOTE: この境界flux積分は現在使用中のintegrator=vl2を前提とする。
  // 多段Runge-Kutta法へ変更する場合は、各stageの重み付きfluxを積算する必要がある。

  // sink境界を外側から内側へ通過した正味質量 [code mass]
  Real flux_mass_local = 0.0;

  // sink内部を目標密度へ戻す際に除去した質量 [code mass]
  // これは診断量であり、中心星質量には加えない。
  Real reset_removed_mass_local = 0.0;

  // Alfven速度floorにより人工的に追加した質量 [code mass]
  Real floor_added_mass_local = 0.0;

  // リセット後にsink内部へ残るガス質量 [code mass]
  Real sink_gas_mass_local = 0.0;

  const Real x0 = 0.5 * (mesh_size.x1min + mesh_size.x1max);
  const Real y0 = 0.5 * (mesh_size.x2min + mesh_size.x2max);
  const Real z0 = 0.5 * (mesh_size.x3min + mesh_size.x3max);

  // --------------------------------------------------------------------------
  // 1. sink内外をまたぐCartesian faceの数値質量fluxを積分する。
  // sink内部セル側から各faceを一度だけ数え、二重計数を避ける。
  // phydro->flux[dir](IDN,...)は単位面積・単位時間当たりの質量flux。
  // --------------------------------------------------------------------------
  for (int n = 0; n < nblocal; ++n) {
    MeshBlock *pmb = my_blocks(n);
    Coordinates *pcoord = pmb->pcoord;
    AthenaArray<Real> &x1flux = pmb->phydro->flux[X1DIR];
    AthenaArray<Real> &x2flux = pmb->phydro->flux[X2DIR];
    AthenaArray<Real> &x3flux = pmb->phydro->flux[X3DIR];

    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          const Real x = pcoord->x1v(i) - x0;
          const Real y = pcoord->x2v(j) - y0;
          const Real z = pcoord->x3v(k) - z0;
          const Real r = std::sqrt(x*x + y*y + z*z);

          if (r >= r_sink) {
            continue;
          }

          // x1負側: +x1向きfluxがsinkへ入るので正。
          const Real xm = pcoord->x1v(i-1) - x0;
          if (std::sqrt(xm*xm + y*y + z*z) >= r_sink) {
            flux_mass_local += x1flux(IDN,k,j,i)
                * pcoord->GetFace1Area(k,j,i) * dt;
          }

          // x1正側: -x1向きfluxがsinkへ入るので符号を反転。
          const Real xp = pcoord->x1v(i+1) - x0;
          if (std::sqrt(xp*xp + y*y + z*z) >= r_sink) {
            flux_mass_local -= x1flux(IDN,k,j,i+1)
                * pcoord->GetFace1Area(k,j,i+1) * dt;
          }

          // x2負側: +x2向きfluxがsinkへ入るので正。
          const Real ym = pcoord->x2v(j-1) - y0;
          if (std::sqrt(x*x + ym*ym + z*z) >= r_sink) {
            flux_mass_local += x2flux(IDN,k,j,i)
                * pcoord->GetFace2Area(k,j,i) * dt;
          }

          // x2正側: -x2向きfluxがsinkへ入るので符号を反転。
          const Real yp = pcoord->x2v(j+1) - y0;
          if (std::sqrt(x*x + yp*yp + z*z) >= r_sink) {
            flux_mass_local -= x2flux(IDN,k,j+1,i)
                * pcoord->GetFace2Area(k,j+1,i) * dt;
          }

          // x3負側: +x3向きfluxがsinkへ入るので正。
          const Real zm = pcoord->x3v(k-1) - z0;
          if (std::sqrt(x*x + y*y + zm*zm) >= r_sink) {
            flux_mass_local += x3flux(IDN,k,j,i)
                * pcoord->GetFace3Area(k,j,i) * dt;
          }

          // x3正側: -x3向きfluxがsinkへ入るので符号を反転。
          const Real zp = pcoord->x3v(k+1) - z0;
          if (std::sqrt(x*x + y*y + zp*zp) >= r_sink) {
            flux_mass_local -= x3flux(IDN,k+1,j,i)
                * pcoord->GetFace3Area(k+1,j,i) * dt;
          }
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // 2. sink内部を数値的な目標状態へ戻す。
  // 磁場は変更せず、密度だけを上げてvA=|B|/sqrt(rho)を制限する。
  // --------------------------------------------------------------------------
  for (int n = 0; n < nblocal; ++n) {
    MeshBlock *pmb = my_blocks(n);
    AthenaArray<Real> &u = pmb->phydro->u;
    AthenaArray<Real> &w = pmb->phydro->w;
#if MAGNETIC_FIELDS_ENABLED
    // face-centered磁場から作られたcell-centered磁場 [code magnetic field]
    // density floorの評価にだけ使い、磁場自体は変更しない。
    AthenaArray<Real> &bcc = pmb->pfield->bcc;
#endif

    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          const Real x = pmb->pcoord->x1v(i) - x0;
          const Real y = pmb->pcoord->x2v(j) - y0;
          const Real z = pmb->pcoord->x3v(k) - z0;
          const Real r = std::sqrt(x*x + y*y + z*z);
          const Real dx = std::min({pmb->pcoord->dx1v(i),
                                    pmb->pcoord->dx2v(j),
                                    pmb->pcoord->dx3v(k)});

          if (r < r_sink) {
            const Real rho_old = u(IDN,k,j,i);  // リセット前密度 [code density]
            Real rho_target = sink_rho_floor;   // sinkに残す密度 [code density]

#if MAGNETIC_FIELDS_ENABLED
            if (sink_va_cap > 0.0) {
              const Real bx = bcc(IB1,k,j,i);  // cell-centered Bx
              const Real by = bcc(IB2,k,j,i);  // cell-centered By
              const Real bz = bcc(IB3,k,j,i);  // cell-centered Bz
              const Real b2 = bx*bx + by*by + bz*bz;  // |B|^2

              // vA=|B|/sqrt(rho)<=sink_va_capに必要な最低密度 [code density]
              const Real rho_va_floor = b2/SQR(sink_va_cap);
              rho_target = std::max(rho_target, rho_va_floor);
            }
#endif

            const Real cell_volume = pmb->pcoord->GetCellVolume(k,j,i);
            if (rho_old > rho_target) {
              reset_removed_mass_local += (rho_old-rho_target)*cell_volume;
            } else if (rho_old < rho_target) {
              floor_added_mass_local += (rho_target-rho_old)*cell_volume;
            }

            // isothermal EOSなのでエネルギー変数はない。
            u(IDN,k,j,i) = rho_target;
            u(IM1,k,j,i) = 0.0;
            u(IM2,k,j,i) = 0.0;
            u(IM3,k,j,i) = 0.0;
            w(IDN,k,j,i) = rho_target;
            w(IVX,k,j,i) = 0.0;
            w(IVY,k,j,i) = 0.0;
            w(IVZ,k,j,i) = 0.0;

            sink_gas_mass_local += rho_target*cell_volume;
            continue;
          }

          // sink表面外側で、sink中心から外向きの速度成分だけを除去する。
          const Real shell_outer = r_sink + sink_no_outflow_cells*dx;
          if (r < shell_outer && r > 0.0) {
            const Real inv_r = 1.0/r;
            const Real vx = w(IVX,k,j,i);
            const Real vy = w(IVY,k,j,i);
            const Real vz = w(IVZ,k,j,i);
            const Real vr = (vx*x + vy*y + vz*z)*inv_r;  // 動径速度

            if (vr > 0.0) {
              const Real vx_new = vx-vr*x*inv_r;
              const Real vy_new = vy-vr*y*inv_r;
              const Real vz_new = vz-vr*z*inv_r;
              w(IVX,k,j,i) = vx_new;
              w(IVY,k,j,i) = vy_new;
              w(IVZ,k,j,i) = vz_new;
              const Real rho = u(IDN,k,j,i);
              u(IM1,k,j,i) = rho*vx_new;
              u(IM2,k,j,i) = rho*vy_new;
              u(IM3,k,j,i) = rho*vz_new;
            }
          }
        }
      }
    }
  }

  // 全MPI rankで、flux・cell reset・floor・sink gas massを合算する。
  Real sink_data[4] = {flux_mass_local, reset_removed_mass_local,
                       floor_added_mass_local, sink_gas_mass_local};
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, sink_data, 4, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  const Real flux_mass_global = sink_data[0];
  const Real reset_removed_mass_global = sink_data[1];
  const Real floor_added_mass_global = sink_data[2];
  const Real sink_gas_mass_global = sink_data[3];

  // 中心星質量はsink境界からの正味流入だけで更新する。
  // 微小な外向き正味fluxで中心星質量が減らないよう0で下限を取る。
  const Real accreted_mass = std::max(flux_mass_global, static_cast<Real>(0.0));
  ruser_mesh_data[0](SINK_MSTAR) += accreted_mass;

  if (dt > 0.0) {
    ruser_mesh_data[0](SINK_MDOT_FLUX) = accreted_mass/dt;
    ruser_mesh_data[0](SINK_MDOT_RESET) = reset_removed_mass_global/dt;
    ruser_mesh_data[0](SINK_MDOT_FLOOR) = floor_added_mass_global/dt;
  } else {
    ruser_mesh_data[0](SINK_MDOT_FLUX) = 0.0;
    ruser_mesh_data[0](SINK_MDOT_RESET) = 0.0;
    ruser_mesh_data[0](SINK_MDOT_FLOOR) = 0.0;
  }
  ruser_mesh_data[0](SINK_MGAS) = sink_gas_mass_global;

  if (Globals::my_rank == 0 && ncycle_out > 0 && ncycle % ncycle_out == 0) {
    std::cout << "Sink: Mstar=" << ruser_mesh_data[0](SINK_MSTAR)
              << " Mdot_flux=" << ruser_mesh_data[0](SINK_MDOT_FLUX)
              << " dM_flux=" << accreted_mass
              << " Mdot_reset=" << ruser_mesh_data[0](SINK_MDOT_RESET)
              << " Mdot_floor=" << ruser_mesh_data[0](SINK_MDOT_FLOOR)
              << " Msink_gas=" << ruser_mesh_data[0](SINK_MGAS)
              << std::endl;
  }
}

// hstへsink関連の大域物理量・数値診断量を出力する。
Real SinkHistory(MeshBlock *pmb, int iout) {
  // Historyは全MeshBlockの返り値を加算するため、大域量はgid=0だけから返す。
  if (pmb->gid != 0) {
    return 0.0;
  }

  if (iout == 0) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MSTAR);       // [code mass]
  }
  if (iout == 1) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MDOT_FLUX);  // [mass/time]
  }
  if (iout == 2) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MDOT_RESET); // [mass/time]
  }
  if (iout == 3) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MDOT_FLOOR); // [mass/time]
  }
  if (iout == 4) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MGAS);       // [code mass]
  }
  return 0.0;
}

int RefinementCondition(MeshBlock *pmb) {

  bool need_refine = false;
  bool need_derefine = true;

  // シンク周辺の強制リファイン（シンク表面周辺ではセルを細かくしておきたい）
  const Real x0 =
      0.5 * (pmb->pmy_mesh->mesh_size.x1min
           + pmb->pmy_mesh->mesh_size.x1max);
  const Real y0 =
      0.5 * (pmb->pmy_mesh->mesh_size.x2min
           + pmb->pmy_mesh->mesh_size.x2max);
  const Real z0 =
      0.5 * (pmb->pmy_mesh->mesh_size.x3min
           + pmb->pmy_mesh->mesh_size.x3max);

  // セル中心ではなく、MeshBlockとシンク球の最短距離で判定する。
  // これにより、基本格子がシンクより粗くても中心ブロックの細分化を開始できる。
  const Real xmin = pmb->block_size.x1min - x0;
  const Real xmax = pmb->block_size.x1max - x0;
  const Real ymin = pmb->block_size.x2min - y0;
  const Real ymax = pmb->block_size.x2max - y0;
  const Real zmin = pmb->block_size.x3min - z0;
  const Real zmax = pmb->block_size.x3max - z0;

  Real dx_block = 0.0;
  Real dy_block = 0.0;
  Real dz_block = 0.0;

  if (xmin > 0.0) {
    dx_block = xmin;
  } else if (xmax < 0.0) {
    dx_block = -xmax;
  }
  if (ymin > 0.0) {
    dy_block = ymin;
  } else if (ymax < 0.0) {
    dy_block = -ymax;
  }
  if (zmin > 0.0) {
    dz_block = zmin;
  } else if (zmax < 0.0) {
    dz_block = -zmax;
  }

  const Real block_min_radius =
      std::sqrt(dx_block*dx_block + dy_block*dy_block + dz_block*dz_block);
  const bool sink_region_in_block =
      use_sink && block_min_radius < r_sink + sink_refine_buffer;

  if (sink_region_in_block) {
    need_refine = true;
    need_derefine = false;
  }

  Real cs_iso = std::sqrt(cs2);

  Real gradmax = 0.0;
  Real njmin   = 1e30;

  for (int k = pmb->ks+1; k <= pmb->ke-1; ++k) {
    for (int j = pmb->js+1; j <= pmb->je-1; ++j) {
      for (int i = pmb->is+1; i <= pmb->ie-1; ++i) {

        Real rho = pmb->phydro->w(IDN, k, j, i);

        Real rho_safe = std::max(rho, 1e-20);

        // 音速
        Real cs = (NON_BAROTROPIC_EOS)
                    ? std::sqrt(pmb->phydro->w(IPR, k, j, i) / rho_safe)
                    : cs_iso;

        // ===== Jeans判定（collapse.cpp風）=====
        Real dx = std::min({
          pmb->pcoord->dx1v(i),
          pmb->pcoord->dx2v(j),
          pmb->pcoord->dx3v(k)
        });

        Real nj = cs / std::sqrt(rho_safe);
        nj *= (2.0 * M_PI / dx);

        // ★minは記録だけ（derefine用）
        njmin = std::min(njmin, nj);

        // ★これが超重要：局所で即refine
        if (use_jeans_refine && nj < jeans_cells) {
          need_refine = true;
        }

        // ===== 密度勾配 =====
        if (use_grad_refine) {
          Real drho = 0.5 * std::abs(
            pmb->phydro->w(IDN, k, j, i+1) -
            pmb->phydro->w(IDN, k, j, i-1)
          );

          Real grad = drho / rho_safe;
          gradmax = std::max(gradmax, grad);
        }
      }
    }
  }

  // ===== gradによるrefine =====
  if (use_grad_refine && gradmax > refine_thr) {
    need_refine = true;
  }

  // ===== derefine条件（少し緩め）=====
  if (use_jeans_refine) {
    if (njmin < jeans_cells * 1.2) {
      need_derefine = false;
    }
  }

  if (use_grad_refine && gradmax > derefine_thr) {  // gradによるderefine（ヒステリシス）
    need_derefine = false;
  }

  if (need_refine) return 1;  // 細かくする
  else if (need_derefine) return -1;  // 粗くする
  else return 0;  // そのまま現状維持
}

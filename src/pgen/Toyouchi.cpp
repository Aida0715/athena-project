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

// シンク表面の外側で逆流を抑制する幅
// 各セルで sink_no_outflow_cells * dx として使用
Real sink_no_outflow_cells = 1.0;

// AMRで強制的に細分化するシンク外側の幅
Real sink_refine_buffer = 1.0;

// 診断量を格納する添字
enum SinkDataIndex {
  SINK_MSTAR = 0, // 現在の中心星質量
  SINK_MDOT  = 1  // 直前のタイムステップの降着率
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
  ruser_mesh_data[0].NewAthenaArray(2);
  ruser_mesh_data[0](SINK_MSTAR) = mstar_init;
  ruser_mesh_data[0](SINK_MDOT)  = 0.0;

  // 履歴出力を登録する
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, SinkHistory, "Mstar");
  EnrollUserHistoryOutput(1, SinkHistory, "Mdot_sink");

  // 初期設定をログに表示
  if (Globals::my_rank == 0) {
    std::cout
        << "Sink region:" << std::endl
        << "  enabled          = " << use_sink << std::endl
        << "  r_sink           = " << r_sink << " code units" << std::endl
        << "  r_sink           = " << r_sink_au << " AU" << std::endl
        << "  rho floor        = " << sink_rho_floor << std::endl
        << "  initial Mstar    = " << mstar_init << " code units"
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

        // --- radial inflow ---
        if (vr0 != 0.0) {
            Real vr_eff = vr0;

            if (r_cyl < r_cut) {
                vr_eff *= r_cyl / r_cut;
            }

            Real inv_r = 1.0 / std::max(r_cyl, 1e-12);
            vx += vr_eff * x * inv_r;
            vy += vr_eff * y * inv_r;
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
      ruser_mesh_data[0](SINK_MDOT) = 0.0;
      return;
    }
  
    Real removed_mass_local = 0.0;
  
    const Real x0 = 0.5 * (mesh_size.x1min + mesh_size.x1max);
    const Real y0 = 0.5 * (mesh_size.x2min + mesh_size.x2max);
    const Real z0 = 0.5 * (mesh_size.x3min + mesh_size.x3max);
  
    // このMPIランクが所有する全MeshBlockを処理する
    for (int n = 0; n < nblocal; ++n) {
      MeshBlock *pmb = my_blocks(n);
  
      AthenaArray<Real> &u = pmb->phydro->u;
      AthenaArray<Real> &w = pmb->phydro->w;
  
      for (int k = pmb->ks; k <= pmb->ke; ++k) {
        for (int j = pmb->js; j <= pmb->je; ++j) {
          for (int i = pmb->is; i <= pmb->ie; ++i) {
            const Real x = pmb->pcoord->x1v(i) - x0;
            const Real y = pmb->pcoord->x2v(j) - y0;
            const Real z = pmb->pcoord->x3v(k) - z0;
  
            const Real r = std::sqrt(x*x + y*y + z*z);
  
            const Real dx = std::min({
                pmb->pcoord->dx1v(i),
                pmb->pcoord->dx2v(j),
                pmb->pcoord->dx3v(k)
            });
  
            // ============================================================
            // 1. シンク内部：ガスを密度フロアまで吸収
            // ============================================================
            if (r < r_sink) {
              const Real rho_old = u(IDN,k,j,i);
              const Real rho_new = sink_rho_floor;
  
              if (rho_old > rho_new) {
                const Real cell_volume =
                    pmb->pcoord->GetCellVolume(k,j,i);
  
                removed_mass_local +=
                    (rho_old - rho_new) * cell_volume;
              }
  
              // isothermal EOSなのでエネルギー変数はない
              u(IDN,k,j,i) = rho_new;
              u(IM1,k,j,i) = 0.0;
              u(IM2,k,j,i) = 0.0;
              u(IM3,k,j,i) = 0.0;
  
              // Conserved変数だけでなくPrimitive変数も更新する
              w(IDN,k,j,i) = rho_new;
              w(IVX,k,j,i) = 0.0;
              w(IVY,k,j,i) = 0.0;
              w(IVZ,k,j,i) = 0.0;
  
              continue;
            }
  
            // ============================================================
            // 2. シンク表面の外側：外向き速度成分を禁止
            // ============================================================
            const Real shell_outer =
                r_sink + sink_no_outflow_cells * dx;
  
            if (r < shell_outer && r > 0.0) {
              const Real inv_r = 1.0 / r;
  
              const Real vx = w(IVX,k,j,i);
              const Real vy = w(IVY,k,j,i);
              const Real vz = w(IVZ,k,j,i);
  
              const Real vr =
                  (vx*x + vy*y + vz*z) * inv_r;
  
              // vr > 0 はシンク中心から外向き
              if (vr > 0.0) {
                const Real vx_new = vx - vr*x*inv_r;
                const Real vy_new = vy - vr*y*inv_r;
                const Real vz_new = vz - vr*z*inv_r;
  
                w(IVX,k,j,i) = vx_new;
                w(IVY,k,j,i) = vy_new;
                w(IVZ,k,j,i) = vz_new;
  
                const Real rho = u(IDN,k,j,i);
  
                u(IM1,k,j,i) = rho * vx_new;
                u(IM2,k,j,i) = rho * vy_new;
                u(IM3,k,j,i) = rho * vz_new;
              }
            }
          }
        }
      }
    }
  
    // 全MPIランクの除去質量を合算する
    Real removed_mass_global = removed_mass_local;
  
  #ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &removed_mass_global, 1,
                  MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  #endif
  
    // 中心星質量を更新
    ruser_mesh_data[0](SINK_MSTAR) += removed_mass_global;
  
    // コード単位での瞬間降着率
    if (dt > 0.0) {
      ruser_mesh_data[0](SINK_MDOT) =
          removed_mass_global / dt;
    } else {
      ruser_mesh_data[0](SINK_MDOT) = 0.0;
    }
  
    if (Globals::my_rank == 0
        && ncycle_out > 0
        && ncycle % ncycle_out == 0) {
      std::cout
          << "Sink: Mstar = "
          << ruser_mesh_data[0](SINK_MSTAR)
          << ", Mdot = "
          << ruser_mesh_data[0](SINK_MDOT)
          << ", dM = "
          << removed_mass_global
          << std::endl;
    }
  }

// hst出力関数を追加する（hstにはコード単位で、Mstar、Mdot_sinkが追加される）
Real SinkHistory(MeshBlock *pmb, int iout) {
  // History出力は全MeshBlockについて和を取る。
  // 同じMstarを全ブロックから返すとブロック数倍になるため、
  // gid == 0 のブロックだけが値を返す。
  if (pmb->gid != 0) {
    return 0.0;
  }

  if (iout == 0) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MSTAR);
  }

  if (iout == 1) {
    return pmb->pmy_mesh->ruser_mesh_data[0](SINK_MDOT);
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

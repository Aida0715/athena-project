//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jeans.cpp
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

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

void CentralGravity(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc,
                    AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar);

int RefinementCondition(MeshBlock *pmb); // AMRのリファインメント宣言

namespace {
// with functions A1,2,3 which compute vector potentials
Real cs2, gam, gm1, gconst;

// AMR関連のグローバル変数（RefinementCondition関数で使用）
bool use_jeans_refine = true;   // Jeans長ベースのリファインを使用するか
bool use_grad_refine = false;   // 密度勾配ベースのリファインを使用するか
Real jeans_cells = 8.0;         // Jeans長を何セルで解像するか
Real refine_thr = 0.3;          // 密度勾配の閾値（use_grad_refine=true時）

//追加 rotation parameters
Real vphi0 = 0.0;   // cylindrical rotational velocity
Real vr0   = 0.0;   // radial velocity

//===== Toyouchi+23 parameters =====
Real epsilon_soft = 0.5;   // gravitational softening
const Real Mstar = 1.0;    // central star mass (code unit:2solar_mass in CGS)

// ===== unit system (code unit <-> cgs) =====
const Real Munit = 4.0e33;   // g (central star mass)
const Real Lunit = 6.7e15;   // cm (bondi radius)
const Real Tunit = 3.34e10;  // s (free fall time)

const Real Vunit = Lunit/Tunit;                 // velocity unit
const Real Rhounit = Munit/(Lunit*Lunit*Lunit); // density unit
const Real Punit = Rhounit*Vunit*Vunit;         // pressure unit

// unit conversion
const Real AU = 1.496e13;

// 中心カットオフ用
Real racc = 0.0;   // accretion radius (cutoff scale)

} // namespace


void Mesh::InitUserMeshData(ParameterInput *pin) {

  epsilon_soft = pin->GetOrAddReal("problem","epsilon",0.5);
  gconst = pin->GetOrAddReal("problem","grav_const",1.0);

  EnrollUserExplicitSourceFunction(CentralGravity);

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
  vphi0 = pin->GetOrAddReal("problem","vphi0", 0.0);
  vr0   = pin->GetOrAddReal("problem","vr0", 0.0);

  // raccを入力ファイルから
  racc = pin->GetOrAddReal("problem", "racc", 0.44);

  //重力定数（入力ファイルから読み込む）
  if (SELF_GRAVITY_ENABLED) {
    SetGravitationalConstant(gconst);
  }

  // AMR関連の設定
  if (adaptive) {
    // AMRリファイン方式の選択（入力ファイルから読み込み）
    // デフォルト：Jeans長ベースのリファインを使用（自己重力系では必須）
    bool use_jeans_refine = pin->GetOrAddBoolean("problem", "use_jeans_refine", true);
    bool use_grad_refine  = pin->GetOrAddBoolean("problem", "use_grad_refine", false);
    
    // 入力パラメータのチェックと警告
    if (!use_jeans_refine && !use_grad_refine) {
      std::cout << "WARNING: AMR enabled but no refinement condition specified!" << std::endl;
      std::cout << "         Enabling Jeans length refinement as default." << std::endl;
      use_jeans_refine = true;
    }
    
    // Jeans長ベースのリファインが有効な場合のパラメータ読み込み
    if (use_jeans_refine) {
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
      ::use_jeans_refine = use_jeans_refine;
      ::jeans_cells = jeans_cells;
      
      std::cout << "AMR: Jeans length refinement enabled" << std::endl;
      std::cout << "     Cells per Jeans length = " << jeans_cells << std::endl;
    }
    
    // 密度勾配ベースのリファインが有効な場合のパラメータ読み込み
    if (use_grad_refine) {
      // 密度勾配の閾値（必須パラメータ）
      refine_thr = pin->GetReal("problem", "refine_thr");
      
      // グローバル変数に保存
      ::use_grad_refine = use_grad_refine;
      ::refine_thr = refine_thr;
      
      std::cout << "AMR: Density gradient refinement enabled" << std::endl;
      std::cout << "     Gradient threshold = " << refine_thr << std::endl;
    }
    
    // 両方のリファイン方式が有効な場合のメッセージ
    if (use_jeans_refine && use_grad_refine) {
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

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Determine mesh center (default sphere center)
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

        phydro->u(IDN,k,j,i) = rho_final;

        // 追加 rotation velocity
        Real vx = 0.0, vy = 0.0, vz = 0.0;

        Real r_cyl = std::sqrt(x*x + y*y);

        // --- cutoff radius ---
        Real r_cut = std::max(racc, 1e-12); // 入力ファイルでracc=0.44としているが念の為下限値を設定

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

        //追加 set momentum
	phydro->u(IM1,k,j,i) = rho_final * vx;
	phydro->u(IM2,k,j,i) = rho_final * vy;
	phydro->u(IM3,k,j,i) = rho_final * vz;

        // set internal energy consistently (use P_final and kinetic energy)
        if (NON_BAROTROPIC_EOS) {
	  Real ke = 0.5 * rho_final *
         	    (vx*vx + vy*vy + vz*vz);
 	  phydro->u(IEN,k,j,i) = P_final/gm1 + ke;
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

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real r2 = x*x + y*y + z*z + epsilon_soft*epsilon_soft;
        Real r  = std::sqrt(r2);

        Real fac = -gconst * Mstar / (r2 * r);

        Real gx = fac * x;
        Real gy = fac * y;
        Real gz = fac * z;

        cons(IM1,k,j,i) += dt * prim(IDN,k,j,i) * gx;
        cons(IM2,k,j,i) += dt * prim(IDN,k,j,i) * gy;
        cons(IM3,k,j,i) += dt * prim(IDN,k,j,i) * gz;
      }
    }
  }
}
  //  pmy_mesh->tlim=pin->SetReal("time","tlim",TWO_PI/omega*2.0);

int RefinementCondition(MeshBlock *pmb) {

  bool need_refine = false;
  bool need_derefine = true;

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

  if (use_grad_refine) {
    if (gradmax > 0.5 * refine_thr) {
      need_derefine = false;
    }
  }

  if (need_refine) return 1;
  if (need_derefine) return -1;
  return 0;
}

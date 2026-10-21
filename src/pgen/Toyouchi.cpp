//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file Toyouchi.cpp
//! \brief Wise19 MMH cloud with fixed stellar gravity and gas self-gravity.

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <string>     // c_str()
#include <fstream>    // 追加
#include <iomanip>    // setprecision
#include <vector>     // 追加
#include <utility>    // 追加
#include <cstdint>    // white noise

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
#include "../units/units.hpp"

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

Real StarMassHistory(MeshBlock *pmb, int iout);

namespace {
// with functions A1,2,3 which compute vector potentials
Real cs2, gam, gm1, gconst;

// 初期密度ノイズ
bool use_density_noise = false;
Real density_noise_percent = 0.0;
int density_noise_seed = 42;

// 64bit整数を混ぜるハッシュ関数。
// unsigned整数のオーバーフローは意図した動作。
std::uint64_t DensityNoiseHash(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

// セルごとの擬似乱数を約[-1, 1)の範囲で返す。
Real DensityNoiseValue(std::uint64_t ix,
                       std::uint64_t iy,
                       std::uint64_t iz,
                       int refinement_level) {
  std::uint64_t h =
      DensityNoiseHash(static_cast<std::uint64_t>(density_noise_seed));

  h = DensityNoiseHash(h ^ static_cast<std::uint64_t>(refinement_level));
  h = DensityNoiseHash(h ^ ix);
  h = DensityNoiseHash(h ^ iy);
  h = DensityNoiseHash(h ^ iz);

  const double uniform01 =
      static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);

  return static_cast<Real>(2.0 * uniform01 - 1.0);
}

// rotation parameters
bool use_radial_inflow = true;  // radial velocity
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
const Real Msun = 1.98847e33;       // g
const Real AU   = 1.495978707e13;   // cm
// Initialized from Mesh::punit before any MeshBlocks are constructed.
Real Munit, Lunit, Tunit;

const Real fit_rmin_au = 139.9720348;
const Real fit_rmax_au = 146698.0663;
Real Vunit, Rhounit;
Real mstar_fixed = 1.0;  // Fixed stellar mass in code units.
Real vr_taper_au = 580.0;
Real core_a, core_b, core_c, core_d, core_mass_msun;
Real taper_velocity_kms, taper_a, taper_b;
extern const Real mass_coeff[8];
extern const Real velocity_coeff[8];
extern Real bz_microgauss;
extern Real bz_code;
extern Real Bunit;
bool metadata_written = false;

std::string JsonEscape(const std::string &input) {
  std::ostringstream escaped;
  for (std::string::const_iterator it = input.begin(); it != input.end(); ++it) {
    switch (*it) {
      case '"': escaped << "\\\""; break;
      case '\\': escaped << "\\\\"; break;
      case '\n': escaped << "\\n"; break;
      case '\r': escaped << "\\r"; break;
      case '\t': escaped << "\\t"; break;
      default: escaped << *it; break;
    }
  }
  return escaped.str();
}

const char *JsonBool(bool value) { return value ? "true" : "false"; }

void WriteToyouchiMetadata(Mesh *pm, ParameterInput *pin,
                           Real sound_speed_kms, Real epsilon_au) {
  if (Globals::my_rank != 0 || metadata_written) return;

  const std::string problem_id = pin->GetString("job", "problem_id");
  const std::string filename = problem_id + ".metadata.json";
  const std::string temporary = filename + ".tmp";
  const Real sound_speed_code = sound_speed_kms * 1.0e5 / Vunit;
  const Real density_floor_code = pin->GetReal("hydro", "dfloor");
  const Real cfl_number = pin->GetReal("time", "cfl_number");
  const Real gamma_input = pin->GetReal("hydro", "gamma");
  const Real pressure_unit = Rhounit * Vunit * Vunit;
  const Real acceleration_unit = Lunit / (Tunit*Tunit);

  std::ofstream out(temporary.c_str(), std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    std::stringstream msg;
    msg << "Could not open Toyouchi metadata file: " << temporary;
    ATHENA_ERROR(msg);
  }
  out << std::setprecision(17);
  out << "{\n"
      << "  \"format\": \"Toyouchi runtime metadata\",\n"
      << "  \"format_version\": 1,\n"
      << "  \"problem_id\": \"" << JsonEscape(problem_id) << "\",\n"
      << "  \"units\": {\n"
      << "    \"unit_system\": \"custom\",\n"
      << "    \"mass_cgs\": " << Munit << ",\n"
      << "    \"length_cgs\": " << Lunit << ",\n"
      << "    \"time_cgs\": " << Tunit << ",\n"
      << "    \"velocity_cgs\": " << Vunit << ",\n"
      << "    \"density_cgs\": " << Rhounit << ",\n"
      << "    \"pressure_cgs\": " << pressure_unit << ",\n"
      << "    \"magnetic_field_gauss\": " << Bunit << ",\n"
      << "    \"acceleration_cgs\": " << acceleration_unit << ",\n"
      << "    \"gravity_constant_code\": " << gconst << "\n"
      << "  },\n"
      << "  \"hydro\": {\n"
      << "    \"eos\": \"" << (NON_BAROTROPIC_EOS ? "adiabatic" : "isothermal") << "\",\n"
      << "    \"gamma_input\": " << gamma_input << ",\n"
      << "    \"sound_speed_kms\": " << sound_speed_kms << ",\n"
      << "    \"sound_speed_code\": " << sound_speed_code << ",\n"
      << "    \"density_floor_code\": " << density_floor_code << ",\n"
      << "    \"cfl_number\": " << cfl_number << "\n"
      << "  },\n"
      << "  \"gravity\": {\n"
      << "    \"central_enabled\": " << JsonBool(use_central_gravity) << ",\n"
      << "    \"nfw_enabled\": " << JsonBool(use_nfw_gravity) << ",\n"
      << "    \"self_gravity_enabled\": " << JsonBool(SELF_GRAVITY_ENABLED) << ",\n"
      << "    \"mstar_code\": " << mstar_fixed << ",\n"
      << "    \"mstar_msun\": " << mstar_fixed*Munit/Msun << ",\n"
      << "    \"softening_au\": " << epsilon_au << ",\n"
      << "    \"softening_code\": " << epsilon_soft << "\n"
      << "  },\n"
      << "  \"initial_conditions\": {\n"
      << "    \"rotation_enabled\": " << JsonBool(use_rotation) << ",\n"
      << "    \"radial_inflow_enabled\": " << JsonBool(use_radial_inflow) << ",\n"
      << "    \"density_noise_enabled\": " << JsonBool(use_density_noise) << ",\n"
      << "    \"density_noise_percent\": " << density_noise_percent << ",\n"
      << "    \"density_noise_seed\": " << density_noise_seed << ",\n"
      << "    \"fit_rmin_au\": " << fit_rmin_au << ",\n"
      << "    \"fit_rmax_au\": " << fit_rmax_au << ",\n"
      << "    \"vr_taper_au\": " << vr_taper_au << ",\n"
      << "    \"mass_fit_coefficients\": [";
  for (int n = 0; n < 8; ++n) out << (n == 0 ? "" : ", ") << mass_coeff[n];
  out << "],\n    \"velocity_fit_coefficients\": [";
  for (int n = 0; n < 8; ++n) out << (n == 0 ? "" : ", ") << velocity_coeff[n];
  out << "],\n"
      << "    \"core_mass_msun\": " << core_mass_msun << ",\n"
      << "    \"core_coefficients\": [" << core_a << ", " << core_b << ", "
      << core_c << ", " << core_d << "],\n"
      << "    \"velocity_taper_coefficients\": [" << taper_velocity_kms << ", "
      << taper_a << ", " << taper_b << "]\n"
      << "  },\n"
      << "  \"magnetic_field\": {\n"
      << "    \"enabled_at_compile_time\": " << JsonBool(MAGNETIC_FIELDS_ENABLED) << ",\n"
      << "    \"initial_bz_microgauss\": " << bz_microgauss << ",\n"
      << "    \"initial_bz_code\": " << bz_code << "\n"
      << "  },\n"
      << "  \"mesh\": {\n"
      << "    \"x1min_code\": " << pm->mesh_size.x1min << ",\n"
      << "    \"x1max_code\": " << pm->mesh_size.x1max << ",\n"
      << "    \"x2min_code\": " << pm->mesh_size.x2min << ",\n"
      << "    \"x2max_code\": " << pm->mesh_size.x2max << ",\n"
      << "    \"x3min_code\": " << pm->mesh_size.x3min << ",\n"
      << "    \"x3max_code\": " << pm->mesh_size.x3max << ",\n"
      << "    \"nx1\": " << pm->mesh_size.nx1 << ",\n"
      << "    \"nx2\": " << pm->mesh_size.nx2 << ",\n"
      << "    \"nx3\": " << pm->mesh_size.nx3 << "\n"
      << "  }\n"
      << "}\n";
  out.close();
  if (!out) {
    std::stringstream msg;
    msg << "Failed while writing Toyouchi metadata file: " << temporary;
    ATHENA_ERROR(msg);
  }
  if (std::rename(temporary.c_str(), filename.c_str()) != 0) {
    std::stringstream msg;
    msg << "Could not finalize Toyouchi metadata file: " << filename;
    ATHENA_ERROR(msg);
  }
  metadata_written = true;
  std::cout << "Toyouchi runtime metadata: " << filename << std::endl;
}

// ===== magnetic field =====
Real bz_microgauss = 0.0;  // input value [microgauss]
Real bz_code = 0.0;        // Athena++へ実際に渡す無次元磁場(z方向)
Real Bunit = 0.0;          // magnetic-field unit [Gauss]

// ===== MMH profile fit =====

// Supplied coefficients, descending powers; physical AU / Msun / km s^-1.
const Real mass_coeff[] = {
  4.948071878736e-02, -1.267318191930e+00, 1.374036982238e+01,
  -8.163479027134e+01, 2.867189886433e+02, -5.954534679897e+02,
  6.808026360863e+02, -3.349246300414e+02};
const Real velocity_coeff[] = {
  -1.177335421707044e-01, 1.426676148040626e+00,
  -6.672185891204651e+00, 1.512617346113443e+01,
  -1.717522594895782e+01, 9.245470477699504e+00,
  -1.550902780998303e+00, 2.664339122790859e-01};

Real Polynomial(const Real *coeff, Real x, int derivative = 0) {
  Real value = 0.0;
  for (int i = 0; i <= 7-derivative; ++i) {
    Real c = coeff[i];
    for (int j = 0; j < derivative; ++j) c *= 7-i-j;
    value = value*x + c;
  }
  return value;
}

Real EvalMencLog10(Real x) { return Polynomial(mass_coeff, x); }
Real EvalMencDerivative(Real x) { return Polynomial(mass_coeff, x, 1); }

Real MencMMHMsun(Real r_code) {
  const Real r_au = r_code * Lunit / AU;
  if (r_au < fit_rmin_au) {
    const Real t = r_au / fit_rmin_au;
    return core_mass_msun*t*t*t*(core_a + t*t*(core_b + t*t*(core_c + core_d*t*t)));
  }
  return std::pow(10.0, EvalMencLog10(std::log10(
      std::min(r_au, fit_rmax_au))));
}

Real MencMMHCode(Real r_code) {
  return MencMMHMsun(r_code) * Msun / Munit;
}

Real RhoMMHCode(Real r_code) {
  const Real r_au = r_code * Lunit / AU;
  if (r_au < fit_rmin_au) {
    const Real t = r_au / fit_rmin_au;
    const Real rm = fit_rmin_au * AU;
    return core_mass_msun*Msun/(4.0*M_PI*rm*rm*rm*Rhounit)
        * (3.0*core_a + t*t*(5.0*core_b + t*t*(7.0*core_c + 9.0*core_d*t*t)));
  }
  // Copy the endpoint density outside the measured domain.
  const Real x = std::log10(std::min(r_au, fit_rmax_au));
  return Msun/(4.0*M_PI*AU*AU*AU*Rhounit)
      * std::pow(10.0, EvalMencLog10(x)-3.0*x)*EvalMencDerivative(x);
}

Real VrMMHCode(Real r_code) {
  const Real r_au = r_code * Lunit / AU;
  if (r_au < vr_taper_au) {
    const Real t = r_au / vr_taper_au;
    // Odd cubic: regular central vector field; value and slope match at taper.
    return taper_velocity_kms*t*(taper_a + taper_b*t*t)*1.0e5/Vunit;
  }
  const Real X = std::log10(MencMMHMsun(r_code));
  return Polynomial(velocity_coeff, X)*1.0e5/Vunit;
}

void InitProfiles() {
  const Real x = std::log10(fit_rmin_au);
  const Real slope = EvalMencDerivative(x);
  const Real curvature = Polynomial(mass_coeff, x, 2)/std::log(10.0);
  core_mass_msun = std::pow(10.0, EvalMencLog10(x));
  // Match M, density and density slope with an even central density.
  // Add a t^3*(1-t^2)^3 mass term to keep the central density positive.
  core_c = (curvature + slope*slope - 8.0*slope + 15.0)/8.0;
  core_b = (slope-3.0)/2.0 - 2.0*core_c;
  core_a = 1.0-core_b-core_c;
  const Real bump = std::max(static_cast<Real>(0.0), 0.25-core_a);
  core_a += bump;
  core_b -= 3.0*bump;
  core_c += 3.0*bump;
  core_d = -bump;
  const Real xt = std::log10(vr_taper_au);
  const Real X = EvalMencLog10(xt);
  const Real mt = std::pow(10.0, X);
  const Real mmax = std::pow(10.0, EvalMencLog10(std::log10(fit_rmax_au)));
  taper_velocity_kms = Polynomial(velocity_coeff, X);
  const Real velocity_slope = Polynomial(velocity_coeff, X, 1)
      * EvalMencDerivative(xt)/(std::log(10.0)*taper_velocity_kms);
  taper_a = (3.0-velocity_slope)/2.0;
  taper_b = (velocity_slope-1.0)/2.0;
  if (!(core_a > 0.0 && core_b >= 0.0 && core_c >= 0.0
        && 5.0*core_b + 7.0*core_c + 9.0*core_d >= 0.0
        && mt >= 1.0 && mmax <= 4679.64346356
        && std::isfinite(taper_velocity_kms) && taper_velocity_kms > 0.0
        && std::isfinite(taper_a) && taper_a > 0.0)) {
    std::stringstream msg;
    msg << "Invalid Wise19 core or velocity taper parameters.";
    ATHENA_ERROR(msg);
  }
}

} // namespace


void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Require explicit units so legacy inputs cannot silently use the ISM defaults.
  if (punit->unit_system != "custom") {
    std::stringstream msg;
    msg << "Toyouchi requires <units> unit_system=custom and mass_cgs, "
           "length_cgs, time_cgs.";
    ATHENA_ERROR(msg);
  }
  Munit = punit->code_mass_cgs;
  Lunit = punit->code_length_cgs;
  Tunit = punit->code_time_cgs;
  if (!(std::isfinite(Munit) && Munit > 0.0
        && std::isfinite(Lunit) && Lunit > 0.0
        && std::isfinite(Tunit) && Tunit > 0.0)) {
    std::stringstream msg;
    msg << "Toyouchi units must be finite and positive.";
    ATHENA_ERROR(msg);
  }
  Vunit = punit->code_velocity_cgs;
  Rhounit = punit->code_density_cgs;
  gconst = punit->grav_const_code;

  const Real sound_speed_kms = pin->GetReal("hydro", "sound_speed_kms");
  const Real epsilon_au = pin->GetReal("problem", "epsilon_au");
  if (!(std::isfinite(sound_speed_kms) && sound_speed_kms > 0.0
        && std::isfinite(epsilon_au) && epsilon_au >= 0.0)) {
    std::stringstream msg;
    msg << "sound_speed_kms must be finite and positive; "
           "epsilon_au must be finite and non-negative.";
    ATHENA_ERROR(msg);
  }
  // EOS constructors and ProblemGenerator read this code-unit value later.
  pin->SetReal("hydro", "iso_sound_speed", sound_speed_kms * 1.0e5 / Vunit);
  epsilon_soft = epsilon_au * AU / Lunit;

  if (Globals::my_rank == 0) {
    punit->PrintCodeUnits();
    std::cout << "Toyouchi converted parameters: G=" << gconst
              << " iso_sound_speed=" << sound_speed_kms * 1.0e5 / Vunit
              << " epsilon=" << epsilon_soft << std::endl;
  }


  // ===== 初期密度ノイズの設定 =====
  use_density_noise =
      pin->GetOrAddBoolean("problem", "use_density_noise", false);

  density_noise_percent =
      pin->GetOrAddReal("problem", "density_noise_percent", 0.0);

  density_noise_seed =
      pin->GetOrAddInteger("problem", "density_noise_seed", 42);

  if (!std::isfinite(density_noise_percent)
      || density_noise_percent < 0.0
      || density_noise_percent >= 100.0) {
    std::stringstream msg;
    msg << "density_noise_percent must satisfy 0 <= value < 100.";
    ATHENA_ERROR(msg);
  }

  if (density_noise_seed < 0) {
    std::stringstream msg;
    msg << "density_noise_seed must be non-negative.";
    ATHENA_ERROR(msg);
  }


  mstar_fixed = pin->GetOrAddReal("problem", "mstar_init", 1.0);
  vr_taper_au = pin->GetOrAddReal("problem", "vr_taper_au", 580.0);
  if (!(std::isfinite(mstar_fixed) && mstar_fixed >= 0.0
        && std::isfinite(vr_taper_au) && vr_taper_au >= fit_rmin_au
        && vr_taper_au <= fit_rmax_au)) {
    std::stringstream msg;
    msg << "Require finite mstar_init>=0 and vr_taper_au within mass fit range.";
    ATHENA_ERROR(msg);
  }
  InitProfiles();
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, StarMassHistory, "Mstar");
  if (Globals::my_rank == 0) {
    std::cout << "Fixed stellar mass=" << mstar_fixed
              << " code mass; Wise19 core radius=" << fit_rmin_au
              << " AU; velocity taper=" << vr_taper_au << " AU" << std::endl;
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
  use_radial_inflow = pin->GetOrAddBoolean("problem", "use_radial_inflow", true);
  use_rotation = pin->GetOrAddBoolean("problem","use_rotation",true);

  //磁場の強さを入力ファイルから読み込む(μG単位)
  bz_microgauss = pin->GetOrAddReal("problem", "bz_microgauss", 0.0);

  // CGS→コード単位の変換基準
  Bunit = punit->code_magneticfield_cgs;

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

  // Units が CGS 定数から換算した重力定数を自己重力にも使用
  if (SELF_GRAVITY_ENABLED) {
    SetGravitationalConstant(gconst);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // 中心星・NFWの動径加速度に加え、力の釣り合いを解析するための
  // 自己重力・圧力勾配・Lorentz加速度のCartesian 3成分を出力する。
  AllocateUserOutputVariables(11);
  SetUserOutputVariableName(0, "gr_star");  // 中心星の動径重力加速度
  SetUserOutputVariableName(1, "gr_nfw");   // NFW haloの動径重力加速度
  SetUserOutputVariableName(2, "gself_x");  // 自己重力加速度x成分
  SetUserOutputVariableName(3, "gself_y");  // 自己重力加速度y成分
  SetUserOutputVariableName(4, "gself_z");  // 自己重力加速度z成分
  SetUserOutputVariableName(5, "apres_x");  // 圧力勾配加速度x成分
  SetUserOutputVariableName(6, "apres_y");  // 圧力勾配加速度y成分
  SetUserOutputVariableName(7, "apres_z");  // 圧力勾配加速度z成分
  SetUserOutputVariableName(8, "amag_x");   // Lorentz加速度x成分
  SetUserOutputVariableName(9, "amag_y");   // Lorentz加速度y成分
  SetUserOutputVariableName(10, "amag_z");  // Lorentz加速度z成分
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

        Real rho_final = RhoMMHCode(r_sph);
        Real P_final;

        // Optional cell-wise density perturbation, including the central gas.
        if (use_density_noise && density_noise_percent > 0.0) {

         // この解像度レベルにおける領域全体のセル番号。
         // ローカルなi,j,kだけを使わず、ブロック位置も反映する。
         const std::uint64_t ix =
             static_cast<std::uint64_t>(loc.lx1)
             * static_cast<std::uint64_t>(block_size.nx1)
             + static_cast<std::uint64_t>(i - is);

         const std::uint64_t iy =
             static_cast<std::uint64_t>(loc.lx2)
             * static_cast<std::uint64_t>(block_size.nx2)
             + static_cast<std::uint64_t>(j - js);

         const std::uint64_t iz =
             static_cast<std::uint64_t>(loc.lx3)
             * static_cast<std::uint64_t>(block_size.nx3)
             + static_cast<std::uint64_t>(k - ks);

         const int refinement_level =
             loc.level - pmy_mesh->root_level;

         const Real noise =
             DensityNoiseValue(ix, iy, iz, refinement_level);

         rho_final *=
             1.0 + 0.01 * density_noise_percent * noise;
       }

        P_final = rho_final * cs * cs;

        phydro->u(IDN,k,j,i) = rho_final;

        // user output variablesを初期化する。
        // 力の診断量はUserWorkBeforeOutput()で出力直前に更新する。
        for (int n = 0; n < 11; ++n) {
          user_out_var(n,k,j,i) = 0.0;
        }

        // 追加 rotation velocity
        Real vx = 0.0, vy = 0.0, vz = 0.0;

        const Real r_cyl = std::sqrt(x*x + y*y);

        Real vphi_eff = 0.0;

        // rotation用Mencは円筒半径Rで評価する
        if (use_rotation && r_cyl > 0.0) {
          // Hold the endpoint rotation speed constant beyond the fitted radius.
          const Real r_rot = std::min(r_cyl, fit_rmax_au*AU/Lunit);
          const Real menc_rot_code = MencMMHCode(r_rot);

          vphi_eff =
              0.5 * std::sqrt(gconst * menc_rot_code / r_rot);
        }

        // 円筒回転をCartesian成分へ変換
        if (r_cyl > 0.0) {
          const Real inv_r_cyl = 1.0 / r_cyl;

          vx += -vphi_eff * y * inv_r_cyl;
          vy +=  vphi_eff * x * inv_r_cyl;
        }

        // 球対称のradial inflow
        if (use_radial_inflow && r_sph > 0.0) {
          const Real vr_code = VrMMHCode(r_sph);
          const Real inv_r_sph = 1.0 / r_sph;

          // fit値は正の大きさなので、マイナスを付けて中心向きにする。
          vx += -vr_code * x * inv_r_sph;
          vy += -vr_code * y * inv_r_sph;
          vz += -vr_code * z * inv_r_sph;
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
        if (r2 == 0.0) {
          pmb->user_out_var(0,k,j,i) = 0.0;
          pmb->user_out_var(1,k,j,i) = 0.0;
          continue;
        }
        Real r  = std::sqrt(r2);

	Real gx = 0.0;
	Real gy = 0.0;
	Real gz = 0.0;

        if (use_central_gravity) {
        const Real mstar = mstar_fixed;
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

        cons(IM1,k,j,i) += dt * prim(IDN,k,j,i) * gx;
        cons(IM2,k,j,i) += dt * prim(IDN,k,j,i) * gy;
        cons(IM3,k,j,i) += dt * prim(IDN,k,j,i) * gz;
      }
    }
  }
}

// VTK出力直前に、力の釣り合いを調べるための診断加速度を計算する。
// これらはcell-centered量の2次中心差分であり、Riemann solverが用いる
  // 離散fluxそのものではない。特に細分化境界では微分誤差に注意する。
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Athena++ changes into the run directory before this callback, so the metadata
  // file is placed beside the VTK output. The rank/process guard makes this a
  // one-time write even though the callback runs once per local MeshBlock.
  WriteToyouchiMetadata(pmy_mesh, pin,
                        pin->GetReal("hydro", "sound_speed_kms"),
                        pin->GetReal("problem", "epsilon_au"));

  AthenaArray<Real> &w = phydro->w;
  Coordinates *coord = pcoord;

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        // 隣接セル中心間の距離 [code length]
        const Real dx = coord->x1v(i+1) - coord->x1v(i-1);
        const Real dy = coord->x2v(j+1) - coord->x2v(j-1);
        const Real dz = coord->x3v(k+1) - coord->x3v(k-1);
        const Real rho = std::max(w(IDN,k,j,i), static_cast<Real>(1.0e-300));

        // 自己重力加速度 g_self=-grad(Phi_self) [code acceleration]
        Real gself_x = 0.0;
        Real gself_y = 0.0;
        Real gself_z = 0.0;
        if (SELF_GRAVITY_ENABLED) {
          AthenaArray<Real> &phi = pgrav->phi;
          gself_x = -(phi(k,j,i+1)-phi(k,j,i-1))/dx;
          gself_y = -(phi(k,j+1,i)-phi(k,j-1,i))/dy;
          gself_z = -(phi(k+1,j,i)-phi(k-1,j,i))/dz;
        }
        user_out_var(2,k,j,i) = gself_x;
        user_out_var(3,k,j,i) = gself_y;
        user_out_var(4,k,j,i) = gself_z;

        // 圧力勾配加速度 a_pres=-grad(P)/rho [code acceleration]
        // 断熱EOSではprimitive pressure、等温EOSではP=cs^2*rhoを使う。
#if NON_BAROTROPIC_EOS
        const Real pxm = w(IPR,k,j,i-1);
        const Real pxp = w(IPR,k,j,i+1);
        const Real pym = w(IPR,k,j-1,i);
        const Real pyp = w(IPR,k,j+1,i);
        const Real pzm = w(IPR,k-1,j,i);
        const Real pzp = w(IPR,k+1,j,i);
#else
        const Real pxm = cs2*w(IDN,k,j,i-1);
        const Real pxp = cs2*w(IDN,k,j,i+1);
        const Real pym = cs2*w(IDN,k,j-1,i);
        const Real pyp = cs2*w(IDN,k,j+1,i);
        const Real pzm = cs2*w(IDN,k-1,j,i);
        const Real pzp = cs2*w(IDN,k+1,j,i);
#endif
        user_out_var(5,k,j,i) = -((pxp-pxm)/dx)/rho;
        user_out_var(6,k,j,i) = -((pyp-pym)/dy)/rho;
        user_out_var(7,k,j,i) = -((pzp-pzm)/dz)/rho;

        // Lorentz加速度 a_mag=[curl(B) x B]/rho [code acceleration]
        // Athena++の磁場規格化では4*piはBの定義へ吸収されている。
        Real amag_x = 0.0;
        Real amag_y = 0.0;
        Real amag_z = 0.0;
#if MAGNETIC_FIELDS_ENABLED
        AthenaArray<Real> &bcc = pfield->bcc;
        const Real bx = bcc(IB1,k,j,i);
        const Real by = bcc(IB2,k,j,i);
        const Real bz = bcc(IB3,k,j,i);

        const Real dbz_dy = (bcc(IB3,k,j+1,i)-bcc(IB3,k,j-1,i))/dy;
        const Real dby_dz = (bcc(IB2,k+1,j,i)-bcc(IB2,k-1,j,i))/dz;
        const Real dbx_dz = (bcc(IB1,k+1,j,i)-bcc(IB1,k-1,j,i))/dz;
        const Real dbz_dx = (bcc(IB3,k,j,i+1)-bcc(IB3,k,j,i-1))/dx;
        const Real dby_dx = (bcc(IB2,k,j,i+1)-bcc(IB2,k,j,i-1))/dx;
        const Real dbx_dy = (bcc(IB1,k,j+1,i)-bcc(IB1,k,j-1,i))/dy;

        const Real jx = dbz_dy-dby_dz;
        const Real jy = dbx_dz-dbz_dx;
        const Real jz = dby_dx-dbx_dy;
        amag_x = (jy*bz-jz*by)/rho;
        amag_y = (jz*bx-jx*bz)/rho;
        amag_z = (jx*by-jy*bx)/rho;
#endif
        user_out_var(8,k,j,i) = amag_x;
        user_out_var(9,k,j,i) = amag_y;
        user_out_var(10,k,j,i) = amag_z;
      }
    }
  }
}
// History sums over blocks; report the fixed global mass only once.
Real StarMassHistory(MeshBlock *pmb, int iout) {
  return pmb->gid == 0 ? mstar_fixed : 0.0;
}

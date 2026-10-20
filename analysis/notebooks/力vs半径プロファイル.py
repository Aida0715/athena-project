# -*- coding: utf-8 -*-
# Converted from 力vs半径プロファイル.ipynb

# %% cell 1
#Toyouchi-test16の可視化に使った
#重力（中心星）＋遠心力＋圧力勾配をプロットするコード

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort

# ===============================
# パラメータ
# ===============================
vtk_dir = os.path.expanduser("~/athena-project/results/Toyouchi-test16")

xc, yc, zc = 0.0, 0.0, 0.0

r_max = 224.0
n_bins = 25

M_UNIT_CGS = 3.97694e33
L_UNIT_CGS = 1.495978707e15
T_UNIT_CGS = 3.15576e10
G_CGS = 6.67430e-8
G = G_CGS * M_UNIT_CGS * T_UNIT_CGS**2 / L_UNIT_CGS**3
AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
RHO_UNIT_CGS = M_UNIT_CGS / L_UNIT_CGS**3

FIT_RMIN_AU = 139.9720348
FIT_RMAX_AU = 146698.0663
MASS_COEFF = np.array([
    4.948071878736e-02, -1.267318191930e+00, 1.374036982238e+01,
    -8.163479027134e+01, 2.867189886433e+02, -5.954534679897e+02,
    6.808026360863e+02, -3.349246300414e+02,
])
_x_core = np.log10(FIT_RMIN_AU)
_dmass_coeff = np.polyder(MASS_COEFF)
_slope_core = np.polyval(_dmass_coeff, _x_core)
_curvature_core = np.polyval(np.polyder(MASS_COEFF, 2), _x_core) / np.log(10.0)
_core_mass_msun = 10.0**np.polyval(MASS_COEFF, _x_core)
_core_c = (_curvature_core + _slope_core**2 - 8.0*_slope_core + 15.0) / 8.0
_core_b = (_slope_core - 3.0) / 2.0 - 2.0*_core_c
_core_a = 1.0 - _core_b - _core_c
_bump = max(0.0, 0.25 - _core_a)
_core_a += _bump
_core_b -= 3.0*_bump
_core_c += 3.0*_bump
_core_d = -_bump


def wise19_profiles_code(r_code):
    r_au = np.asarray(r_code, dtype=float) * LENGTH_UNIT_AU
    clipped = np.clip(r_au, FIT_RMIN_AU, FIT_RMAX_AU)
    x = np.log10(clipped)
    mass_msun = 10.0**np.polyval(MASS_COEFF, x)
    rho_code = (
        MSUN_CGS / (4.0*np.pi*AU_CGS**3*RHO_UNIT_CGS)
        * 10.0**(np.polyval(MASS_COEFF, x) - 3.0*x)
        * np.polyval(_dmass_coeff, x)
    )
    inner = r_au < FIT_RMIN_AU
    if np.any(inner):
        t = np.maximum(r_au[inner], 0.0) / FIT_RMIN_AU
        mass_msun[inner] = _core_mass_msun*t**3 * (
            _core_a + t**2*(_core_b + t**2*(_core_c + _core_d*t**2))
        )
        rm = FIT_RMIN_AU * AU_CGS
        rho_code[inner] = (
            _core_mass_msun*MSUN_CGS/(4.0*np.pi*rm**3*RHO_UNIT_CGS)
            * (3.0*_core_a + t**2*(5.0*_core_b + t**2*(7.0*_core_c + 9.0*_core_d*t**2)))
        )
    return mass_msun*MSUN_CGS/M_UNIT_CGS, rho_code
Mstar = 1.0
cs = 17.346823586722763  # 8.22321048489829 km/s (T=1e4 K, mu=1.22)

z_cut = 2.0  # midplane

# 出力
output_dir1 = "./force_profiles1"
output_dir2 = "./force_profiles2"
os.makedirs(output_dir1, exist_ok=True)
os.makedirs(output_dir2, exist_ok=True)

# ===============================
# VTK整理
# ===============================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Toyouchi.block") and f.endswith(".vtk")
])

step_dict = {}
for f in vtk_files:
    step = f.split('.')[-2]
    step_dict.setdefault(step, []).append(f)

print(f"Found {len(step_dict)} timesteps")

# ===============================
# ビン
# ===============================
bin_edges = np.logspace(np.log10(1.5), np.log10(r_max), n_bins + 1)
bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])

# ===============================
# timestep loop
# ===============================
for step, files in step_dict.items():

    print(f"Processing step {step}")

    blocks = [pv.read(f) for f in files]
    grid = pv.MultiBlock(blocks).combine()

    centers = grid.cell_centers().points

    # ===== 物理量 =====
    rho = grid["rho"]
    vel = grid["vel"]

    vx = vel[:, 0]
    vy = vel[:, 1]
    vz = vel[:, 2]

    # ===== 座標 =====
    dx = centers[:, 0] - xc
    dy = centers[:, 1] - yc
    dz = centers[:, 2] - zc

    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    r_cyl = np.sqrt(dx*dx + dy*dy)

    # 単位ベクトル
    r_hat_x = dx / np.maximum(r, 1e-12)
    r_hat_y = dy / np.maximum(r, 1e-12)
    r_hat_z = dz / np.maximum(r, 1e-12)

    # midplane制限
    midplane = np.abs(dz) < z_cut

    # ===============================
    # 遠心力
    # ===============================
    vphi = (-dy * vx + dx * vy) / np.maximum(r_cyl, 1e-12)
    F_cent = vphi**2 / np.maximum(r_cyl, 1e-12)

    # ===============================
    # 重力
    # ===============================

    # 中心星
    F_star = -G * Mstar / np.maximum(r, 1e-12)**2

    # 自己重力（phiから計算）
    if "phi" in grid.array_names:

        phi = grid["phi"]

        # 勾配（cellベース）
        grad = grid.compute_derivative(scalars="phi", gradient=True)
        grad_phi = grad["gradient"]

        gx = -grad_phi[:, 0]
        gy = -grad_phi[:, 1]
        gz = -grad_phi[:, 2]

        F_self = gx*r_hat_x + gy*r_hat_y + gz*r_hat_z

    else:
        print("WARNING: phi not found → self-gravity ignored")
        F_self = 0.0

    F_grav = F_star + F_self

    # ===============================
    # 圧力勾配（簡易 radial bin 微分）
    # ===============================
    # 一旦ビン平均した後に微分する

    valid = (
        np.isfinite(rho)
        & (rho > 0)
        & (r > 0)
        & (r <= r_max)
        & midplane
    )

    r_valid = r[valid]
    rho_valid = rho[valid]

    bin_index = np.digitize(r_valid, bin_edges) - 1

    rho_shell = np.zeros(n_bins)

    for i in range(n_bins):
        mask = bin_index == i
        if np.any(mask):
            rho_shell[i] = rho_valid[mask].mean()
        else:
            rho_shell[i] = np.nan
            
    # --- ln rho ---
    lnrho = np.log(np.maximum(rho_shell, 1e-30))

    # --- d ln rho / dr ---
    dlnrho_dr = np.gradient(lnrho, bin_centers)

    # --- 圧力勾配力 ---
    F_pres_shell = - (cs**2)*dlnrho_dr

    # ===============================
    # 各力のビン平均
    # ===============================
    Fgrav_prof = np.zeros(n_bins)
    Fcent_prof = np.zeros(n_bins)

    for i in range(n_bins):
        mask = (bin_index == i)
        if np.any(mask):
            Fgrav_prof[i] = F_grav[valid][mask].mean()
            Fcent_prof[i] = F_cent[valid][mask].mean()
        else:
            Fgrav_prof[i] = np.nan
            Fcent_prof[i] = np.nan

    # 圧力はshell定義
    Fpres_prof = F_pres_shell

    # ===============================
    # 合力(シミュレーション)
    # ===============================
    F_tot = Fgrav_prof + Fcent_prof + Fpres_prof
    
    # ===============================
    # 各力(理論)
    # ===============================
    r_th = bin_centers
    # --- Wise et al. (2019) 密度と包有質量（Toyouchi.cpp と同一） ---
    M_enc, rho_th = wise19_profiles_code(r_th)

    # --- 回転 ---
    vphi_th = 0.5 * np.sqrt(G * M_enc / np.maximum(r_th, 1e-12))

    # --- 各力 ---
    Fself_th = -G * M_enc / np.maximum(r_th, 1e-12)**2   # ガス自己重力
    Fgrav_th = -G * Mstar / np.maximum(r_th, 1e-12)**2 + Fself_th   # 中心星＋ガス自己重力
    Fcent_th = vphi_th**2 / np.maximum(r_th, 1e-12)

    # 圧力（lnρ微分）
    lnrho_th = np.log(rho_th)
    dlnrho_dr_th = np.gradient(lnrho_th, r_th)
    Fpres_th = - (cs**2) * dlnrho_dr_th

    # --- 合力 ---
    Ftot_th = Fgrav_th + Fcent_th + Fpres_th

    # ===============================
    # プロット①：各力
    # ===============================
    plt.figure(figsize=(7,5))

    plt.plot(bin_centers, Fgrav_prof, label="Gravity", lw=1.2)
    plt.plot(bin_centers, Fcent_prof, label="Centrifugal", lw=1.2)
    plt.plot(bin_centers, Fpres_prof, label="Pressure", lw=1.2)

    plt.axhline(0, color="gray", lw=0.8)

    plt.xscale("log")
    plt.yscale("symlog", linthresh=1e-6)

    plt.xlabel("r (×100 AU)")
    plt.ylabel("Force (radial)")

    plt.title(f"Force components (step {step})")

    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)

    png1 = os.path.join(output_dir1,
                        f"force_components_{step}.png")

    plt.savefig(png1, dpi=150, bbox_inches="tight")
    plt.close()

    # ===============================
    # プロット②：合力のみ
    # ===============================
    plt.figure(figsize=(7,5))

    # --- シミュ ---
    plt.plot(bin_centers, F_tot,
             "k", lw=1.3, label="Simulation")

    # --- 理論 ---
    plt.plot(r_th, Ftot_th,
             "r--", lw=1.0, label="Theory")

    plt.axhline(0, color="gray", lw=0.8)

    plt.xscale("log")
    plt.yscale("symlog", linthresh=1e-6)

    plt.xlabel("r (×100 AU)")
    plt.ylabel("Force (radial)")

    plt.title(f"Total force (step {step})")

    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)

    png2 = os.path.join(output_dir2,
                        f"force_total_{step}.png")

    plt.savefig(png2, dpi=150, bbox_inches="tight")
    plt.close()

print("Done.")

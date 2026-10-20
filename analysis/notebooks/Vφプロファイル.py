# -*- coding: utf-8 -*-
# Converted from Vφプロファイル.ipynb

# %% cell 1
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort

vtk_dir = os.path.expanduser("~/athena-project/results/〇〇")

# ===============================
# パラメータ（cppと一致）
# ===============================
xc, yc, zc = 0.0, 0.0, 0.0

r_max = 244.0
n_bins = 25

M_UNIT_CGS = 3.97694e33
L_UNIT_CGS = 1.495978707e15
T_UNIT_CGS = 3.15576e10
G_CGS = 6.67430e-8
G = G_CGS * M_UNIT_CGS * T_UNIT_CGS**2 / L_UNIT_CGS**3
AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS

# Wise et al. (2019) 包有質量フィット（Toyouchi.cpp と同一）
FIT_RMIN_AU = 139.9720348
FIT_RMAX_AU = 146698.0663
MASS_COEFF = np.array([
    4.948071878736e-02, -1.267318191930e+00, 1.374036982238e+01,
    -8.163479027134e+01, 2.867189886433e+02, -5.954534679897e+02,
    6.808026360863e+02, -3.349246300414e+02,
])

_x_core = np.log10(FIT_RMIN_AU)
_slope_core = np.polyval(np.polyder(MASS_COEFF), _x_core)
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


def menc_wise19_code(r_code):
    r_au = np.asarray(r_code, dtype=float) * LENGTH_UNIT_AU
    clipped = np.clip(r_au, FIT_RMIN_AU, FIT_RMAX_AU)
    mass_msun = 10.0**np.polyval(MASS_COEFF, np.log10(clipped))
    inner = r_au < FIT_RMIN_AU
    if np.any(inner):
        t = np.maximum(r_au[inner], 0.0) / FIT_RMIN_AU
        mass_msun[inner] = _core_mass_msun*t**3 * (
            _core_a + t**2*(_core_b + t**2*(_core_c + _core_d*t**2))
        )
    return mass_msun * MSUN_CGS / M_UNIT_CGS

z_cut = 2.0

# 出力
output_dir_log = "./vphi_profiles_log"
output_dir_lin = "./vphi_profiles_linear"
os.makedirs(output_dir_log, exist_ok=True)
os.makedirs(output_dir_lin, exist_ok=True)

# ======================================
# VTK ファイル整理
# ======================================
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

# ======================================
# ビン
# ======================================
bin_edges = np.logspace(np.log10(1.5), np.log10(r_max), n_bins + 1)
bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])

# ======================================
# 理論 Vphi（cppと一致）
# ======================================
def theoretical_vphi(r):
    M_gas = menc_wise19_code(r)
    return 0.5 * np.sqrt(G * M_gas / np.maximum(r, 1e-12))

# ======================================
# timestep loop
# ======================================
for step, files in step_dict.items():

    print(f"Processing step {step}")

    blocks = [pv.read(f) for f in files]
    grid = pv.MultiBlock(blocks).combine()

    centers = grid.cell_centers().points

    # 速度
    vel = grid["vel"]
    vx = vel[:, 0]
    vy = vel[:, 1]

    # 座標
    dx = centers[:, 0] - xc
    dy = centers[:, 1] - yc
    dz = centers[:, 2] - zc

    r_cyl = np.sqrt(dx*dx + dy*dy)

    # midplane抽出
    midplane = np.abs(dz) < z_cut

    # Vphi
    vphi = (-dy * vx + dx * vy) / np.maximum(r_cyl, 1e-12)

    valid = (
        np.isfinite(vphi)
        & (r_cyl > 0)
        & (r_cyl <= r_max)
        & midplane
    )

    r = r_cyl[valid]
    vphi = vphi[valid]

    # ======================================
    # ビン平均
    # ======================================
    prof = np.zeros(n_bins)

    bin_index = np.digitize(r, bin_edges) - 1

    for i in range(n_bins):
        mask = bin_index == i

        if np.any(mask):
            prof[i] = vphi[mask].mean()
        else:
            prof[i] = np.nan

    # ======================================
    # 理論
    # ======================================
    theory = theoretical_vphi(bin_centers)

    # ======================================
    # ① LOGプロット
    # ======================================
    plt.figure(figsize=(6,4))

    plt.plot(bin_centers, prof, lw=1.5, label="Simulation")
    plt.plot(bin_centers, theory, "--", lw=1.2,
             label=r"Theory (gas-only 0.5 $v_{\rm kep}$)")

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("r (×100 AU)")
    plt.ylabel(r"$v_\phi$")
    plt.title(f"Rotation profile LOG (step {step})")

    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

    png = os.path.join(output_dir_log, f"vphi_log_{step}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close()

    # ======================================
    # ② LINEARプロット
    # ======================================
    plt.figure(figsize=(6,4))

    plt.plot(bin_centers, prof, lw=1.5, label="Simulation")
    plt.plot(bin_centers, theory, "--", lw=1.2,
             label=r"Theory (gas-only 0.5 $v_{\rm kep}$)")

    plt.xscale("linear")
    plt.yscale("linear")

    plt.xlabel("r (×100 AU)")
    plt.ylabel(r"$v_\phi$")
    plt.title(f"Rotation profile LINEAR (step {step})")

    plt.grid(True, ls="--", alpha=0.5)
    plt.legend()

    png = os.path.join(output_dir_lin, f"vphi_linear_{step}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close()

print("Done.")

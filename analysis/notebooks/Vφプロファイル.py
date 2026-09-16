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

G = 1.0

racc = 0.44
r_cut = max(racc, 1e-12)

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

    M_gas = 4.0 * np.pi * 0.11 * r**1.25 / 1.25

    vphi = 0.5 * np.sqrt(G * M_gas / np.maximum(r, r_cut))

    mask = r < r_cut
    vphi[mask] *= (r[mask] / r_cut)**2

    return vphi

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

    plt.xlabel("r(×450AU)")
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

    plt.xlabel("r(×450AU)")
    plt.ylabel(r"$v_\phi$")
    plt.title(f"Rotation profile LINEAR (step {step})")

    plt.grid(True, ls="--", alpha=0.5)
    plt.legend()

    png = os.path.join(output_dir_lin, f"vphi_linear_{step}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close()

print("Done.")


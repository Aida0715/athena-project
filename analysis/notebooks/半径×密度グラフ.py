# -*- coding: utf-8 -*-
# Converted from 半径×密度グラフ.ipynb

# %% cell 1
#Toyouchi-test4を解析するときに使った
#Toyouchi+23再現実験では今後これを使う

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort

vtk_dir = os.path.expanduser(
    "~/athena-project/results/〇〇"
)

output_dir1 = "./radial_density_profiles"
output_dir2 = "./radial_density_scatter"
os.makedirs(output_dir1, exist_ok=True)
os.makedirs(output_dir2, exist_ok=True)

# =====================================
# 設定
# =====================================
xc, yc, zc = 0.0, 0.0, 0.0
r_max = 20000.0
n_bins = 30

# =====================================
# VTKファイル整理
# =====================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Toyouchi.block")
    and f.endswith(".vtk")
])

step_dict = {}

for f in vtk_files:
    step = os.path.basename(f).split(".")[-2]
    step_dict.setdefault(step, []).append(f)

print("Found", len(step_dict), "timesteps")

# =====================================
# 理論密度
# =====================================
def toyouchi_profile(r):
    return r**(-1.75)

# =====================================
# timestep loop
# =====================================
for step, files in step_dict.items():

    print("\n=========================")
    print("step =", step)
    print("=========================")

    all_r = []
    all_rho = []

    # -------------------------------
    # blockごとに読む
    # -------------------------------
    for f in files:

        grid = pv.read(f)

        centers = grid.cell_centers().points
        rho = grid["rho"]

        x = centers[:,0] - xc
        y = centers[:,1] - yc
        z = centers[:,2] - zc

        r = np.sqrt(x*x + y*y + z*z)

        valid = (
            np.isfinite(rho)
            & (rho > 0.0)
            & np.isfinite(r)
            & (r <= r_max)
        )

        all_r.append(r[valid])
        all_rho.append(rho[valid])

    # 全block結合
    r = np.concatenate(all_r)
    rho = np.concatenate(all_rho)

    print("Ncell =", len(r))
    print("min r =", np.min(r))
    print("max r =", np.max(r))

    # ----------------------------------
    # 散布図を保存
    # ----------------------------------
    plt.figure(figsize=(6,4))

    plt.scatter(r, rho, s=0.3)

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("r")
    plt.ylabel("rho")
    plt.title(f"Scatter step {step}")

    plt.grid(alpha=0.3)

    plt.savefig(
        os.path.join(
            output_dir2,
            f"scatter_{step}.png"
        ),
        dpi=150,
        bbox_inches="tight"
    )

    plt.close()

    # ----------------------------------
    # ビンを自動決定
    # ----------------------------------
    rmin = np.min(r)

    bin_edges = np.logspace(
        np.log10(rmin),
        np.log10(r_max),
        n_bins+1
    )

    bin_centers = np.sqrt(
        bin_edges[:-1]*bin_edges[1:]
    )

    rho_shell = np.zeros(n_bins)

    counts = np.zeros(n_bins)

    ind = np.digitize(r, bin_edges)-1

    for i in range(n_bins):

        mask = (ind == i)

        counts[i] = np.sum(mask)

        if counts[i] > 0:
            rho_shell[i] = np.mean(rho[mask])
        else:
            rho_shell[i] = np.nan

    print("cells/bin =")
    print(counts)

    # ----------------------------------
    # 理論曲線
    # ----------------------------------
    theory = toyouchi_profile(bin_centers)

    scale = np.nanmedian(rho_shell/theory)
    theory *= scale

    # ----------------------------------
    # プロット
    # ----------------------------------
    plt.figure(figsize=(7,5))

    plt.plot(
        bin_centers,
        rho_shell,
        lw=1.5,
        label="Simulation"
    )

    plt.plot(
        bin_centers,
        theory,
        "--",
        lw=1.0,
        label=r"$\rho \propto r^{-1.75}$"
    )

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("r(×450AU)")
    plt.ylabel("Density")

    plt.title(
        f"Spherical density profile step {step}"
    )

    plt.grid(True, which="both",
             ls="--", alpha=0.4)

    plt.legend()

    plt.savefig(
        os.path.join(
            output_dir1,
            f"profile_{step}.png"
        ),
        dpi=200,
        bbox_inches="tight"
    )

    plt.close()

print("\nDone.")


# -*- coding: utf-8 -*-
# Converted from 半径×圧力グラフ.ipynb

# %% cell 1
#半径vs圧力グラフ　出力されるconsデータから作っている（jeans_3d-test14を可視化したときに使った）

import os
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import natsort

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")

# 球・計算領域
xc, yc, zc = 6.0, 6.0, 6.0     # 球中心
r_max = 6.0                   # 評価する最大半径
n_bins = 200                  # 球殻数（滑らかさはここ）

gamma = 1.20                  # ★ Athena++ input と必ず一致させる

# 出力
output_dir = "./radial_pressure_profiles"
os.makedirs(output_dir, exist_ok=True)

# ======================================
# VTK ファイルを timestep ごとに整理
# ======================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

step_dict = {}
for f in vtk_files:
    step = f.split('.')[-2]   # 00000 など
    step_dict.setdefault(step, []).append(f)

print(f"Found {len(step_dict)} timesteps")

# ======================================
# 球殻ビン定義（r=0 を含む）
# ======================================
bin_edges = np.linspace(0.0, r_max, n_bins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# ======================================
# 各 timestep の処理
# ======================================
for step, files in step_dict.items():

    # --- 全 block を結合 ---
    blocks = [pv.read(f) for f in files]
    grid = pv.MultiBlock(blocks).combine()

    print(f"[step {step}] bounds = {grid.bounds}")

    # --- セル中心 ---
    centers = grid.cell_centers().points
    x = centers[:, 0]
    y = centers[:, 1]
    z = centers[:, 2]

    # --- 半径 ---
    dx = x - xc
    dy = y - yc
    dz = z - zc
    r = np.sqrt(dx*dx + dy*dy + dz*dz)

    # --- 保存量 ---
    rho  = grid.cell_data["dens"]
    Etot = grid.cell_data["Etot"]
    phi  = grid.cell_data["phi"]
    mom  = grid.cell_data["mom"]   # (N,3)

    # --- 速度 ---
    vx = mom[:, 0] / rho
    vy = mom[:, 1] / rho
    vz = mom[:, 2] / rho
    v2 = vx*vx + vy*vy + vz*vz

    # --- 圧力（最重要） ---
    P = (gamma - 1.0) * (Etot - 0.5 * rho * v2 - rho * phi)

    # --- 有効セルのみ ---
    valid = (
        np.isfinite(P) &
        np.isfinite(rho) &
        (rho > 0.0) &
        (P > 0.0) &
        (r <= r_max)
    )

    r = r[valid]
    P = P[valid]

    # ======================================
    # 球殻平均
    # ======================================
    P_shell = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=int)

    bin_index = np.digitize(r, bin_edges) - 1

    for i in range(n_bins):
        mask = bin_index == i
        if np.any(mask):
            P_shell[i] = P[mask].mean()
            counts[i] = np.count_nonzero(mask)
        else:
            P_shell[i] = np.nan

    print(
        f"  center bin: r~{bin_centers[0]:.3e}, "
        f"P={P_shell[0]:.6e}, N={counts[0]}"
    )

    # ======================================
    # プロット（log スケール）
    # ======================================
    plt.figure(figsize=(6, 4))

    plt.plot(
        bin_centers,
        P_shell,
        lw=0.4,          # 線は細め
        color="C1"
    )

    plt.scatter(
        bin_centers,
        P_shell,
        s=2,
        color="C1"
    )

    plt.yscale("log")
    plt.xlabel("r")
    plt.ylabel("Pressure")
    plt.xlim(0, r_max)
    plt.ylim(0.000001,2.0)
    plt.title(f"Spherically averaged pressure (step {step})")
    plt.grid(True, which="both", alpha=0.3)

    png = os.path.join(output_dir, f"radial_pressure_shell_{step}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved {png}")

print("Done.")


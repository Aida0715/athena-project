# -*- coding: utf-8 -*-
# Converted from 過去作品.ipynb

# %% cell 1
#もう使わないかもしれないが念の為とっておく

# %% cell 2
#密度vs半径グラフ（jeans_3d-test14を可視化するときに使った）
#レーンエムデン用

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")

# 球・計算領域
xc, yc, zc = 6.0, 6.0, 6.0   # 球中心
r_max = 6.0                 # 球半径（Rsphere=5 だが少し余裕）
n_bins = 200                # 球殻数（滑らかさはここで決まる）

# 出力
output_dir = "./radial_density_profiles"
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
# 球殻ビン定義（r=0 を必ず含む）
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

    # --- セル中心と密度 ---
    centers = grid.cell_centers().points
    rho = grid["dens"]

    # --- 半径 ---
    dx = centers[:, 0] - xc
    dy = centers[:, 1] - yc
    dz = centers[:, 2] - zc
    r = np.sqrt(dx*dx + dy*dy + dz*dz)

    # --- 有効領域のみ ---
    valid = np.isfinite(rho) & (rho > 0) & (r <= r_max)
    r = r[valid]
    rho = rho[valid]

    # ======================================
    # 球殻平均
    # ======================================
    rho_shell = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=int)

    bin_index = np.digitize(r, bin_edges) - 1

    for i in range(n_bins):
        mask = bin_index == i
        if np.any(mask):
            rho_shell[i] = rho[mask].mean()
            counts[i] = np.count_nonzero(mask)
        else:
            rho_shell[i] = np.nan

    # --- 中心ビンの確認 ---
    print(f"  center bin: r~{bin_centers[0]:.3e}, "
          f"rho={rho_shell[0]:.6f}, N={counts[0]}")

    # ======================================
    # プロット
    # ======================================
    plt.figure(figsize=(6, 4))
    plt.plot(bin_centers, rho_shell, lw=2, color="C0")
    plt.scatter(bin_centers, rho_shell, s=10, color="C0")

    plt.xlabel("r")
    plt.ylabel("Density")
    plt.xlim(0, 2.0)
    plt.ylim(-0.1,1.0)
    plt.title(f"Spherically averaged density (step {step})")
    plt.grid(True)

    png = os.path.join(output_dir, f"radial_density_shell_{step}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved {png}")

print("Done.")

# %% cell 3
#半径vs密度　縦軸logスケール（jeans_3d-test14を可視化するときに使った）
#レーンエムデン用

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
import matplotlib.animation as animation
from PIL import Image

vtk_dir = os.path.expanduser("~/athenapp/results/jeans_3d-test27")

# 球・計算領域
xc, yc, zc = 6.0, 6.0, 6.0   # 球中心
r_max = 6.0                 # 球半径（Rsphere=5 だが少し余裕）
n_bins = 200                # 球殻数（滑らかさはここで決まる）

# 出力
output_dir = "./radial_density_profiles_log"
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
# 球殻ビン定義（r=0 を必ず含む）
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

    # --- セル中心と密度 ---
    centers = grid.cell_centers().points
    rho = grid["dens"]

    # --- 半径 ---
    dx = centers[:, 0] - xc
    dy = centers[:, 1] - yc
    dz = centers[:, 2] - zc
    r = np.sqrt(dx*dx + dy*dy + dz*dz)

    # --- 有効領域のみ ---
    valid = np.isfinite(rho) & (rho > 0) & (r <= r_max)
    r = r[valid]
    rho = rho[valid]

    # ======================================
    # 球殻平均
    # ======================================
    rho_shell = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=int)

    bin_index = np.digitize(r, bin_edges) - 1

    for i in range(n_bins):
        mask = bin_index == i
        if np.any(mask):
            rho_shell[i] = rho[mask].mean()
            counts[i] = np.count_nonzero(mask)
        else:
            rho_shell[i] = np.nan

    # --- 中心ビンの確認 ---
    print(f"  center bin: r~{bin_centers[0]:.3e}, "
          f"rho={rho_shell[0]:.6f}, N={counts[0]}")

    # ======================================
    # プロット
    # ======================================
    rho_floor = 1e-8   # ← 見たい精度に応じて調整（1e-8 などでもOK）

    rho_plot = rho_shell.copy()
    rho_plot[(rho_plot <= rho_floor) | ~np.isfinite(rho_plot)] = np.nan

    plt.figure(figsize=(6, 4))
    plt.plot(bin_centers, rho_plot, lw=0.6, color="C0")
    plt.scatter(bin_centers, rho_plot, s=3, color="C0", alpha=0.4)

    plt.xlabel("r")
    plt.ylabel("Density")
    plt.xlim(0, 6.0)
    plt.yscale("log")
    plt.ylim(rho_floor, 1.5)

    plt.title(f"Spherically averaged density (File {step})")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    png = os.path.join(output_dir, f"radial_density_log_{step}.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close()
    
print("Done.")

# %% cell 4
#ｘ−ｙ平面の密度map作成用コード（jeans_3d-test14の可視化に使ったやつ）

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from PIL import Image
import matplotlib.animation as animation

#密度画像（x-y平面密度）

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_density_maps"
os.makedirs(output_dir, exist_ok=True)

# 計算領域・球パラメータ
Lbox = 12.0
xc, yc = 6.0, 6.0
z_slice = 6.0          # 箱中心
Rsphere = 5.0

# 全 block の VTK ファイルを自動取得（block0〜blockXX）
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

# ブロック番号ごとにまとめる（x-y平面密度）

block_dict = {}
for f in vtk_files:
    basename = os.path.basename(f)
    block_id = int(basename.split('.')[1].replace('block',''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内のみで vmin / vmax を決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):
    points_all = []
    dens_all = []

    for f in step_files:
        grid = pv.read(f)
        xy_slice = grid.slice(normal='z', origin=(0,0,z_slice))
        centers = xy_slice.cell_centers().points
        dens = xy_slice['dens']

        points_all.append(centers)
        dens_all.append(dens)

    points_all = np.vstack(points_all)
    dens_all = np.hstack(dens_all)

    r = np.sqrt((points_all[:,0]-xc)**2 + (points_all[:,1]-yc)**2)
    mask_inside = r <= Rsphere
    dens_inside_all.append(dens_all[mask_inside])

dens_inside_all = np.concatenate(dens_inside_all)

# 極端値を除外
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"Density scale (inside sphere): vmin={vmin}, vmax={vmax}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    points_all = []
    dens_all = []

    for f in step_files:
        grid = pv.read(f)
        xy_slice = grid.slice(normal='z', origin=(0,0,z_slice))
        centers = xy_slice.cell_centers().points
        dens = xy_slice['dens']

        points_all.append(centers)
        dens_all.append(dens)

    points_all = np.vstack(points_all)
    dens_all = np.hstack(dens_all)

    xs_unique = np.unique(points_all[:,0])
    ys_unique = np.unique(points_all[:,1])
    nx, ny = len(xs_unique), len(ys_unique)

    dens_2d = np.full((ny, nx), np.nan)

    for i, x in enumerate(xs_unique):
        for j, y in enumerate(ys_unique):
            mask = (np.isclose(points_all[:,0], x)) & \
                   (np.isclose(points_all[:,1], y))
            if np.any(mask):
                dens_2d[j, i] = dens_all[mask][0]

    # 球内・球外マスク
    X, Y = np.meshgrid(xs_unique, ys_unique)
    r = np.sqrt((X-xc)**2 + (Y-yc)**2)

    dens_in  = np.ma.masked_where(r > Rsphere, dens_2d)
    dens_out = np.ma.masked_where(r <= Rsphere, dens_2d)

    # --------------------------------------------------------
    # 描画
    # --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal')

    # 球外：背景（単色）
    ax.pcolormesh(xs_unique, ys_unique, dens_out,
                  cmap='Greys', vmin=0, vmax=1)

    # 球内：密度勾配（対数）
    im = ax.pcolormesh(xs_unique, ys_unique, dens_in,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            yc + Rsphere*np.sin(theta),
            'c--', lw=1.2)

    ax.set_xlim(0, Lbox)
    ax.set_ylim(0, Lbox)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    filename = os.path.basename(step_files[0])
    num = filename.split('.')[-2]
    ax.set_title(f"file={filename}")

    fig.colorbar(im, ax=ax, label="Density (inside sphere)")

    png_file = os.path.join(output_dir, f"dens_{num}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

# %% cell 5
#ｘ−ｙ平面の密度map作成用コード（jeans_3d-test16の可視化に使ったやつ）
#球外へのmassの漏れ出しが視覚的にわかりやすい

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from PIL import Image

# ============================================================
# 設定
# ============================================================
vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_density_maps"
os.makedirs(output_dir, exist_ok=True)

Lbox = 8.0
xc, yc = 4.0, 4.0
z_slice = 4.0
Rsphere = 1.25
rho_bg = 1e-3   # 入力ファイルと合わせる

# ============================================================
# VTK ファイル整理（block ごと）
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    block_id = int(os.path.basename(f).split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内密度から vmin / vmax を決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        slc  = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
        pts_all.append(slc.cell_centers().points)
        dens_all.append(slc['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    xs = np.unique(pts_all[:,0])
    ys = np.unique(pts_all[:,1])

    dens_2d = np.full((len(ys), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(xs == p[0])[0][0]
        iy = np.where(ys == p[1])[0][0]
        dens_2d[iy, ix] = d

    X, Y = np.meshgrid(xs, ys)
    r2d  = np.sqrt((X-xc)**2 + (Y-yc)**2)

    dens_inside_all.append(dens_2d[r2d <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)

vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"[INFO] Density scale: vmin={vmin:.3e}, vmax={vmax:.3e}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        slc  = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
        pts_all.append(slc.cell_centers().points)
        dens_all.append(slc['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    xs = np.unique(pts_all[:,0])
    ys = np.unique(pts_all[:,1])

    dens_2d = np.full((len(ys), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(xs == p[0])[0][0]
        iy = np.where(ys == p[1])[0][0]
        dens_2d[iy, ix] = d

    X, Y = np.meshgrid(xs, ys)
    r2d  = np.sqrt((X-xc)**2 + (Y-yc)**2)

    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal')

    # ===============================
    # ① 背景：全領域を linear 表示
    # ===============================
    ax.pcolormesh(xs, ys, dens_2d,
                  cmap='Greys',
                  vmin=rho_bg,
                  vmax=vmin)

    # ===============================
    # ② 球内：LogNorm で上書き
    # ===============================
    mask_sphere = r2d <= Rsphere
    dens_sphere = np.ma.masked_where(~mask_sphere, dens_2d)

    im = ax.pcolormesh(xs, ys, dens_sphere,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            yc + Rsphere*np.sin(theta),
            'c--', lw=1.2)

    ax.set_xlim(0, Lbox)
    ax.set_ylim(0, Lbox)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    fname = os.path.basename(step_files[0])
    num   = fname.split('.')[-2]
    ax.set_title(f"z={z_slice}, γ=2.0, β=0.02")

    fig.colorbar(im, ax=ax, label="Density")

    png = os.path.join(output_dir, f"dens_{num}.png")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {png}")

# %% cell 6
#可視化のために球外に漏れ出た密度をいじって色付けしている。球外の色付けは実際に漏れ出た密度の絶対値ではないので注意
#スパイラルアームなどが立つとき、より見やすい可視化コード（jeans_3d-test24を可視化したときに使った）
#x-y平面用

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
import matplotlib.animation as animation
from PIL import Image
from matplotlib.colors import LogNorm, ListedColormap

# ============================================================
# 設定
# ============================================================
vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_density_maps"
os.makedirs(output_dir, exist_ok=True)

Lbox = 12.0
xc, yc = 6.0, 6.0
z_slice = 6.0
Rsphere = 5.0

rho_bg   = 1e-15            # 初期背景密度
rho_leak = 5.0 * rho_bg    # 漏れ出た質量の判定閾値（固定）

# 漏れ出た質量用：単色・強調ブルー
leak_cmap = ListedColormap(["deepskyblue"])

# ============================================================
# VTK ファイル整理（block ごと）
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    block_id = int(os.path.basename(f).split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内密度から vmin / vmax を決定（全ステップ共通）
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        slc  = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
        pts_all.append(slc.cell_centers().points)
        dens_all.append(slc['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    xs = np.unique(pts_all[:, 0])
    ys = np.unique(pts_all[:, 1])

    dens_2d = np.full((len(ys), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(xs == p[0])[0][0]
        iy = np.where(ys == p[1])[0][0]
        dens_2d[iy, ix] = d

    X, Y = np.meshgrid(xs, ys)
    r2d  = np.sqrt((X - xc)**2 + (Y - yc)**2)

    dens_inside_all.append(dens_2d[r2d <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)

vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"[INFO] Sphere density scale: vmin={vmin:.3e}, vmax={vmax:.3e}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        slc  = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
        pts_all.append(slc.cell_centers().points)
        dens_all.append(slc['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    xs = np.unique(pts_all[:, 0])
    ys = np.unique(pts_all[:, 1])

    dens_2d = np.full((len(ys), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(xs == p[0])[0][0]
        iy = np.where(ys == p[1])[0][0]
        dens_2d[iy, ix] = d

    X, Y = np.meshgrid(xs, ys)
    r2d  = np.sqrt((X - xc)**2 + (Y - yc)**2)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')

    # ========================================================
    # ① 球外背景（薄グレー）
    # ========================================================
    mask_bg = (r2d > Rsphere) & (dens_2d <= rho_leak)
    dens_bg = np.ma.masked_where(~mask_bg, dens_2d)

    ax.pcolormesh(xs, ys, dens_bg,
                  cmap='Greys',
                  vmin=rho_bg,
                  vmax=rho_leak,
                  alpha=0.35)

    # ========================================================
    # ② 漏れ出た質量（単色・濃いブルー）
    # ========================================================
    mask_leak = (r2d > Rsphere) & (dens_2d > rho_leak)
    dens_leak = np.ma.masked_where(~mask_leak, dens_2d)

    ax.pcolormesh(xs, ys, dens_leak,
                  cmap=leak_cmap,
                  alpha=0.45)

    # 境界を等高線で強調（スパイラル検出用）
    ax.contour(X, Y, dens_2d,
               levels=[rho_leak],
               colors='dodgerblue',
               linewidths=1.5)

    # ========================================================
    # ③ 球内（LogNorm）
    # ========================================================
    mask_sphere = r2d <= Rsphere
    dens_sphere = np.ma.masked_where(~mask_sphere, dens_2d)

    im = ax.pcolormesh(xs, ys, dens_sphere,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            yc + Rsphere*np.sin(theta),
            'c--', lw=1.2)

    ax.set_xlim(0, Lbox)
    ax.set_ylim(0, Lbox)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    fname = os.path.basename(step_files[0])
    num   = fname.split('.')[-2]
    ax.set_title(f"z={z_slice}, γ=1.2, β=0.02")

    fig.colorbar(im, ax=ax, label="Density (sphere)")

    png = os.path.join(output_dir, f"dens_{num}.png")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {png}")

# %% cell 7
#可視化のために球外に漏れ出た密度をいじって色付けしている。球外の色付けは実際に漏れ出た密度の絶対値ではないので注意
#スパイラルアームなどが立つとき、より見やすい可視化コード（jeans_3d-test24を可視化したときに使った）
#x-z平面用



import os
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import natsort
from matplotlib.colors import LogNorm, ListedColormap

# ============================================================
# 設定
# ============================================================
vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xz_density_maps"
os.makedirs(output_dir, exist_ok=True)

Lbox = 12.0
xc, yc, zc = 6.0, 6.0, 6.0
y_slice = 6.0          # ★ y でスライス
Rsphere = 5.0

rho_bg   = 1e-15
rho_leak = 5.0 * rho_bg

# 漏れ出た質量用：単色ブルー
leak_cmap = ListedColormap(["deepskyblue"])

# ============================================================
# VTK ファイル整理（block ごと）
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    block_id = int(os.path.basename(f).split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内密度から vmin / vmax を決定（全ステップ共通）
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        slc  = grid.slice(normal='y', origin=(0.0, y_slice, 0.0))
        pts_all.append(slc.cell_centers().points)
        dens_all.append(slc['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    xs = np.unique(pts_all[:, 0])
    zs = np.unique(pts_all[:, 2])

    dens_2d = np.full((len(zs), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(xs == p[0])[0][0]
        iz = np.where(zs == p[2])[0][0]
        dens_2d[iz, ix] = d

    X, Z = np.meshgrid(xs, zs)
    r2d  = np.sqrt((X - xc)**2 + (Z - zc)**2)

    dens_inside_all.append(dens_2d[r2d <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"[INFO] Sphere density scale: vmin={vmin:.3e}, vmax={vmax:.3e}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        slc  = grid.slice(normal='y', origin=(0.0, y_slice, 0.0))
        pts_all.append(slc.cell_centers().points)
        dens_all.append(slc['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    xs = np.unique(pts_all[:, 0])
    zs = np.unique(pts_all[:, 2])

    dens_2d = np.full((len(zs), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(xs == p[0])[0][0]
        iz = np.where(zs == p[2])[0][0]
        dens_2d[iz, ix] = d

    X, Z = np.meshgrid(xs, zs)
    r2d  = np.sqrt((X - xc)**2 + (Z - zc)**2)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')

    # ========================================================
    # ① 球外背景（薄グレー）
    # ========================================================
    mask_bg = (r2d > Rsphere) & (dens_2d <= rho_leak)
    dens_bg = np.ma.masked_where(~mask_bg, dens_2d)

    ax.pcolormesh(xs, zs, dens_bg,
                  cmap='Greys',
                  vmin=rho_bg,
                  vmax=rho_leak,
                  alpha=0.35)

    # ========================================================
    # ② 漏れ出た質量（ブルー）
    # ========================================================
    mask_leak = (r2d > Rsphere) & (dens_2d > rho_leak)
    dens_leak = np.ma.masked_where(~mask_leak, dens_2d)

    ax.pcolormesh(xs, zs, dens_leak,
                  cmap=leak_cmap,
                  alpha=0.45)

    ax.contour(X, Z, dens_2d,
               levels=[rho_leak],
               colors='dodgerblue',
               linewidths=1.5)

    # ========================================================
    # ③ 球内（LogNorm）
    # ========================================================
    mask_sphere = r2d <= Rsphere
    dens_sphere = np.ma.masked_where(~mask_sphere, dens_2d)

    im = ax.pcolormesh(xs, zs, dens_sphere,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            zc + Rsphere*np.sin(theta),
            'c--', lw=1.2)

    ax.set_xlim(0, Lbox)
    ax.set_ylim(0, Lbox)
    ax.set_xlabel("X")
    ax.set_ylabel("Z")

    fname = os.path.basename(step_files[0])
    num   = fname.split('.')[-2]
    ax.set_title(f"y={y_slice}, γ=1.2, β=0.02")

    fig.colorbar(im, ax=ax, label="Density (sphere)")

    png = os.path.join(output_dir, f"dens_xz_{num}.png")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {png}")

# %% cell 8
#Toyouchi-test4を可視化するとき使った
#Toyouchi+23再現ではこれを使うと良い
#x-y平面用

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm

# ============================================================
# 設定
# ============================================================
vtk_dir = os.path.expanduser("~/athenapp/results/your_result_directory")  # ★要修正★
output_dir = "./xy_density_maps"
os.makedirs(output_dir, exist_ok=True)

# 計算領域の設定（一辺448の立方体、中心が原点）
Lbox = 448.0
# 座標範囲: -224 から 224 まで
x_min, x_max = -Lbox/2, Lbox/2  # -224 から 224
y_min, y_max = -Lbox/2, Lbox/2
z_slice = 0.0  # 真ん中で切断（原点を通る面）

# 中心は原点
xc, yc = 0.0, 0.0

# プロットの解像度を上げるための設定
dpi = 300  # 高解像度で保存
figsize = (12, 10)  # 図のサイズを大きく

# ============================================================
# デバッグ: VTKファイルの確認
# ============================================================
print(f"[DEBUG] Looking for VTK files in: {vtk_dir}")
print(f"[DEBUG] Directory exists: {os.path.exists(vtk_dir)}")

if os.path.exists(vtk_dir):
    all_files = os.listdir(vtk_dir)
    vtk_files_candidates = [f for f in all_files if f.endswith(".vtk")]
    print(f"[DEBUG] Total .vtk files found: {len(vtk_files_candidates)}")
    if len(vtk_files_candidates) > 0:
        print(f"[DEBUG] First 5 VTK files: {vtk_files_candidates[:5]}")
else:
    print(f"[ERROR] Directory does not exist: {vtk_dir}")
    exit()

# ============================================================
# VTK ファイル整理（block ごと）
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.endswith(".vtk")
])

if not vtk_files:
    print("[ERROR] VTK files not found!")
    exit()

print(f"[INFO] Found {len(vtk_files)} VTK files")

# blockごとにグループ化
block_dict = {}
for f in vtk_files:
    try:
        basename = os.path.basename(f)
        # 様々なパターンに対応
        if '.block' in basename:
            if '.block_' in basename:
                block_id = int(basename.split('.block_')[1].split('.')[0])
            else:
                block_id = int(basename.split('.block')[1].split('.')[0])
        else:
            import re
            numbers = re.findall(r'\d+', basename)
            if numbers:
                block_id = int(numbers[0])
            else:
                block_id = 0
        block_dict.setdefault(block_id, []).append(f)
    except Exception as e:
        print(f"[WARNING] Could not parse block ID from {f}: {e}")
        block_dict.setdefault(0, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]
print(f"[INFO] Number of blocks: {len(all_blocks_sorted)}")

# タイムステップ数を確認
n_timesteps = min(len(block_files) for block_files in all_blocks_sorted)
print(f"[INFO] Number of timesteps: {n_timesteps}")

if n_timesteps == 0:
    print("[ERROR] No timesteps found!")
    exit()

# ============================================================
# 全ステップの密度からvmin/vmaxを決定
# ============================================================
all_densities = []
print("[INFO] Collecting density data for colormap scaling...")

sample_steps = min(10, n_timesteps)

for step_idx in range(sample_steps):
    dens_list = []
    for block_files in all_blocks_sorted:
        if step_idx < len(block_files):
            f = block_files[step_idx]
            try:
                grid = pv.read(f)
                slc = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
                dens_list.append(slc['dens'])
            except Exception as e:
                print(f"[WARNING] Could not read file {f}: {e}")
    
    if dens_list:
        dens_all = np.hstack(dens_list)
        all_densities.append(dens_all)

if all_densities:
    all_densities = np.concatenate(all_densities)
    positive_dens = all_densities[all_densities > 0]
    
    if len(positive_dens) > 0:
        # より細かい構造を見るためにパーセンタイルを調整
        vmin = np.percentile(positive_dens, 2)   # 下限を2%に
        vmax = np.percentile(positive_dens, 98)  # 上限を98%に
        print(f"[INFO] Density scale: vmin={vmin:.3e}, vmax={vmax:.3e}")
    else:
        vmin, vmax = 1e-10, 1e-5
        print("[WARNING] No positive densities found, using default scale")
else:
    print("[ERROR] No density data could be read!")
    exit()

# ============================================================
# 各タイムステップ描画（高解像度版）
# ============================================================
print("[INFO] Generating high-resolution density maps...")

for step_idx in range(n_timesteps):
    
    # 全ブロックの点と密度を収集
    pts_all = []
    dens_all = []
    
    for block_files in all_blocks_sorted:
        if step_idx < len(block_files):
            f = block_files[step_idx]
            try:
                grid = pv.read(f)
                slc = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
                pts_all.append(slc.cell_centers().points)
                dens_all.append(slc['dens'])
            except Exception as e:
                continue
    
    if not pts_all:
        print(f"[WARNING] No data for timestep {step_idx}, skipping...")
        continue
    
    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    
    # グリッド座標の取得（ソートして一意な値を取得）
    xs = np.sort(np.unique(pts_all[:, 0]))
    ys = np.sort(np.unique(pts_all[:, 1]))
    
    # 2D配列にマッピング（より滑らかな表示のために補間を有効化）
    dens_2d = np.full((len(ys), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(np.isclose(xs, p[0], rtol=1e-6, atol=1e-8))[0][0]
        iy = np.where(np.isclose(ys, p[1], rtol=1e-6, atol=1e-8))[0][0]
        dens_2d[iy, ix] = d
    
    # より滑らかな表示のために線形補間を適用
    from scipy.interpolate import griddata
    # 元のグリッド点
    points_2d = np.column_stack([pts_all[:, 0], pts_all[:, 1]])
    # より細かいグリッドを作成（解像度を2倍に）
    xs_fine = np.linspace(x_min, x_max, len(xs) * 2)
    ys_fine = np.linspace(y_min, y_max, len(ys) * 2)
    X_fine, Y_fine = np.meshgrid(xs_fine, ys_fine)
    # 線形補間
    dens_2d_fine = griddata(points_2d, dens_all, (X_fine, Y_fine), method='linear')
    
    # プロット
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_aspect('equal')
    
    # 密度のヒートマップ（高解像度、対数スケール）
    im = ax.pcolormesh(X_fine, Y_fine, dens_2d_fine,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading='auto',
                       rasterized=True)  # ファイルサイズ最適化
    
    # 中心に原点マーク
    ax.plot(xc, yc, 'r+', markersize=12, markeredgewidth=2, 
            label='Center (Origin)', zorder=10)
    
    # 半径を表示するための円を追加（オプション）
    radii = [50, 100, 150, 200]
    for r in radii:
        if r <= Lbox/2:
            theta = np.linspace(0, 2*np.pi, 500)
            ax.plot(xc + r*np.cos(theta), yc + r*np.sin(theta), 
                   'w--', alpha=0.3, linewidth=0.5, zorder=5)
            ax.text(xc + r*0.7, yc + r*0.7, f'r={r}', 
                   color='white', fontsize=8, alpha=0.5, zorder=5)
    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("X", fontsize=14)
    ax.set_ylabel("Y", fontsize=14)
    
    ax.set_title(f"Density Distribution (x-y plane at z=0)\n"
                 f"ρ ∝ r^-1.75, Domain: [{x_min:.0f}, {x_max:.0f}]", 
                 fontsize=14, fontweight='bold')
    
    # カラーバー（より詳細な表示）
    cbar = fig.colorbar(im, ax=ax, label="Density", extend='both')
    cbar.set_label("Density", fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    
    # グリッドを追加（構造を明確に）
    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    
    ax.legend(loc='upper right', framealpha=0.7)
    
    plt.tight_layout()
    
    # 保存（高解像度）
    png = os.path.join(output_dir, f"density_xy_timestep_{step_idx:04d}.png")
    fig.savefig(png, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    if (step_idx + 1) % 5 == 0:
        print(f"[INFO] Processed {step_idx + 1}/{n_timesteps} time steps...")

print(f"[INFO] All high-resolution density maps saved to {output_dir}")
print(f"[INFO] Image size: {figsize[0]*dpi:.0f} x {figsize[1]*dpi:.0f} pixels")

# %% cell 9
#Toyouchi-test4を可視化するのに使った
#半径100以内を拡大表示して全体図の右側に並べて表示する(AMR_OFFバージョン)
#x-y平面用

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from collections import defaultdict

# ============================================================
# 設定
# ============================================================
vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_density_maps"
os.makedirs(output_dir, exist_ok=True)

# 計算領域の設定
Lbox = 448.0
x_min, x_max = -Lbox/2, Lbox/2
y_min, y_max = -Lbox/2, Lbox/2
z_slice = 0.0

xc, yc = 0.0, 0.0
high_res_radius = 100.0
dpi = 200
figsize = (16, 8)

# ============================================================
# VTK ファイル整理（改良版）
# ============================================================
print("[INFO] Organizing VTK files...")

timestep_dict = defaultdict(list)
all_timesteps = set()

for f in os.listdir(vtk_dir):
    if f.endswith(".vtk") and f.startswith("Jeans.block"):
        # ファイル名: Jeans.blockXXX.out2.YYYYY.vtk
        try:
            # タイムステップ番号を抽出（YYYYYの部分）
            # .out2. の後ろから .vtk の前まで
            timestep_str = f.split('.out2.')[1].split('.')[0]
            timestep = int(timestep_str)
            timestep_dict[timestep].append(os.path.join(vtk_dir, f))
            all_timesteps.add(timestep)
        except (IndexError, ValueError) as e:
            print(f"[WARNING] Could not parse timestep from: {f}")
            continue

# タイムステップをソート
timesteps = sorted(all_timesteps)
print(f"[INFO] Found {len(timesteps)} timesteps")
print(f"[INFO] Timestep range: {timesteps[0]} - {timesteps[-1]}")
print(f"[INFO] Files per timestep: {len(timestep_dict[timesteps[0]])}")

# 全タイムステップのリストを確認
print(f"[INFO] First 10 timesteps: {timesteps[:10]}")
print(f"[INFO] Last 10 timesteps: {timesteps[-10:]}")

# ============================================================
# 密度スケールの決定
# ============================================================
print("[INFO] Determining density scale...")

# 全タイムステップから代表的なものをサンプリング
sample_indices = list(range(0, len(timesteps), len(timesteps)//10))[:10]
sample_timesteps = [timesteps[i] for i in sample_indices]
print(f"[INFO] Sampling timesteps: {sample_timesteps}")

all_densities = []

for ts in sample_timesteps:
    dens_list = []
    for f in timestep_dict[ts]:
        try:
            grid = pv.read(f)
            slc = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
            dens_list.append(slc['dens'])
        except Exception as e:
            continue
    if dens_list:
        all_densities.append(np.hstack(dens_list))
        print(f"[INFO] Sampled timestep {ts}: {len(dens_list[-1])} points")

if all_densities:
    all_densities = np.concatenate(all_densities)
    positive_dens = all_densities[all_densities > 0]
    
    if len(positive_dens) > 0:
        vmin = np.percentile(positive_dens, 1)
        vmax = np.percentile(positive_dens, 99)
        print(f"[INFO] Density scale: vmin={vmin:.3e}, vmax={vmax:.3e}")
        print(f"[INFO] Density range: {positive_dens.min():.3e} - {positive_dens.max():.3e}")
    else:
        vmin, vmax = 1e-6, 1e-2
        print(f"[INFO] Using default density scale: {vmin:.3e} - {vmax:.3e}")
else:
    vmin, vmax = 1e-6, 1e-2
    print(f"[INFO] Using default density scale: {vmin:.3e} - {vmax:.3e}")

# ============================================================
# 各タイムステップ描画（全タイムステップ）
# ============================================================
print("[INFO] Generating density maps for all timesteps...")

# 全タイムステップをプロット
timesteps_to_plot = timesteps  # 全タイムステップ

for ts_idx, ts in enumerate(timesteps_to_plot):
    print(f"[INFO] Processing timestep {ts} ({ts_idx+1}/{len(timesteps_to_plot)})")
    
    # 全ブロックの点と密度を収集
    pts_all = []
    dens_all = []
    
    for f in timestep_dict[ts]:
        try:
            grid = pv.read(f)
            slc = grid.slice(normal='z', origin=(0.0, 0.0, z_slice))
            # cell_centersを使うとセル中心の座標が得られる
            pts_all.append(slc.cell_centers().points)
            dens_all.append(slc['dens'])
        except Exception as e:
            continue
    
    if not pts_all:
        print(f"[WARNING] No data for timestep {ts}, skipping...")
        continue
    
    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    
    # デバッグ: 最初の数ステップでデータ確認
    if ts_idx < 3:
        print(f"[DEBUG] Timestep {ts}: {len(pts_all)} points")
        print(f"[DEBUG] X range: {pts_all[:,0].min():.1f} - {pts_all[:,0].max():.1f}")
        print(f"[DEBUG] Y range: {pts_all[:,1].min():.1f} - {pts_all[:,1].max():.1f}")
        print(f"[DEBUG] Density range: {dens_all.min():.3e} - {dens_all.max():.3e}")
    
    # 内側領域のデータ抽出
    r_all = np.sqrt(pts_all[:, 0]**2 + pts_all[:, 1]**2)
    mask_inner = r_all <= high_res_radius
    
    # グリッド生成
    xs = np.sort(np.unique(pts_all[:, 0]))
    ys = np.sort(np.unique(pts_all[:, 1]))
    
    # 2D配列に変換
    dens_2d = np.full((len(ys), len(xs)), np.nan)
    for p, d in zip(pts_all, dens_all):
        ix = np.where(np.isclose(xs, p[0], rtol=1e-5, atol=1e-8))[0]
        iy = np.where(np.isclose(ys, p[1], rtol=1e-5, atol=1e-8))[0]
        if len(ix) > 0 and len(iy) > 0:
            dens_2d[iy[0], ix[0]] = d
    
    X, Y = np.meshgrid(xs, ys)
    
    # 内側領域用の高解像度グリッド
    if mask_inner.sum() > 0:
        pts_inner = pts_all[mask_inner]
        dens_inner = dens_all[mask_inner]
        xs_inner = np.sort(np.unique(pts_inner[:, 0]))
        ys_inner = np.sort(np.unique(pts_inner[:, 1]))
        
        if len(xs_inner) > 1 and len(ys_inner) > 1:
            xs_fine = np.linspace(xs_inner.min(), xs_inner.max(), len(xs_inner) * 2)
            ys_fine = np.linspace(ys_inner.min(), ys_inner.max(), len(ys_inner) * 2)
            X_inner, Y_inner = np.meshgrid(xs_fine, ys_fine)
            
            points_2d = np.column_stack([pts_inner[:, 0], pts_inner[:, 1]])
            dens_inner_fine = griddata(points_2d, dens_inner, (X_inner, Y_inner), method='linear')
        else:
            X_inner, Y_inner, dens_inner_fine = None, None, None
    else:
        X_inner, Y_inner, dens_inner_fine = None, None, None
    
    # ====================================================
    # プロット
    # ====================================================
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[4, 4, 0.3], wspace=0.3)
    
    # 全体図
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_aspect('equal')
    
    im1 = ax1.pcolormesh(X, Y, dens_2d,
                         cmap='inferno',
                         norm=LogNorm(vmin=vmin, vmax=vmax),
                         shading='auto')
    
    # 高解像度領域を示す四角形
    rect = patches.Rectangle(
        (-high_res_radius, -high_res_radius),
        2*high_res_radius, 2*high_res_radius,
        linewidth=2, edgecolor='cyan', facecolor='none',
        linestyle='--', label=f'Zoom region (r={high_res_radius})'
    )
    ax1.add_patch(rect)
    
    # 中心マーク
    ax1.plot(xc, yc, 'r+', markersize=12, markeredgewidth=2, label='Center')
    
    # 半径の円
    for r in [50, 100, 150, 200]:
        if r <= Lbox/2:
            theta = np.linspace(0, 2*np.pi, 200)
            ax1.plot(xc + r*np.cos(theta), yc + r*np.sin(theta), 
                    'w--', alpha=0.3, linewidth=0.8)
            if r in [100, 200]:
                ax1.text(xc + r*0.7, yc + r*0.7, f'r={r}', 
                        color='white', fontsize=8, alpha=0.6)
    
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    ax1.set_xlabel("X", fontsize=12)
    ax1.set_ylabel("Y", fontsize=12)
    ax1.set_title(f"Full Domain (448³)", fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax1.legend(loc='upper right', framealpha=0.7, fontsize=8)
    
    # 拡大図
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_aspect('equal')
    
    if X_inner is not None and dens_inner_fine is not None:
        im2 = ax2.pcolormesh(X_inner, Y_inner, dens_inner_fine,
                             cmap='inferno',
                             norm=LogNorm(vmin=vmin, vmax=vmax),
                             shading='auto')
        
        ax2.plot(xc, yc, 'r+', markersize=12, markeredgewidth=2, label='Center')
        
        for r in [20, 40, 60, 80, 100]:
            theta = np.linspace(0, 2*np.pi, 200)
            ax2.plot(xc + r*np.cos(theta), yc + r*np.sin(theta), 
                    'w--', alpha=0.5, linewidth=0.8)
            ax2.text(xc + r*0.7, yc + r*0.7, f'r={r}', 
                    color='white', fontsize=9, alpha=0.7,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5))
        
        ax2.set_xlim(-high_res_radius, high_res_radius)
        ax2.set_ylim(-high_res_radius, high_res_radius)
        ax2.set_xlabel("X", fontsize=12)
        ax2.set_ylabel("Y", fontsize=12)
        ax2.set_title(f"Zoomed Region (r ≤ {high_res_radius})\nHigh Resolution", 
                     fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        ax2.legend(loc='upper right', framealpha=0.7, fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'No data in inner region', 
                transform=ax2.transAxes, ha='center', va='center', fontsize=12)
        ax2.set_xlim(-high_res_radius, high_res_radius)
        ax2.set_ylim(-high_res_radius, high_res_radius)
    
    # カラーバー
    cax = fig.add_subplot(gs[0, 2])
    cbar = fig.colorbar(im1, cax=cax, label="Density", extend='both')
    cbar.set_label("Density", fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    
    # タイトル
    fig.suptitle(f"Density Distribution (x-y plane at z=0)\n"
                 f"ρ ∝ r^-1.75 | Timestep: {ts:05d} AMR_OFF", 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.05, right=0.95, wspace=0.25)
    
    # 保存
    png = os.path.join(output_dir, f"density_xy_timestep_{ts:05d}_zoom.png")
    fig.savefig(png, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    if (ts_idx + 1) % 10 == 0:
        print(f"[INFO] Processed {ts_idx + 1}/{len(timesteps_to_plot)} timesteps")

print(f"[INFO] All density maps saved to {output_dir}")
print(f"[INFO] Total images: {len(timesteps_to_plot)}")

# %% cell 10
#Toyouchi-test5を可視化するのに使った
#半径100以内を拡大表示して全体図の右側に並べて表示する(AMR_ONバージョン)
#x-y平面用

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from collections import defaultdict

# ============================================================
# 設定
# ============================================================
vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")  # AMRあり
output_dir = "./xy_density_maps"
os.makedirs(output_dir, exist_ok=True)

# 計算領域の設定
Lbox = 448.0
x_min, x_max = -Lbox/2, Lbox/2
y_min, y_max = -Lbox/2, Lbox/2
z_slice = 0.0

xc, yc = 0.0, 0.0
high_res_radius = 100.0
dpi = 200
figsize = (16, 8)

# ============================================================
# VTK ファイル整理
# ============================================================
print("[INFO] Organizing VTK files...")

timestep_dict = defaultdict(list)
all_timesteps = set()

for f in os.listdir(vtk_dir):
    if f.endswith(".vtk") and f.startswith("Jeans.block"):
        try:
            timestep_str = f.split('.out2.')[1].split('.')[0]
            timestep = int(timestep_str)
            timestep_dict[timestep].append(os.path.join(vtk_dir, f))
            all_timesteps.add(timestep)
        except (IndexError, ValueError):
            continue

timesteps = sorted(all_timesteps)
print(f"[INFO] Found {len(timesteps)} timesteps")
print(f"[INFO] Timestep range: {timesteps[0]} - {timesteps[-1]}")
print(f"[INFO] Files per timestep: {len(timestep_dict[timesteps[0]])}")

# ============================================================
# 密度スケールの決定
# ============================================================
print("[INFO] Determining density scale...")

# 全タイムステップから均等にサンプリング（最大10ステップ）
sample_indices = np.linspace(0, len(timesteps)-1, min(10, len(timesteps))).astype(int)
sample_timesteps = [timesteps[i] for i in sample_indices]
print(f"[INFO] Sampling timesteps: {sample_timesteps}")

all_densities = []

for ts in sample_timesteps:
    dens_list = []
    for f in timestep_dict[ts]:
        try:
            grid = pv.read(f)
            dens_list.append(grid['dens'])
        except Exception as e:
            continue
    if dens_list:
        all_densities.append(np.hstack(dens_list))

if all_densities:
    all_densities = np.concatenate(all_densities)
    positive_dens = all_densities[all_densities > 0]
    
    if len(positive_dens) > 0:
        vmin = np.percentile(positive_dens, 1)
        vmax = np.percentile(positive_dens, 99)
        print(f"[INFO] Density scale: vmin={vmin:.3e}, vmax={vmax:.3e}")
        print(f"[INFO] Density range: {positive_dens.min():.3e} - {positive_dens.max():.3e}")
    else:
        vmin, vmax = 1e-6, 1e-2
        print(f"[INFO] Using default density scale: {vmin:.3e} - {vmax:.3e}")
else:
    vmin, vmax = 1e-6, 1e-2
    print(f"[INFO] Using default density scale: {vmin:.3e} - {vmax:.3e}")

# ============================================================
# 各タイムステップ描画（全タイムステップ）
# ============================================================
print("[INFO] Generating density maps for AMR data...")

# 高解像度でプロットするためのグリッド解像度
grid_resolution = 800  # 800x800グリッド（メモリ節約のため少し下げる）

# 全タイムステップを処理
for ts_idx, ts in enumerate(timesteps):
    print(f"[INFO] Processing timestep {ts} ({ts_idx+1}/{len(timesteps)})")
    
    # 全ブロックの点と密度を収集
    points_list = []
    dens_list = []
    
    for f in timestep_dict[ts]:
        try:
            grid = pv.read(f)
            # AMRデータはセル中心の座標を取得
            points = grid.cell_centers().points
            dens = grid['dens']
            
            # z=0付近のデータのみ抽出（z_sliceの近傍）
            z_tolerance = 5.0  # 許容範囲
            mask_z = np.abs(points[:, 2]) <= z_tolerance
            points_list.append(points[mask_z])
            dens_list.append(dens[mask_z])
        except Exception as e:
            continue
    
    if not points_list:
        print(f"[WARNING] No data for timestep {ts}, skipping...")
        continue
    
    pts_all = np.vstack(points_list)
    dens_all = np.hstack(dens_list)
    
    # 最初の数ステップと一定間隔でデバッグ情報を表示
    if ts_idx < 3 or ts_idx % 20 == 0:
        print(f"[DEBUG] Timestep {ts}: {len(pts_all)} points")
        print(f"[DEBUG] X range: {pts_all[:,0].min():.1f} - {pts_all[:,0].max():.1f}")
        print(f"[DEBUG] Y range: {pts_all[:,1].min():.1f} - {pts_all[:,1].max():.1f}")
        print(f"[DEBUG] Density range: {dens_all.min():.3e} - {dens_all.max():.3e}")
    
    # ====================================================
    # 全体図用：高解像度グリッドに補間
    # ====================================================
    # 均一グリッドを作成
    xs_full = np.linspace(x_min, x_max, grid_resolution)
    ys_full = np.linspace(y_min, y_max, grid_resolution)
    X_full, Y_full = np.meshgrid(xs_full, ys_full)
    
    # グリッドポイントに補間（cubicの方が滑らかだが計算重いのでlinear）
    points_2d = np.column_stack([pts_all[:, 0], pts_all[:, 1]])
    dens_full = griddata(points_2d, dens_all, (X_full, Y_full), method='linear')
    
    # ====================================================
    # 内側領域用：より高解像度
    # ====================================================
    mask_inner = np.sqrt(pts_all[:, 0]**2 + pts_all[:, 1]**2) <= high_res_radius
    pts_inner = pts_all[mask_inner]
    dens_inner = dens_all[mask_inner]
    
    if len(pts_inner) > 100:  # 十分なデータがある場合
        inner_resolution = 600  # 内側は高解像度
        xs_inner = np.linspace(-high_res_radius, high_res_radius, inner_resolution)
        ys_inner = np.linspace(-high_res_radius, high_res_radius, inner_resolution)
        X_inner, Y_inner = np.meshgrid(xs_inner, ys_inner)
        
        points_inner_2d = np.column_stack([pts_inner[:, 0], pts_inner[:, 1]])
        dens_inner_fine = griddata(points_inner_2d, dens_inner, (X_inner, Y_inner), method='linear')
    else:
        X_inner, Y_inner, dens_inner_fine = None, None, None
    
    # ====================================================
    # プロット
    # ====================================================
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[4, 4, 0.3], wspace=0.3)
    
    # 全体図
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_aspect('equal')
    
    im1 = ax1.pcolormesh(X_full, Y_full, dens_full,
                         cmap='inferno',
                         norm=LogNorm(vmin=vmin, vmax=vmax),
                         shading='auto',
                         rasterized=True)  # メモリ節約
    
    # 高解像度領域を示す四角形
    rect = patches.Rectangle(
        (-high_res_radius, -high_res_radius),
        2*high_res_radius, 2*high_res_radius,
        linewidth=2, edgecolor='cyan', facecolor='none',
        linestyle='--', label=f'Zoom region (r={high_res_radius})'
    )
    ax1.add_patch(rect)
    
    ax1.plot(xc, yc, 'r+', markersize=12, markeredgewidth=2, label='Center')
    
    for r in [50, 100, 150, 200]:
        if r <= Lbox/2:
            theta = np.linspace(0, 2*np.pi, 200)
            ax1.plot(xc + r*np.cos(theta), yc + r*np.sin(theta), 
                    'w--', alpha=0.3, linewidth=0.8)
            if r in [100, 200]:
                ax1.text(xc + r*0.7, yc + r*0.7, f'r={r}', 
                        color='white', fontsize=8, alpha=0.6)
    
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    ax1.set_xlabel("X", fontsize=12)
    ax1.set_ylabel("Y", fontsize=12)
    ax1.set_title(f"Full Domain (448³) - AMR", fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax1.legend(loc='upper right', framealpha=0.7, fontsize=8)
    
    # 拡大図
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_aspect('equal')
    
    if X_inner is not None and dens_inner_fine is not None:
        im2 = ax2.pcolormesh(X_inner, Y_inner, dens_inner_fine,
                             cmap='inferno',
                             norm=LogNorm(vmin=vmin, vmax=vmax),
                             shading='auto',
                             rasterized=True)
        
        ax2.plot(xc, yc, 'r+', markersize=12, markeredgewidth=2, label='Center')
        
        for r in [20, 40, 60, 80, 100]:
            theta = np.linspace(0, 2*np.pi, 200)
            ax2.plot(xc + r*np.cos(theta), yc + r*np.sin(theta), 
                    'w--', alpha=0.5, linewidth=0.8)
            ax2.text(xc + r*0.7, yc + r*0.7, f'r={r}', 
                    color='white', fontsize=9, alpha=0.7,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5))
        
        ax2.set_xlim(-high_res_radius, high_res_radius)
        ax2.set_ylim(-high_res_radius, high_res_radius)
        ax2.set_xlabel("X", fontsize=12)
        ax2.set_ylabel("Y", fontsize=12)
        ax2.set_title(f"Zoomed Region (r ≤ {high_res_radius})\nAMR Data", 
                     fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        ax2.legend(loc='upper right', framealpha=0.7, fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'No data in inner region', 
                transform=ax2.transAxes, ha='center', va='center', fontsize=12)
        ax2.set_xlim(-high_res_radius, high_res_radius)
        ax2.set_ylim(-high_res_radius, high_res_radius)
    
    # カラーバー
    cax = fig.add_subplot(gs[0, 2])
    cbar = fig.colorbar(im1, cax=cax, label="Density", extend='both')
    cbar.set_label("Density", fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    
    # タイトル
    fig.suptitle(f"Density Distribution (x-y plane at z=0)\n"
                 f"ρ ∝ r^-1.75 | Timestep: {ts:05d} | AMR_ON", 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.05, right=0.95, wspace=0.25)
    
    # 保存
    png = os.path.join(output_dir, f"density_xy_timestep_{ts:05d}_zoom_amr.png")
    fig.savefig(png, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    # 進捗表示
    if (ts_idx + 1) % 10 == 0:
        print(f"[INFO] Processed {ts_idx + 1}/{len(timesteps)} timesteps")

print(f"[INFO] All AMR density maps saved to {output_dir}")
print(f"[INFO] Total images: {len(timesteps)}")

# %% cell 11
#ｘ−ｙ平面のベクトルon密度mapの画像作成用コード（jeans_3d-test14の可視化に使ったやつ）

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from PIL import Image
import matplotlib.animation as animation

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_vector_field"
os.makedirs(output_dir, exist_ok=True)

# 計算領域
Lbox = 12.0
x_min, x_max = 0.0, 12.0
y_min, y_max = 0.0, 12.0

# 球
Rsphere = 5.0
xc, yc = 6.0, 6.0

# 中央スライス
z_slice = 6.0

# ============================================================
# VTK ファイル取得
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    basename = os.path.basename(f)
    block_id = int(basename.split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内のみで密度スケール決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):
    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))
        pts_all.append(sl.cell_centers().points)
        dens_all.append(sl['dens'])

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    r = np.sqrt((pts_all[:, 0] - xc)**2 + (pts_all[:, 1] - yc)**2)
    dens_inside_all.append(dens_all[r <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"Density scale (inside sphere): vmin={vmin}, vmax={vmax}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all, vx_all, vy_all = [], [], [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))

        centers = sl.cell_centers().points
        dens = sl['dens']
        mom = sl['mom']

        pts_all.append(centers)
        dens_all.append(dens)
        vx_all.append(mom[:, 0] / dens)
        vy_all.append(mom[:, 1] / dens)

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    vx_all = np.hstack(vx_all)
    vy_all = np.hstack(vy_all)

    xs = np.unique(pts_all[:, 0])
    ys = np.unique(pts_all[:, 1])
    nx, ny = len(xs), len(ys)

    dens_2d = np.full((ny, nx), np.nan)
    vx_2d = np.full((ny, nx), np.nan)
    vy_2d = np.full((ny, nx), np.nan)

    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            m = (np.isclose(pts_all[:, 0], x)) & (np.isclose(pts_all[:, 1], y))
            if np.any(m):
                dens_2d[j, i] = dens_all[m][0]
                vx_2d[j, i] = vx_all[m][0]
                vy_2d[j, i] = vy_all[m][0]

    # 球内・球外マスク
    X, Y = np.meshgrid(xs, ys)
    r = np.sqrt((X - xc)**2 + (Y - yc)**2)

    dens_in  = np.ma.masked_where(r > Rsphere, dens_2d)
    dens_out = np.ma.masked_where(r <= Rsphere, dens_2d)

    # ========================================================
    # 描画
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')

    # 球外：やや濃い青背景
    ax.pcolormesh(xs, ys, dens_out,
                  color='#0b3c5d')

    # 球内：密度（対数）
    im = ax.pcolormesh(xs, ys, dens_in,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # ベクトル
    skip = max(1, nx // 80)
    ax.quiver(xs[::skip], ys[::skip],
              vx_2d[::skip, ::skip],
              vy_2d[::skip, ::skip],
              color='white',
              scale=12,
              linewidth=0.15,
              width=0.002)

    # 球境界（細め・実線）
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            yc + Rsphere*np.sin(theta),
            color='cyan',
            lw=0.6,
            linestyle='-')

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    filename = os.path.basename(step_files[0])
    num = filename.split('.')[-2]
    ax.set_title(f"file={filename}  (Rsphere=5)") #タイトルは適切に変えると良い

    fig.colorbar(im, ax=ax, label="Density (inside sphere)")

    png_file = os.path.join(output_dir, f"vector_{num}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

# %% cell 12
#ｘ−z平面のベクトルon密度mapの画像作成用コード（jeans_3d-test27の可視化に使ったやつ）

import os
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import natsort
from matplotlib.colors import LogNorm
from PIL import Image
import matplotlib.animation as animation

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xz_vector_field"
os.makedirs(output_dir, exist_ok=True)

# 計算領域
Lbox = 12.0
x_min, x_max = 0.0, 12.0
z_min, z_max = 0.0, 12.0

# 球
Rsphere = 5.0
xc, zc = 6.0, 6.0

# 中央スライス（y）
y_slice = 6.0

# ============================================================
# VTK ファイル取得
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    block_id = int(os.path.basename(f).split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内のみで密度スケール決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='y', origin=(0, y_slice, 0))
        centers = sl.cell_centers().points

        pts_all.append(centers)
        dens_all.append(sl['dens'])

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    r = np.sqrt((pts_all[:, 0] - xc)**2 +
                (pts_all[:, 2] - zc)**2)

    dens_inside_all.append(dens_all[r <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"Density scale (inside sphere): vmin={vmin:.3e}, vmax={vmax:.3e}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all, vx_all, vz_all = [], [], [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='y', origin=(0, y_slice, 0))

        centers = sl.cell_centers().points
        dens = sl['dens']
        mom  = sl['mom']

        pts_all.append(centers)
        dens_all.append(dens)
        vx_all.append(mom[:, 0] / dens)
        vz_all.append(mom[:, 2] / dens)

    pts_all  = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    vx_all   = np.hstack(vx_all)
    vz_all   = np.hstack(vz_all)

    xs = np.unique(pts_all[:, 0])
    zs = np.unique(pts_all[:, 2])
    nx, nz = len(xs), len(zs)

    dens_2d = np.full((nz, nx), np.nan)
    vx_2d   = np.full((nz, nx), np.nan)
    vz_2d   = np.full((nz, nx), np.nan)

    for i, x in enumerate(xs):
        for k, z in enumerate(zs):
            m = (np.isclose(pts_all[:, 0], x)) & \
                (np.isclose(pts_all[:, 2], z))
            if np.any(m):
                dens_2d[k, i] = dens_all[m][0]
                vx_2d[k, i]   = vx_all[m][0]
                vz_2d[k, i]   = vz_all[m][0]

    # 球内・球外マスク
    X, Z = np.meshgrid(xs, zs)
    r = np.sqrt((X - xc)**2 + (Z - zc)**2)

    dens_in  = np.ma.masked_where(r > Rsphere, dens_2d)
    dens_out = np.ma.masked_where(r <= Rsphere, dens_2d)

    # ========================================================
    # 描画
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')

    # 球外：濃い青背景
    ax.pcolormesh(xs, zs, dens_out,
                  color='#0b3c5d')

    # 球内：密度（対数）
    im = ax.pcolormesh(xs, zs, dens_in,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # ベクトル
    skip = max(1, nx // 80)
    ax.quiver(xs[::skip], zs[::skip],
              vx_2d[::skip, ::skip],
              vz_2d[::skip, ::skip],
              color='white',
              scale=12,
              linewidth=0.15,
              width=0.002)

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            zc + Rsphere*np.sin(theta),
            color='cyan',
            lw=0.6)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_xlabel("X")
    ax.set_ylabel("Z")

    filename = os.path.basename(step_files[0])
    num = filename.split('.')[-2]
    ax.set_title(f"file={filename}  (γ=1.2,β=0.07)") #タイトルは適切に変えると良い

    fig.colorbar(im, ax=ax, label="Density (inside sphere)")

    png_file = os.path.join(output_dir, f"vector_{num}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

# %% cell 13
# ｘ−ｙ平面のベクトル on 密度 map 可視化コード（ベクトル密度 1/3 版）

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from PIL import Image
import matplotlib.animation as animation

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_vector_field"
os.makedirs(output_dir, exist_ok=True)

# 計算領域
Lbox = 12.0
x_min, x_max = 0.0, 12.0
y_min, y_max = 0.0, 12.0

# 球
Rsphere = 5.0
xc, yc = 6.0, 6.0

# 中央スライス
z_slice = 6.0

# ============================================================
# VTK ファイル取得
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    basename = os.path.basename(f)
    block_id = int(basename.split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内のみで密度スケール決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):
    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))
        pts_all.append(sl.cell_centers().points)
        dens_all.append(sl['dens'])

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    r = np.sqrt((pts_all[:, 0] - xc)**2 + (pts_all[:, 1] - yc)**2)
    dens_inside_all.append(dens_all[r <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"Density scale (inside sphere): vmin={vmin}, vmax={vmax}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all, vx_all, vy_all = [], [], [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))

        centers = sl.cell_centers().points
        dens = sl['dens']
        mom = sl['mom']

        pts_all.append(centers)
        dens_all.append(dens)
        vx_all.append(mom[:, 0] / dens)
        vy_all.append(mom[:, 1] / dens)

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    vx_all = np.hstack(vx_all)
    vy_all = np.hstack(vy_all)

    xs = np.unique(pts_all[:, 0])
    ys = np.unique(pts_all[:, 1])
    nx, ny = len(xs), len(ys)

    dens_2d = np.full((ny, nx), np.nan)
    vx_2d = np.full((ny, nx), np.nan)
    vy_2d = np.full((ny, nx), np.nan)

    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            m = (np.isclose(pts_all[:, 0], x)) & (np.isclose(pts_all[:, 1], y))
            if np.any(m):
                dens_2d[j, i] = dens_all[m][0]
                vx_2d[j, i] = vx_all[m][0]
                vy_2d[j, i] = vy_all[m][0]

    # 球内・球外マスク
    X, Y = np.meshgrid(xs, ys)
    r = np.sqrt((X - xc)**2 + (Y - yc)**2)

    dens_in  = np.ma.masked_where(r > Rsphere, dens_2d)
    dens_out = np.ma.masked_where(r <= Rsphere, dens_2d)

    # ========================================================
    # 描画
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')

    # 球外：背景
    ax.pcolormesh(xs, ys, dens_out, color='#0b3c5d')

    # 球内：密度（対数）
    im = ax.pcolormesh(xs, ys, dens_in,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # ベクトル（密度を 1/3 に）
    base_skip = max(1, nx // 80)
    skip = base_skip * 3

    ax.quiver(xs[::skip], ys[::skip],
              vx_2d[::skip, ::skip],
              vy_2d[::skip, ::skip],
              color='white',
              scale=12,
              linewidth=0.15,
              width=0.002)

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            yc + Rsphere*np.sin(theta),
            color='cyan',
            lw=0.6)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    filename = os.path.basename(step_files[0])
    num = filename.split('.')[-2]
    ax.set_title(f"file={filename}  (Rsphere=5)")

    fig.colorbar(im, ax=ax, label="Density (inside sphere)")

    png_file = os.path.join(output_dir, f"vector_{num}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

# %% cell 14
# x–z 平面のベクトル on 密度 map 可視化コード（ベクトル密度 1/3）


import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import natsort
from matplotlib.colors import LogNorm
from PIL import Image
import matplotlib.animation as animation

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xz_vector_field"
os.makedirs(output_dir, exist_ok=True)

# 計算領域
Lbox = 12.0
x_min, x_max = 0.0, 12.0
z_min, z_max = 0.0, 12.0

# 球
Rsphere = 5.0
xc, zc = 6.0, 6.0

# 中央スライス
y_slice = 6.0

# ============================================================
# VTK ファイル取得
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    basename = os.path.basename(f)
    block_id = int(basename.split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内のみで密度スケール決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):
    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='y', origin=(0, y_slice, 0))
        pts_all.append(sl.cell_centers().points)
        dens_all.append(sl['dens'])

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    r = np.sqrt((pts_all[:, 0] - xc)**2 + (pts_all[:, 2] - zc)**2)
    dens_inside_all.append(dens_all[r <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

print(f"Density scale (inside sphere): vmin={vmin}, vmax={vmax}")

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all, vx_all, vz_all = [], [], [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='y', origin=(0, y_slice, 0))

        centers = sl.cell_centers().points
        dens = sl['dens']
        mom = sl['mom']

        pts_all.append(centers)
        dens_all.append(dens)
        vx_all.append(mom[:, 0] / dens)
        vz_all.append(mom[:, 2] / dens)

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    vx_all = np.hstack(vx_all)
    vz_all = np.hstack(vz_all)

    xs = np.unique(pts_all[:, 0])
    zs = np.unique(pts_all[:, 2])
    nx, nz = len(xs), len(zs)

    dens_2d = np.full((nz, nx), np.nan)
    vx_2d   = np.full((nz, nx), np.nan)
    vz_2d   = np.full((nz, nx), np.nan)

    for i, x in enumerate(xs):
        for j, z in enumerate(zs):
            m = (np.isclose(pts_all[:, 0], x)) & (np.isclose(pts_all[:, 2], z))
            if np.any(m):
                dens_2d[j, i] = dens_all[m][0]
                vx_2d[j, i]   = vx_all[m][0]
                vz_2d[j, i]   = vz_all[m][0]

    # 球内・球外マスク
    X, Z = np.meshgrid(xs, zs)
    r = np.sqrt((X - xc)**2 + (Z - zc)**2)

    dens_in  = np.ma.masked_where(r > Rsphere, dens_2d)
    dens_out = np.ma.masked_where(r <= Rsphere, dens_2d)

    # ========================================================
    # 描画
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')

    # 球外：背景
    ax.pcolormesh(xs, zs, dens_out, color='#0b3c5d')

    # 球内：密度（対数）
    im = ax.pcolormesh(xs, zs, dens_in,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # ベクトル（密度 1/3）
    base_skip = max(1, nx // 80)
    skip = base_skip * 3

    ax.quiver(xs[::skip], zs[::skip],
              vx_2d[::skip, ::skip],
              vz_2d[::skip, ::skip],
              color='white',
              scale=12,
              linewidth=0.15,
              width=0.002)

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            zc + Rsphere*np.sin(theta),
            color='cyan',
            lw=0.6)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_xlabel("X")
    ax.set_ylabel("Z")

    filename = os.path.basename(step_files[0])
    num = filename.split('.')[-2]
    ax.set_title(f"file={filename}  (Rsphere=5)")

    fig.colorbar(im, ax=ax, label="Density (inside sphere)")

    png_file = os.path.join(output_dir, f"vector_{num}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

# %% cell 15
#x-y平面に流線を描画するコード（jeans_3d-test26の可視化で使った）

import os
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import natsort
import matplotlib.animation as animation
from PIL import Image

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")
output_dir = "./xy_streamline"
os.makedirs(output_dir, exist_ok=True)

# 計算領域
Lbox = 12.0
x_min, x_max = 0.0, 12.0
y_min, y_max = 0.0, 12.0

# 球
Rsphere = 5.0
xc, yc = 6.0, 6.0

# 中央スライス
z_slice = 6.0

# ============================================================
# VTK ファイル取得
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Jeans.block") and f.endswith(".vtk")
])

block_dict = {}
for f in vtk_files:
    block_id = int(os.path.basename(f).split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]

# ============================================================
# 球内のみで密度スケール決定
# ============================================================
dens_inside_all = []

for step_files in zip(*all_blocks_sorted):
    pts_all, dens_all = [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))
        pts_all.append(sl.cell_centers().points)
        dens_all.append(sl['dens'])

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)

    r = np.sqrt((pts_all[:, 0]-xc)**2 + (pts_all[:, 1]-yc)**2)
    dens_inside_all.append(dens_all[r <= Rsphere])

dens_inside_all = np.concatenate(dens_inside_all)
vmin = np.percentile(dens_inside_all, 5)
vmax = np.percentile(dens_inside_all, 95)

# ============================================================
# 各タイムステップ描画
# ============================================================
for step_files in zip(*all_blocks_sorted):

    pts_all, dens_all, vx_all, vy_all = [], [], [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))

        centers = sl.cell_centers().points
        dens = sl['dens']
        mom  = sl['mom']

        pts_all.append(centers)
        dens_all.append(dens)
        vx_all.append(mom[:, 0] / dens)
        vy_all.append(mom[:, 1] / dens)

    pts_all = np.vstack(pts_all)
    dens_all = np.hstack(dens_all)
    vx_all = np.hstack(vx_all)
    vy_all = np.hstack(vy_all)

    xs = np.unique(pts_all[:, 0])
    ys = np.unique(pts_all[:, 1])
    nx, ny = len(xs), len(ys)

    dens_2d = np.full((ny, nx), np.nan)
    vx_2d   = np.full((ny, nx), np.nan)
    vy_2d   = np.full((ny, nx), np.nan)

    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            m = (np.isclose(pts_all[:, 0], x)) & (np.isclose(pts_all[:, 1], y))
            if np.any(m):
                dens_2d[j, i] = dens_all[m][0]
                vx_2d[j, i]   = vx_all[m][0]
                vy_2d[j, i]   = vy_all[m][0]

    X, Y = np.meshgrid(xs, ys)
    r = np.sqrt((X-xc)**2 + (Y-yc)**2)

    dens_in  = np.ma.masked_where(r > Rsphere, dens_2d)
    dens_out = np.ma.masked_where(r <= Rsphere, dens_2d)

    # ========================================================
    # 描画
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')

    # 球外
    ax.pcolormesh(xs, ys, dens_out, color='#0b3c5d')

    # 球内密度
    im = ax.pcolormesh(xs, ys, dens_in,
                       cmap='inferno',
                       norm=LogNorm(vmin=vmin, vmax=vmax))

    # ===== 流線 =====
    ax.streamplot(
        xs, ys,
        vx_2d, vy_2d,
        color='white',
        density=1.2,
        linewidth=0.7,
        arrowsize=0.8
    )

    # 球境界
    theta = np.linspace(0, 2*np.pi, 400)
    ax.plot(xc + Rsphere*np.cos(theta),
            yc + Rsphere*np.sin(theta),
            color='cyan', lw=0.6)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    filename = os.path.basename(step_files[0])
    num = filename.split('.')[-2]
    ax.set_title(f"Streamlines (γ=1.2, β=0.05)")

    fig.colorbar(im, ax=ax, label="Density (inside sphere)")

    png_file = os.path.join(output_dir, f"stream_{num}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

# %% cell 16
#Toyouchi-test12を可視化するときに使った、流線on密度mapを作成するコード
#Toyouchi+23の再現は今後これを使用すると良い

import os
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import natsort
from scipy.interpolate import griddata

# ============================================================
# ユーザー設定
# ============================================================
vtk_dir = os.path.expanduser("~/athena-project/results/〇〇")  # 適宜変更
output_dir = "./x-y_streamline"
os.makedirs(output_dir, exist_ok=True)

# 計算領域（code unit）
x_min, x_max = -224, 224
y_min, y_max = -224, 224

# スライスするz座標（2D計算なら0.0でOK）
z_slice = 0.0

# 可視化する変数（Athena++の出力に合わせる）
var_name = 'rho'      # 密度

# 流線の密度（数値を大きくすると流線が増える）
stream_density = 1.5
stream_linewidth = 0.35   # 線幅
stream_color = 'white'    # 流線の色

# カラーマップの範囲（自動設定する場合はNone）
vmin_user = None
vmax_user = None

# 内挿グリッドの解像度（数値を大きくすると高解像度）
grid_resolution = 200  # 200x200の均一格子上に内挿

# ============================================================
# VTKファイルの読み込み
# ============================================================
vtk_files = natsort.natsorted([
    os.path.join(vtk_dir, f)
    for f in os.listdir(vtk_dir)
    if f.startswith("Toyouchi.block") and f.endswith(".vtk")
])

print(f"Found {len(vtk_files)} VTK files")

# ブロックごとにグループ化
block_dict = {}
for f in vtk_files:
    block_id = int(os.path.basename(f).split('.')[1].replace('block', ''))
    block_dict.setdefault(block_id, []).append(f)

all_blocks_sorted = [block_dict[b] for b in sorted(block_dict.keys())]
print(f"Number of blocks: {len(all_blocks_sorted)}")
print(f"Number of time steps: {len(all_blocks_sorted[0]) if all_blocks_sorted else 0}")

if not all_blocks_sorted:
    print("No VTK files found. Check the directory path.")
    exit()

# ============================================================
# 密度スケールを自動決定（全タイムステップから）
# ============================================================
if vmin_user is None or vmax_user is None:
    var_all_steps = []
    for step_idx, step_files in enumerate(zip(*all_blocks_sorted)):
        print(f"Scanning time step {step_idx+1}/{len(all_blocks_sorted[0])} for color scale...")
        for f in step_files:
            grid = pv.read(f)
            sl = grid.slice(normal='z', origin=(0, 0, z_slice))
            var_all_steps.append(sl[var_name])
    
    var_all_steps = np.hstack(var_all_steps)
    vmin = vmin_user if vmin_user else np.percentile(var_all_steps, 2)
    vmax = vmax_user if vmax_user else np.percentile(var_all_steps, 98)
else:
    vmin, vmax = vmin_user, vmax_user

print(f"Color scale: vmin={vmin:.3e}, vmax={vmax:.3e}")

# ============================================================
# 各タイムステップで描画
# ============================================================
for step_idx, step_files in enumerate(zip(*all_blocks_sorted)):
    print(f"Processing time step {step_idx+1}/{len(all_blocks_sorted[0])}...")

    pts_all, var_all, vx_all, vy_all = [], [], [], []

    for f in step_files:
        grid = pv.read(f)
        sl = grid.slice(normal='z', origin=(0, 0, z_slice))

        centers = sl.cell_centers().points
        var = sl['rho']
        vel = sl['vel']
        vx = vel[:, 0]
        vy = vel[:, 1]

        pts_all.append(centers)
        var_all.append(var)
        vx_all.append(vx)
        vy_all.append(vy)

    pts_all = np.vstack(pts_all)
    var_all = np.hstack(var_all)
    vx_all = np.hstack(vx_all)
    vy_all = np.hstack(vy_all)

    # 計算領域内の点のみを使用
    mask = (pts_all[:, 0] >= x_min) & (pts_all[:, 0] <= x_max) & \
           (pts_all[:, 1] >= y_min) & (pts_all[:, 1] <= y_max)
    
    pts_all = pts_all[mask]
    var_all = var_all[mask]
    vx_all = vx_all[mask]
    vy_all = vy_all[mask]

    # 等間隔なグリッドを作成
    xi = np.linspace(x_min, x_max, grid_resolution)
    yi = np.linspace(y_min, y_max, grid_resolution)
    X, Y = np.meshgrid(xi, yi)

    # 2次元の点群として内挿（x, y座標のみ使用）
    print("  Interpolating data to regular grid...")
    points_2d = pts_all[:, :2]  # x, y座標のみ抽出
    
    try:
        # method='linear'で内挿
        var_2d = griddata(points_2d, var_all, (X, Y), method='linear', fill_value=np.nan)
        vx_2d = griddata(points_2d, vx_all, (X, Y), method='linear', fill_value=np.nan)
        vy_2d = griddata(points_2d, vy_all, (X, Y), method='linear', fill_value=np.nan)
    except Exception as e:
        print(f"  Linear interpolation failed: {e}")
        print("  Falling back to nearest interpolation...")
        # フォールバック：nearest法を使用
        var_2d = griddata(points_2d, var_all, (X, Y), method='nearest', fill_value=np.nan)
        vx_2d = griddata(points_2d, vx_all, (X, Y), method='nearest', fill_value=np.nan)
        vy_2d = griddata(points_2d, vy_all, (X, Y), method='nearest', fill_value=np.nan)

    # NaNをマスク（領域外のデータを非表示に）
    var_2d = np.ma.masked_invalid(var_2d)
    vx_2d = np.ma.masked_invalid(vx_2d)
    vy_2d = np.ma.masked_invalid(vy_2d)

    # データがすべてNaNの場合のチェック
    if np.all(var_2d.mask):
        print(f"  Warning: No valid data after interpolation for step {step_idx}")
        continue

    # ========================================================
    # 描画
    # ========================================================
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    # 密度分布（対数スケール）
    try:
        im = ax.pcolormesh(X, Y, var_2d,
                           cmap='plasma',
                           norm=LogNorm(vmin=vmin, vmax=vmax),
                           shading='auto')
    except Exception as e:
        print(f"  pcolormesh failed: {e}")
        # フォールバック：imshowを使用
        im = ax.imshow(var_2d.T, origin='lower', extent=[x_min, x_max, y_min, y_max],
                       cmap='plasma', norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    # ===== 流線 =====
    # 速度場の極端な値をクリップ
    vmag = np.sqrt(vx_2d**2 + vy_2d**2)
    valid_vmag = vmag[~np.ma.masked_invalid(vmag).mask]
    
    if len(valid_vmag) > 0:
        vmax_clip = np.percentile(valid_vmag, 99)
        
        vx_clipped = np.clip(vx_2d.filled(0), -vmax_clip, vmax_clip)
        vy_clipped = np.clip(vy_2d.filled(0), -vmax_clip, vmax_clip)
        
        # 流線描画
        try:
            ax.streamplot(
                xi, yi,
                vx_clipped, vy_clipped,
                color=stream_color,
                density=stream_density,
                linewidth=stream_linewidth,
                arrowsize=0.6,
                arrowstyle='->',
                minlength=0.1,
                broken_streamlines=False
            )
        except Exception as e:
            print(f"  Streamplot failed: {e}")
    else:
        print("  Warning: No valid velocity data for streamlines")

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("X (×450AU)", fontsize=12)
    ax.set_ylabel("Y (×450AU)", fontsize=12)

    # タイムステップ情報
    ax.set_title(f"Density & Streamlines (t = step {step_idx})", fontsize=12)

    cbar = fig.colorbar(im, ax=ax, label=f"Density (code unit)")
    cbar.ax.tick_params(labelsize=10)

    png_file = os.path.join(output_dir, f"stream_{step_idx:04d}.png")
    fig.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved {png_file}")

print("All processing completed!")

# %% cell 17
#計算領域の各境界面での流出境界条件をチェックするコード
#出力結果がすべての面で０であればいい

import pyvista as pv
import numpy as np
import glob
import re
import os

# ============================
# VTKファイル取得
# ============================

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")

vtk_files = sorted(glob.glob(os.path.join(vtk_dir, "Jeans.block*.out2.*.vtk")))

print("Total VTK files:", len(vtk_files))

# timestep番号取得
steps = sorted(set(re.search(r"out2\.(\d+)", f).group(1) for f in vtk_files))

last_step = steps[-1]

files = sorted(glob.glob(os.path.join(vtk_dir, f"Jeans.block*.out2.{last_step}.vtk")))

print("\nFiles used for last timestep:")
for f in files[:5]:
    print(f)
print("...")

# ============================
# MeshBlock数取得
# ============================

block_ids = sorted(set(int(re.search(r"block(\d+)", f).group(1)) for f in vtk_files))

nblock = len(block_ids)

print("\nMeshBlocks:", nblock)

# ============================
# 全ブロックの座標範囲取得
# ============================

xmin = []
xmax = []
ymin = []
ymax = []
zmin = []
zmax = []

datasets = []

for f in files:

    grid = pv.read(f)

    bounds = grid.bounds
    x0,x1,y0,y1,z0,z1 = bounds

    xmin.append(x0)
    xmax.append(x1)
    ymin.append(y0)
    ymax.append(y1)
    zmin.append(z0)
    zmax.append(z1)

    datasets.append(grid)

# ============================
# 計算領域のグローバル境界
# ============================

XMIN = min(xmin)
XMAX = max(xmax)

YMIN = min(ymin)
YMAX = max(ymax)

ZMIN = min(zmin)
ZMAX = max(zmax)

print("\nGlobal domain:")
print("x:", XMIN, XMAX)
print("y:", YMIN, YMAX)
print("z:", ZMIN, ZMAX)

# ============================
# inflow判定
# ============================

tol = 1e-10

count = {
    "left":0,
    "right":0,
    "front":0,
    "back":0,
    "bottom":0,
    "top":0
}

for grid in datasets:

    rho = grid["dens"]
    mom = grid["mom"]

    vx = mom[:,0] / rho
    vy = mom[:,1] / rho
    vz = mom[:,2] / rho

    # セル中心座標
    centers = grid.cell_centers().points

    x = centers[:,0]
    y = centers[:,1]
    z = centers[:,2]

    # ============================
    # left boundary
    # ============================

    mask = np.abs(x - XMIN) < tol
    count["left"] += np.sum(vx[mask] > 0)  #流出条件v・n>0を満たさない（つまり外部からの流入）セルを数えている。
    # →これが０であれば流出境界条件を満たしていることになる
    # ============================
    # right boundary
    # ============================

    mask = np.abs(x - XMAX) < tol
    count["right"] += np.sum(vx[mask] < 0)

    # ============================
    # front
    # ============================

    mask = np.abs(y - YMIN) < tol
    count["front"] += np.sum(vy[mask] > 0)

    # ============================
    # back
    # ============================

    mask = np.abs(y - YMAX) < tol
    count["back"] += np.sum(vy[mask] < 0)

    # ============================
    # bottom
    # ============================

    mask = np.abs(z - ZMIN) < tol
    count["bottom"] += np.sum(vz[mask] > 0)

    # ============================
    # top
    # ============================

    mask = np.abs(z - ZMAX) < tol
    count["top"] += np.sum(vz[mask] < 0)

# ============================
# 結果表示
# ============================

print("\n=== GLOBAL inflow cells ===")

print("left  :", count["left"])
print("right :", count["right"])

print("front :", count["front"])
print("back  :", count["back"])

print("bottom:", count["bottom"])
print("top   :", count["top"])

# %% cell 18
# ============================================================
# AMR密度マップ（x-y平面・AMRグリッドあり）
#
# 表示単位
#   距離：AU
#   質量：M_sun
#   密度：M_sun AU^-3
#
# 仕様
#   1. z=0と交差するAMRブロックのみ使用
#   2. 各ブロックからz=0に最も近い1層だけ抽出
#   3. 全体図とズーム図の両方にAMR境界を表示
#   4. 線形補間＋最近傍補間で穴埋め
#   5. ズーム図はbilinear表示
#   6. ズーム半径：2×10^4～5×10^4 AU
#   7. ズーム色範囲：シンク外密度の5～98パーセンタイル
#   8. 最大・最小密度：シンク外かつズーム領域内で検索
#   9. 密度マップ自体にはシンク内部も表示
#  10. VTKヘッダーのコード時刻をyrへ変換
#  11. ファイル名はタイムステップ順
# ============================================================

import os
import re
from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogFormatterSciNotation
from scipy.interpolate import griddata


# ============================================================
# 入出力設定
# ============================================================
vtk_dir = os.path.expanduser(
    "~/athena-project/results/〇〇"
)

output_dir = "./xy_density_maps_with_grid"
os.makedirs(output_dir, exist_ok=True)


# ============================================================
# 単位定義
# ============================================================
# Toyouchi.cppのコード単位
M_UNIT_CGS = 4.0e33   # g
L_UNIT_CGS = 7.03e15  # cm (L0 = rb / 2)
T_UNIT_CGS = 3.61e10  # s

# 物理定数
AU_CGS = 1.495978707e13
MSUN_CGS = 1.98847e33
YEAR_CGS = 365.25 * 24.0 * 3600.0

# コード単位 → 表示単位
LENGTH_UNIT_AU = L_UNIT_CGS / AU_CGS
MASS_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS
TIME_UNIT_YR = T_UNIT_CGS / YEAR_CGS

# code density → M_sun AU^-3
DENSITY_UNIT = (
    MASS_UNIT_MSUN / LENGTH_UNIT_AU**3
)

print("[INFO] Unit conversion factors:")
print(
    f"  1 code length  = "
    f"{LENGTH_UNIT_AU:.6e} AU"
)
print(
    f"  1 code mass    = "
    f"{MASS_UNIT_MSUN:.6e} M_sun"
)
print(
    f"  1 code time    = "
    f"{TIME_UNIT_YR:.6e} yr"
)
print(
    f"  1 code density = "
    f"{DENSITY_UNIT:.6e} M_sun AU^-3"
)


# ============================================================
# 計算領域
# ============================================================
# 新しい入力ファイル：
# x1, x2, x3 = -224 ～ +224 code length
LBOX_CODE = 448.0

LBOX_AU = (
    LBOX_CODE * LENGTH_UNIT_AU
)

DOMAIN_MIN = -LBOX_AU / 2.0
DOMAIN_MAX = LBOX_AU / 2.0


# ============================================================
# AMR設定
# ============================================================
BASE_NX = 64
MAX_AMR_LEVEL = 5

BASE_DX_CODE = (
    LBOX_CODE / BASE_NX
)

BASE_DX_AU = (
    BASE_DX_CODE * LENGTH_UNIT_AU
)

LEVEL_COLORS = {
    0: "gray",
    1: "blue",
    2: "green",
    3: "black",
    4: "red",
    5: "purple",
}


# ============================================================
# 描画・ズーム設定
# ============================================================
FULL_RESOLUTION = 800
ZOOM_RESOLUTION = 800

# 自動ズーム倍率
ZOOM_RADIUS_FACTOR = 20.0

# 新しい計算領域用
# AU単位で直接指定
ZOOM_RADIUS_MIN_AU = 2.0e4
ZOOM_RADIUS_MAX_AU = 5.0e4

# 全体図の色範囲
FULL_PERCENTILES = (1.0, 99.0)

# ズーム図の色範囲
ZOOM_PERCENTILES = (5.0, 98.0)

FULL_CMAP = "inferno"
ZOOM_CMAP = "turbo"
ZOOM_INTERPOLATION = "bilinear"

dpi = 200
figsize = (19, 8)


# ============================================================
# シンクマスク設定
# ============================================================
# inputファイルのr_sink_auと一致させる
SINK_RADIUS_AU = 1000.0

# シンク境界付近も除外したい場合は
# 1.1～1.2程度に変更する
SINK_MASK_FACTOR = 1.0

SINK_MASK_RADIUS_AU = (
    SINK_MASK_FACTOR * SINK_RADIUS_AU
)

print("[INFO] Domain and mask settings:")
print(
    f"  Domain width = {LBOX_AU:.6e} AU"
)
print(
    f"  Level 0 dx   = {BASE_DX_AU:.6e} AU"
)
print(
    f"  Sink mask    = "
    f"r < {SINK_MASK_RADIUS_AU:.6e} AU"
)
print(
    f"  Zoom range   = "
    f"{ZOOM_RADIUS_MIN_AU:.6e}–"
    f"{ZOOM_RADIUS_MAX_AU:.6e} AU"
)


# ============================================================
# VTKヘッダーからコード時刻を取得
# ============================================================
def read_vtk_time_code(filename):
    with open(filename, "rb") as vtk_file:
        header = vtk_file.read(512).decode(
            "ascii",
            errors="ignore",
        )

    match = re.search(
        r"time\s*=\s*"
        r"([+-]?(?:\d+\.?\d*|\.\d+)"
        r"(?:[eE][+-]?\d+)?)",
        header,
    )

    if match is None:
        raise ValueError(
            f"Could not find time in VTK header: "
            f"{filename}"
        )

    return float(match.group(1))


# ============================================================
# 対数カラースケール
# ============================================================
def calculate_log_limits(
    values,
    percentiles,
    fallback=None,
):
    values = np.asarray(values).ravel()

    positive_values = values[
        np.isfinite(values)
        & (values > 0.0)
    ]

    if len(positive_values) == 0:
        if fallback is not None:
            return fallback

        return 1.0e-17, 1.0e-12

    vmin, vmax = np.percentile(
        positive_values,
        percentiles,
    )

    if not np.isfinite(vmin) or vmin <= 0.0:
        vmin = np.min(positive_values)

    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = np.max(positive_values)

    if vmax <= vmin:
        vmin *= 0.5
        vmax *= 2.0

    # 色範囲が狭すぎる場合は最低1桁確保
    if vmax / vmin < 10.0:
        center = np.sqrt(vmin * vmax)

        vmin = center / np.sqrt(10.0)
        vmax = center * np.sqrt(10.0)

    return vmin, vmax


# ============================================================
# z=0に最も近いセル層を抽出
# ============================================================
def extract_xy_midplane(
    grid,
    density_name,
):
    """
    z=0と交差するAMRブロックから、
    z=0に最も近いセル中心面を1層だけ抽出する。

    Returns
    -------
    points_au
        (x, y, z)座標。単位はAU。
    density
        密度。単位はM_sun AU^-3。
    """
    bounds = grid.bounds

    # bounds = (xmin, xmax, ymin, ymax, zmin, zmax)
    if not (
        bounds[4] <= 0.0 <= bounds[5]
    ):
        return None, None

    points_code = (
        grid.cell_centers().points
    )

    if len(points_code) == 0:
        return None, None

    density_code = np.asarray(
        grid[density_name]
    ).reshape(-1)

    if len(density_code) != len(points_code):
        raise ValueError(
            "Density size does not match "
            "cell-center size."
        )

    z_values = points_code[:, 2]

    z_nearest = z_values[
        np.argmin(np.abs(z_values))
    ]

    tolerance = (
        1.0e-10
        * max(np.max(np.abs(z_values)), 1.0)
    )

    midplane_mask = np.isclose(
        z_values,
        z_nearest,
        rtol=1.0e-10,
        atol=tolerance,
    )

    if not np.any(midplane_mask):
        return None, None

    points_au = (
        points_code[midplane_mask]
        * LENGTH_UNIT_AU
    )

    density = (
        density_code[midplane_mask]
        * DENSITY_UNIT
    )

    return points_au, density


# ============================================================
# AMRブロック情報
# ============================================================
def extract_amr_block(grid):
    """
    z=0と交差するAMRブロックについて、
    x-y境界、AMRレベル、セル幅を取得する。
    """
    bounds = grid.bounds

    if not (
        bounds[4] <= 0.0 <= bounds[5]
    ):
        return None

    points_code = (
        grid.cell_centers().points
    )

    if len(points_code) == 0:
        return None

    x_values = np.unique(
        np.sort(points_code[:, 0])
    )

    if len(x_values) > 1:
        differences = np.diff(x_values)
        differences = differences[
            differences > 0.0
        ]

        if len(differences) > 0:
            dx_code = np.min(differences)
        else:
            dx_code = BASE_DX_CODE
    else:
        dx_code = BASE_DX_CODE

    dx_au = (
        dx_code * LENGTH_UNIT_AU
    )

    level = int(
        np.clip(
            np.round(
                np.log2(
                    BASE_DX_AU / dx_au
                )
            ),
            0,
            MAX_AMR_LEVEL,
        )
    )

    return {
        "bounds": (
            bounds[0] * LENGTH_UNIT_AU,
            bounds[1] * LENGTH_UNIT_AU,
            bounds[2] * LENGTH_UNIT_AU,
            bounds[3] * LENGTH_UNIT_AU,
        ),
        "level": level,
        "dx": dx_au,
    }


# ============================================================
# AMR境界描画
# ============================================================
def draw_amr_blocks(
    ax,
    blocks,
    linewidth=1.0,
):
    """
    境界線の重複を減らすため、
    各ブロックの右辺と上辺だけを描く。
    """
    sorted_blocks = sorted(
        blocks,
        key=lambda block: block["level"],
    )

    for block in sorted_blocks:
        x0, x1, y0, y1 = block["bounds"]
        level = block["level"]

        color = LEVEL_COLORS.get(
            level,
            "white",
        )

        # 右辺
        ax.plot(
            [x1, x1],
            [y0, y1],
            color=color,
            linewidth=linewidth,
            alpha=0.8,
        )

        # 上辺
        ax.plot(
            [x0, x1],
            [y1, y1],
            color=color,
            linewidth=linewidth,
            alpha=0.8,
        )


# ============================================================
# 線形補間＋最近傍補間
# ============================================================
def smooth_interpolation(
    points_2d,
    values,
    X,
    Y,
):
    """
    線形補間を基本とし、
    線形補間できない外縁だけを最近傍補間で埋める。
    """
    linear = griddata(
        points_2d,
        values,
        (X, Y),
        method="linear",
    )

    nearest = griddata(
        points_2d,
        values,
        (X, Y),
        method="nearest",
    )

    result = np.where(
        np.isfinite(linear),
        linear,
        nearest,
    )

    invalid = (
        ~np.isfinite(result)
        | (result <= 0.0)
    )

    result[invalid] = np.nan

    return result


# ============================================================
# VTKファイル整理
# ============================================================
if not os.path.isdir(vtk_dir):
    raise FileNotFoundError(
        f"VTK directory does not exist: {vtk_dir}"
    )

files_by_step = defaultdict(list)

for filename in os.listdir(vtk_dir):
    if not (
        filename.startswith("Toyouchi.block")
        and filename.endswith(".vtk")
    ):
        continue

    match = re.search(
        r"(?:prim\.)?out2\.(\d+)",
        filename,
    )

    if match is not None:
        timestep = int(match.group(1))

        files_by_step[timestep].append(
            os.path.join(
                vtk_dir,
                filename,
            )
        )

steps = sorted(files_by_step)

if not steps:
    raise RuntimeError(
        f"No VTK files found in {vtk_dir}"
    )

print(
    f"[INFO] Found {len(steps)} timesteps"
)
print(
    f"[INFO] Timestep range: "
    f"{steps[0]}–{steps[-1]}"
)


# ============================================================
# 密度変数名検出
# ============================================================
test_file = files_by_step[steps[0]][0]
test_grid = pv.read(test_file)

print(
    f"[INFO] Available arrays: "
    f"{test_grid.array_names}"
)

density_name = next(
    (
        name
        for name in [
            "dens",
            "density",
            "rho",
            "prim_dens",
            "prim_density",
        ]
        if name in test_grid.array_names
    ),
    None,
)

if density_name is None:
    raise RuntimeError(
        "No density array was found."
    )

print(
    f"[INFO] Density variable: "
    f"{density_name}"
)


# ============================================================
# 全体図共通カラースケール
# シンク内部を除外した密度から計算
# ============================================================
sample_indices = np.linspace(
    0,
    len(steps) - 1,
    min(10, len(steps)),
).astype(int)

sample_steps = [
    steps[index]
    for index in sample_indices
]

sample_density = []

for step in sample_steps:
    for filename in files_by_step[step]:
        try:
            grid = pv.read(filename)

            sample_points, density_sample = (
                extract_xy_midplane(
                    grid,
                    density_name,
                )
            )

            if sample_points is None:
                continue

            sample_radius_spherical = (
                np.linalg.norm(
                    sample_points,
                    axis=1,
                )
            )

            sample_outside_sink = (
                sample_radius_spherical
                >= SINK_MASK_RADIUS_AU
            )

            if np.any(sample_outside_sink):
                sample_density.append(
                    density_sample[
                        sample_outside_sink
                    ]
                )

        except Exception as error:
            print(
                f"[WARNING] Could not sample "
                f"{filename}: {error}"
            )

if not sample_density:
    raise RuntimeError(
        "No valid density data outside "
        "the sink were found."
    )

sample_density = np.concatenate(
    sample_density
)

full_vmin, full_vmax = (
    calculate_log_limits(
        sample_density,
        FULL_PERCENTILES,
    )
)

full_norm = LogNorm(
    vmin=full_vmin,
    vmax=full_vmax,
)

print(
    f"[INFO] Global density range "
    f"(outside sink): "
    f"{full_vmin:.3e}–"
    f"{full_vmax:.3e} M_sun AU^-3"
)


# ============================================================
# メイン処理
# ============================================================
for timestep_index, step in enumerate(steps):
    print(
        f"[INFO] Processing timestep "
        f"{step:05d} "
        f"({timestep_index + 1}/"
        f"{len(steps)})"
    )

    point_arrays = []
    density_arrays = []
    amr_blocks = []

    for filename in files_by_step[step]:
        try:
            grid = pv.read(filename)

            points, density = (
                extract_xy_midplane(
                    grid,
                    density_name,
                )
            )

            block = extract_amr_block(grid)

            if points is not None:
                point_arrays.append(points)
                density_arrays.append(density)

            if block is not None:
                amr_blocks.append(block)

        except Exception as error:
            print(
                f"[WARNING] Failed to process "
                f"{filename}: {error}"
            )

    if not point_arrays:
        print(
            f"[WARNING] No x-y midplane data "
            f"at timestep {step:05d}"
        )
        continue

    points = np.vstack(point_arrays)
    density = np.hstack(density_arrays)

    valid = (
        np.all(
            np.isfinite(points),
            axis=1,
        )
        & np.isfinite(density)
        & (density > 0.0)
    )

    points = points[valid]
    density = density[valid]

    if len(density) == 0:
        print(
            f"[WARNING] No valid density data "
            f"at timestep {step:05d}"
        )
        continue

    # ========================================================
    # 時刻
    # ========================================================
    representative_file = (
        files_by_step[step][0]
    )

    time_code = read_vtk_time_code(
        representative_file
    )

    time_yr = (
        time_code * TIME_UNIT_YR
    )

    # ========================================================
    # 平面内半径と球半径
    # ========================================================
    radius_xy = np.hypot(
        points[:, 0],
        points[:, 1],
    )

    # シンク判定には球半径を使用
    radius_spherical = np.linalg.norm(
        points,
        axis=1,
    )

    # ========================================================
    # ズーム半径
    # ========================================================
    nearest_cell_radius = np.min(
        radius_xy
    )

    zoom_radius = np.clip(
        ZOOM_RADIUS_FACTOR
        * nearest_cell_radius,
        ZOOM_RADIUS_MIN_AU,
        ZOOM_RADIUS_MAX_AU,
    )

    # 実データ領域を超えないように制限
    data_half_width = min(
        np.max(np.abs(points[:, 0])),
        np.max(np.abs(points[:, 1])),
    )

    zoom_radius = min(
        zoom_radius,
        data_half_width,
    )

    print(
        f"[INFO] Zoom radius = "
        f"{zoom_radius:.3e} AU"
    )

    # ========================================================
    # 全体図補間
    # ========================================================
    axis_full = np.linspace(
        DOMAIN_MIN,
        DOMAIN_MAX,
        FULL_RESOLUTION,
    )

    X_full, Y_full = np.meshgrid(
        axis_full,
        axis_full,
    )

    density_full = smooth_interpolation(
        points[:, [0, 1]],
        density,
        X_full,
        Y_full,
    )

    # ========================================================
    # ズーム領域
    # ========================================================
    zoom_mask = (
        radius_xy <= zoom_radius
    )

    points_zoom = points[zoom_mask]
    density_zoom_raw = density[zoom_mask]

    radius_spherical_zoom = (
        radius_spherical[zoom_mask]
    )

    if len(points_zoom) <= 10:
        print(
            f"[WARNING] Insufficient zoom data "
            f"at timestep {step:05d}"
        )
        continue

    # ========================================================
    # シンク外部を選ぶ診断用マスク
    # ========================================================
    outside_sink_mask = (
        radius_spherical_zoom
        >= SINK_MASK_RADIUS_AU
    )

    if not np.any(outside_sink_mask):
        print(
            f"[WARNING] No cells outside sink "
            f"at timestep {step:05d}"
        )
        continue

    points_diagnostic = points_zoom[
        outside_sink_mask
    ]

    density_diagnostic = density_zoom_raw[
        outside_sink_mask
    ]

    radius_diagnostic = (
        radius_spherical_zoom[
            outside_sink_mask
        ]
    )

    # ========================================================
    # ズーム補間
    # 描画にはシンク内部も含める
    # ========================================================
    axis_zoom = np.linspace(
        -zoom_radius,
        zoom_radius,
        ZOOM_RESOLUTION,
    )

    X_zoom, Y_zoom = np.meshgrid(
        axis_zoom,
        axis_zoom,
    )

    density_zoom = smooth_interpolation(
        points_zoom[:, [0, 1]],
        density_zoom_raw,
        X_zoom,
        Y_zoom,
    )

    # ========================================================
    # ズーム専用カラースケール
    # シンク外部の密度のみから計算
    # ========================================================
    zoom_vmin, zoom_vmax = (
        calculate_log_limits(
            density_diagnostic,
            ZOOM_PERCENTILES,
            fallback=(
                full_vmin,
                full_vmax,
            ),
        )
    )

    zoom_norm = LogNorm(
        vmin=zoom_vmin,
        vmax=zoom_vmax,
    )

    # ========================================================
    # 最大・最小密度
    # シンク外かつズーム領域内だけを検索
    # ========================================================
    maximum_index = np.argmax(
        density_diagnostic
    )

    minimum_index = np.argmin(
        density_diagnostic
    )

    rho_max = density_diagnostic[
        maximum_index
    ]

    rho_min = density_diagnostic[
        minimum_index
    ]

    r_max = radius_diagnostic[
        maximum_index
    ]

    r_min = radius_diagnostic[
        minimum_index
    ]

    position_max = points_diagnostic[
        maximum_index
    ]

    position_min = points_diagnostic[
        minimum_index
    ]

    print(
        f"[INFO] Density extrema outside sink "
        f"(r >= {SINK_MASK_RADIUS_AU:.3e} AU):"
    )

    print(
        f"  rho_max = {rho_max:.3e} "
        f"M_sun AU^-3 at "
        f"({position_max[0]:.3e}, "
        f"{position_max[1]:.3e}, "
        f"{position_max[2]:.3e}) AU, "
        f"r = {r_max:.3e} AU"
    )

    print(
        f"  rho_min = {rho_min:.3e} "
        f"M_sun AU^-3 at "
        f"({position_min[0]:.3e}, "
        f"{position_min[1]:.3e}, "
        f"{position_min[2]:.3e}) AU, "
        f"r = {r_min:.3e} AU"
    )

    # ========================================================
    # Figure
    # ========================================================
    fig = plt.figure(figsize=figsize)

    gs = GridSpec(
        1,
        5,
        width_ratios=[
            4.0,
            0.25,
            4.0,
            0.25,
            1.6,
        ],
        wspace=0.35,
    )

    # ========================================================
    # 全体図
    # ========================================================
    ax_full = fig.add_subplot(gs[0, 0])
    ax_full.set_aspect("equal")

    image_full = ax_full.pcolormesh(
        X_full,
        Y_full,
        density_full,
        cmap=FULL_CMAP,
        norm=full_norm,
        shading="auto",
        rasterized=True,
    )

    draw_amr_blocks(
        ax_full,
        amr_blocks,
        linewidth=1.0,
    )

    # ズーム範囲
    zoom_rectangle = patches.Rectangle(
        (
            -zoom_radius,
            -zoom_radius,
        ),
        2.0 * zoom_radius,
        2.0 * zoom_radius,
        edgecolor="cyan",
        facecolor="none",
        linestyle="--",
        linewidth=2.0,
    )

    ax_full.add_patch(
        zoom_rectangle
    )

    ax_full.plot(
        0.0,
        0.0,
        "r+",
        markersize=12,
        markeredgewidth=2,
    )

    ax_full.set_xlim(
        DOMAIN_MIN,
        DOMAIN_MAX,
    )

    ax_full.set_ylim(
        DOMAIN_MIN,
        DOMAIN_MAX,
    )

    ax_full.set_xlabel(
        r"$x\ [{\rm AU}]$",
        fontsize=12,
    )

    ax_full.set_ylabel(
        r"$y\ [{\rm AU}]$",
        fontsize=12,
    )

    ax_full.set_title(
        "Full Domain",
        fontsize=12,
        fontweight="bold",
    )

    ax_full.ticklabel_format(
        axis="both",
        style="scientific",
        scilimits=(0, 0),
    )

    ax_full.grid(False)

    # AMR凡例
    legend_elements = []

    for level in range(
        MAX_AMR_LEVEL + 1
    ):
        dx_au = (
            BASE_DX_AU / 2**level
        )

        legend_elements.append(
            mpatches.Patch(
                facecolor="none",
                edgecolor=(
                    LEVEL_COLORS[level]
                ),
                linewidth=1.5,
                label=(
                    f"Level {level}, "
                    f"dx={dx_au:.2e} AU"
                ),
            )
        )

    legend_elements.append(
        mpatches.Patch(
            facecolor="none",
            edgecolor="cyan",
            linestyle="--",
            linewidth=1.5,
            label="Zoom region",
        )
    )

    ax_full.legend(
        handles=legend_elements,
        fontsize=7,
        loc="upper right",
        framealpha=0.75,
    )

    # ========================================================
    # 全体図カラーバー
    # ========================================================
    cax_full = fig.add_subplot(gs[0, 1])

    cbar_full = fig.colorbar(
        image_full,
        cax=cax_full,
        extend="both",
    )

    cbar_full.ax.yaxis.set_major_formatter(
        LogFormatterSciNotation()
    )

    cbar_full.set_label(
        r"$\rho\ "
        r"[M_\odot\,{\rm AU}^{-3}]$"
        "\nGlobal scale",
        fontsize=10,
    )

    # ========================================================
    # ズーム図
    # ========================================================
    ax_zoom = fig.add_subplot(gs[0, 2])
    ax_zoom.set_aspect("equal")

    image_zoom = ax_zoom.imshow(
        density_zoom,
        origin="lower",
        extent=[
            -zoom_radius,
            zoom_radius,
            -zoom_radius,
            zoom_radius,
        ],
        cmap=ZOOM_CMAP,
        norm=zoom_norm,
        interpolation=ZOOM_INTERPOLATION,
        aspect="equal",
        rasterized=True,
    )

    # ズーム範囲と交差するAMRブロックを抽出
    zoom_blocks = []

    for block in amr_blocks:
        x0, x1, y0, y1 = block["bounds"]

        intersects_zoom = (
            x0 <= zoom_radius
            and x1 >= -zoom_radius
            and y0 <= zoom_radius
            and y1 >= -zoom_radius
        )

        if not intersects_zoom:
            continue

        clipped_x0 = max(
            x0,
            -zoom_radius,
        )
        clipped_x1 = min(
            x1,
            zoom_radius,
        )
        clipped_y0 = max(
            y0,
            -zoom_radius,
        )
        clipped_y1 = min(
            y1,
            zoom_radius,
        )

        if (
            clipped_x1 > clipped_x0
            and clipped_y1 > clipped_y0
        ):
            zoom_blocks.append(
                {
                    "bounds": (
                        clipped_x0,
                        clipped_x1,
                        clipped_y0,
                        clipped_y1,
                    ),
                    "level": block["level"],
                }
            )

    draw_amr_blocks(
        ax_zoom,
        zoom_blocks,
        linewidth=1.0,
    )

    ax_zoom.plot(
        0.0,
        0.0,
        "r+",
        markersize=12,
        markeredgewidth=2,
        label="Center",
    )

    ax_zoom.set_xlim(
        -zoom_radius,
        zoom_radius,
    )

    ax_zoom.set_ylim(
        -zoom_radius,
        zoom_radius,
    )

    ax_zoom.set_xlabel(
        r"$x\ [{\rm AU}]$",
        fontsize=12,
    )

    ax_zoom.set_ylabel(
        r"$y\ [{\rm AU}]$",
        fontsize=12,
    )

    ax_zoom.set_title(
        "Zoomed Region\n"
        f"R = {zoom_radius:.3e} AU",
        fontsize=12,
        fontweight="bold",
    )

    ax_zoom.ticklabel_format(
        axis="both",
        style="scientific",
        scilimits=(0, 0),
    )

    ax_zoom.grid(False)

    ax_zoom.legend(
        loc="upper right",
        fontsize=8,
        framealpha=0.7,
    )

    # ========================================================
    # ズーム図カラーバー
    # ========================================================
    cax_zoom = fig.add_subplot(gs[0, 3])

    cbar_zoom = fig.colorbar(
        image_zoom,
        cax=cax_zoom,
        extend="both",
    )

    cbar_zoom.ax.yaxis.set_major_formatter(
        LogFormatterSciNotation()
    )

    cbar_zoom.set_label(
        r"$\rho\ "
        r"[M_\odot\,{\rm AU}^{-3}]$"
        "\n"
        f"Zoom {ZOOM_PERCENTILES[0]:.0f}–"
        f"{ZOOM_PERCENTILES[1]:.0f} percentile"
        "\n(outside sink)",
        fontsize=10,
    )

    # ========================================================
    # 情報欄
    # ========================================================
    ax_info = fig.add_subplot(gs[0, 4])
    ax_info.axis("off")

    information_text = (
        "Time\n"
        f"{time_yr:.6e} yr\n\n"

        "Code time\n"
        f"{time_code:.6e}\n\n"

        "Zoom radius\n"
        f"{zoom_radius:.3e} AU\n\n"

        "Density search region\n"
        f"r >= {SINK_MASK_RADIUS_AU:.3e} AU\n"
        "(outside sink)\n\n"

        "Maximum density\n"
        f"{rho_max:.3e}\n"
        r"$M_\odot\,{\rm AU}^{-3}$"
        "\n"
        f"r = {r_max:.3e} AU\n\n"

        "Minimum density\n"
        f"{rho_min:.3e}\n"
        r"$M_\odot\,{\rm AU}^{-3}$"
        "\n"
        f"r = {r_min:.3e} AU\n\n"

        "Zoom color range\n"
        f"{zoom_vmin:.2e}\n"
        "to\n"
        f"{zoom_vmax:.2e}"
    )

    ax_info.text(
        0.04,
        0.95,
        information_text,
        transform=ax_info.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round",
            facecolor="whitesmoke",
            edgecolor="black",
            alpha=0.9,
        ),
    )

    # ========================================================
    # タイトル
    # ========================================================
    fig.suptitle(
        "AMR Density Map: x-y Midplane\n"
        f"t = {time_yr:.6e} yr | "
        f"Variable: {density_name}",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    plt.subplots_adjust(
        top=0.88,
        bottom=0.10,
        left=0.05,
        right=0.97,
    )

    # ========================================================
    # タイムステップ順のファイル名で保存
    # ========================================================
    png = os.path.join(
        output_dir,
        (
            f"density_xy_timestep_"
            f"{step:05d}_"
            f"time_{time_yr:.6e}yr_"
            "amr_blocks.png"
        ),
    )

    fig.savefig(
        png,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    print(
        f"[INFO] Saved: "
        f"step={step:05d}, "
        f"time={time_yr:.6e} yr"
    )


print(
    f"[INFO] All AMR density maps "
    f"saved to: {output_dir}"
)

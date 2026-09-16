# -*- coding: utf-8 -*-
# Converted from 質量.ipynb

# %% cell 1
# 質量保存を CGS 単位で確認するコード

import os

import matplotlib.pyplot as plt
import numpy as np


# ========= 入力ファイル =========
hst_file = os.path.expanduser(
    "~/athena-project/results/〇〇/Toyouchi.hst"
)

# ========= 出力ディレクトリ =========
output_dir = "./total_mass"
os.makedirs(output_dir, exist_ok=True)

# ========= 単位定義 =========
# Toyouchi.cppのコード単位
T_UNIT_CGS = 3.61e10  # 1 code time [s] (L0 = rb / 2)
M_UNIT_CGS = 4.0e33   # 1 code mass [g]

# 物理定数
YEAR_CGS = 365.25 * 24.0 * 3600.0  # s
MYR_CGS = 1.0e6 * YEAR_CGS         # s
MSUN_CGS = 1.98847e33               # g

# コード単位 → 表示単位
T_UNIT_MYR = T_UNIT_CGS / MYR_CGS
M_UNIT_MSUN = M_UNIT_CGS / MSUN_CGS

print("[INFO] Unit conversion factors:")
print(
    f"  1 code time = {T_UNIT_MYR:.6e} Myr"
)
print(
    f"  1 code mass = {M_UNIT_MSUN:.6e} M_sun"
)

# ========= データ読み込み =========
# コメント行（#）を無視
data = np.loadtxt(
    hst_file,
    comments="#",
)

if data.ndim == 1:
    data = data[np.newaxis, :]

if data.shape[1] < 3:
    raise ValueError(
        "The history file has fewer than three columns."
    )

time_code = data[:, 0]
mass_code = data[:, 2]


# ========= コード単位から変換 =========
time_myr = (
    time_code * T_UNIT_MYR
)

mass_msun = (
    mass_code * M_UNIT_MSUN
)


# ========= 質量変化率 =========
mass_init_msun = mass_msun[0]
mass_final_msun = mass_msun[-1]

if (
    not np.isfinite(mass_init_msun)
    or mass_init_msun == 0.0
):
    raise ValueError(
        "Initial mass is zero or non-finite, "
        "so the mass-change rate cannot be calculated."
    )

mass_change = (
    (mass_final_msun - mass_init_msun)
    / mass_init_msun
    * 100.0
)

print(
    f"Initial mass = "
    f"{mass_init_msun:.6e} M_sun"
)
print(
    f"Final mass   = "
    f"{mass_final_msun:.6e} M_sun"
)
print(
    f"Mass change  = "
    f"{mass_change:+.3f} %"
)


# ========= mass消失チェック =========
# 元コードの1e-10 code massをM_sunへ変換
mass_threshold_msun = (
    1.0e-10 * M_UNIT_MSUN
)

bad = (
    (mass_msun < mass_threshold_msun)
    | (~np.isfinite(mass_msun))
)

if np.any(bad):
    print(
        "⚠️ massが異常になったタイムステップ:"
    )

    for t_myr, m_msun in zip(
        time_myr[bad],
        mass_msun[bad],
    ):
        print(
            f"  t = {t_myr:.6e} Myr, "
            f"mass = {m_msun:.6e} M_sun"
        )
else:
    print(
        "✅ massは全ステップで正常"
    )


# ========= プロット =========
fig, ax = plt.subplots(
    figsize=(6, 4),
)

ax.plot(
    time_myr,
    mass_msun,
    marker="o",
    linestyle="-",
    markersize=3,
)

ax.set_xlabel(
    "Time [Myr]"
)
ax.set_ylabel(
    r"Total mass [$M_\odot$]"
)
ax.set_yscale(
    "log"
)
ax.set_title(
    "Total Mass vs Time"
)
ax.grid(
    True,
    which="both",
    linestyle="--",
    alpha=0.5,
)

fig.text(
    0.5,
    0.02,
    (
        f"Mass change = {mass_change:+.3f}% "
        f"(Initial = {mass_init_msun:.3e} "
        r"$M_\odot$"
        f", Final = {mass_final_msun:.3e} "
        r"$M_\odot$"
        ")"
    ),
    ha="center",
    fontsize=10,
    bbox=dict(
        boxstyle="round",
        facecolor="whitesmoke",
        edgecolor="black",
        alpha=0.8,
    ),
)

fig.tight_layout(
    rect=[0, 0.07, 1, 1]
)


# ========= 保存 =========
output_file = os.path.join(
    output_dir,
    "mass_vs_time_msun_myr.png",
)

fig.savefig(
    output_file,
    dpi=200,
    bbox_inches="tight",
)

print(
    f"[INFO] Saved: {output_file}"
)

plt.show()

# %% cell 2
#各タイムステップで計算領域の質量保存（dM/dt+flux=0）を確認するコード

import pyvista as pv
import numpy as np
import glob
import re
import os

# =====================================================
# 設定
# =====================================================

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")

xmin,xmax = 0.0,8.0  #計算領域に合わせる
ymin,ymax = 0.0,8.0
zmin,zmax = 0.0,8.0

# =====================================================
# VTK ファイル取得
# =====================================================

vtk_files = sorted(glob.glob(os.path.join(vtk_dir,"Jeans.block*.out2.*.vtk")))

print("Total VTK files:",len(vtk_files))

# timestep 抽出
steps = sorted(set(re.search(r"out2\.(\d+)",f).group(1) for f in vtk_files))

print("Detected timesteps:",len(steps))

# =====================================================
# 保存用配列
# =====================================================

time_list = []
mass_list = []
flux_list = []

# =====================================================
# timestep ループ
# =====================================================

for step in steps:

    files = sorted(glob.glob(os.path.join(vtk_dir,f"Jeans.block*.out2.{step}.vtk")))

    total_mass = 0.0
    total_flux = 0.0

    for f in files:

        grid = pv.read(f)

        rho = grid.cell_data["dens"]
        mom = grid.cell_data["mom"]

        vx = mom[:,0] / rho
        vy = mom[:,1] / rho
        vz = mom[:,2] / rho

        centers = grid.cell_centers().points

        x = centers[:,0]
        y = centers[:,1]
        z = centers[:,2]

        bounds = grid.bounds

        dx = (bounds[1]-bounds[0]) / grid.dimensions[0]
        dy = (bounds[3]-bounds[2]) / grid.dimensions[1]
        dz = (bounds[5]-bounds[4]) / grid.dimensions[2]

        cell_vol = dx*dy*dz

        # ------------------------------------------------
        # total mass
        # ------------------------------------------------

        total_mass += np.sum(rho * cell_vol)

        # ------------------------------------------------
        # flux
        # ------------------------------------------------

        mask = np.isclose(x,xmin)
        total_flux += np.sum(rho[mask]*vx[mask]*dy*dz)

        mask = np.isclose(x,xmax)
        total_flux += np.sum(rho[mask]*vx[mask]*dy*dz)

        mask = np.isclose(y,ymin)
        total_flux += np.sum(rho[mask]*vy[mask]*dx*dz)

        mask = np.isclose(y,ymax)
        total_flux += np.sum(rho[mask]*vy[mask]*dx*dz)

        mask = np.isclose(z,zmin)
        total_flux += np.sum(rho[mask]*vz[mask]*dx*dy)

        mask = np.isclose(z,zmax)
        total_flux += np.sum(rho[mask]*vz[mask]*dx*dy)

    time_list.append(int(step))
    mass_list.append(total_mass)
    flux_list.append(total_flux)

# =====================================================
# dM/dt 計算
# =====================================================

mass_list = np.array(mass_list)
flux_list = np.array(flux_list)
time_list = np.array(time_list)

dMdt = np.zeros_like(mass_list)

for i in range(1,len(mass_list)):

    dMdt[i] = mass_list[i] - mass_list[i-1]

# =====================================================
# 結果表示
# =====================================================

print("\n=== Mass evolution ===")

for i in range(len(time_list)):

    print(f"step {time_list[i]:5d}  M={mass_list[i]:.6e}  flux={flux_list[i]:.6e}  dM={dMdt[i]:.6e}")

# =====================================================
# 保存則チェック
# =====================================================

print("\n=== Conservation check ===")

error = dMdt + flux_list

print("mean error:",np.mean(np.abs(error)))
print("max  error:",np.max(np.abs(error)))

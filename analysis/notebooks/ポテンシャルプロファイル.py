# -*- coding: utf-8 -*-
# Converted from ポテンシャルプロファイル.ipynb

# %% cell 1
#重力境界条件が正しく機能しているかをチェックするコード
#multipoleの場合、遠方で∝-1/rになればOK

import pyvista as pv
import numpy as np
import glob
import os
import re
import matplotlib.pyplot as plt

# ============================================
# directory
# ============================================

vtk_dir = os.path.expanduser("~/athenapp/results/〇〇")

vtk_files = sorted(glob.glob(os.path.join(vtk_dir,"Jeans.block*.out2.00000.vtk")))

print("files:",len(vtk_files))

# ============================================
# center
# ============================================

xc = 4.0  #設定により変更
yc = 4.0
zc = 4.0

# ============================================
# storage
# ============================================

r_all = []
phi_all = []

# ============================================
# read blocks
# ============================================

for f in vtk_files:

    grid = pv.read(f)

    phi = grid.cell_data["phi"]

    centers = grid.cell_centers().points

    x = centers[:,0]
    y = centers[:,1]
    z = centers[:,2]

    r = np.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)

    r_all.append(r)
    phi_all.append(phi)

r_all = np.concatenate(r_all)
phi_all = np.concatenate(phi_all)

# ============================================
# sort
# ============================================

idx = np.argsort(r_all)

r_all = r_all[idx]
phi_all = phi_all[idx]

# ============================================
# plot
# ============================================

plt.figure(figsize=(6,5))

plt.scatter(r_all,phi_all,s=1,label="simulation")

# 1/r reference
r = np.linspace(0.1,6,200)
plt.plot(r,-1/r,"r",label="~1/r")

plt.xlabel("r")
plt.ylabel("phi")

plt.legend()
plt.tight_layout()
plt.show()

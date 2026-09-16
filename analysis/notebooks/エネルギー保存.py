# -*- coding: utf-8 -*-
# Converted from エネルギー保存.ipynb

# %% cell 1
#hstファイルを読み各エネルギーの保存を確認するコード

import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================
# hst ファイル読み込み
# ============================================

file = "/home/aian/athenapp/results/〇〇/〇〇.hst"

data = np.loadtxt(file)

time = data[:,0]
dt   = data[:,1]
mass = data[:,2]

momx = data[:,3]
momy = data[:,4]
momz = data[:,5]

KE1  = data[:,6]
KE2  = data[:,7]
KE3  = data[:,8]

Etot = data[:,9]
Egrav = data[:,10]

# 出力
output_dir = "./energy_profiles"
os.makedirs(output_dir, exist_ok=True)

# ============================================
# エネルギー計算
# ============================================

KE = KE1 + KE2 + KE3

#内部エネルギー
Eint = Etot - KE

# 真の全エネルギー
E_total = Etot + Egrav

# 初期値
E0 = E_total[0]

# 相対誤差
E_err = (E_total - E0) / abs(E0)  #初期エネルギーからのズレ（数値計算によるエネルギー非保存の指標）

# ============================================
# 基本情報
# ============================================

print("Initial mass:", mass[0])
print("Final mass:", mass[-1])

print("\nInitial total energy:", E_total[0])
print("Final total energy:", E_total[-1])

print("\nRelative energy error:", E_err[-1])

# ============================================
# 図1：主要エネルギー
# ============================================

plt.figure(figsize=(10,6))

plt.plot(time, Egrav, label="Gravitational Energy")
plt.plot(time, Eint, label="Internal Energy")
plt.plot(time, E_total, label="Total Energy")

plt.xlabel("Time(kyr)")
plt.ylabel("Energy")
plt.xlim(0, max(time))

plt.legend()
plt.grid()

plt.tight_layout()
png = os.path.join(output_dir, f"main_energy.png")
plt.savefig(png, dpi=300)

plt.show()

# ============================================
# 図2：運動エネルギー
# ============================================

plt.figure(figsize=(10,5))

plt.plot(time, KE)

plt.xlabel("Time(kyr)")
plt.ylabel("Kinetic Energy")
plt.xlim(0, max(time))

plt.yscale("log")

plt.grid()
plt.tight_layout()
png = os.path.join(output_dir, f"kinetic_energy.png")
plt.savefig(png, dpi=300)

plt.show()

# ============================================
# 図3：Virial 指標
# ============================================

virial = KE / np.abs(Egrav)

plt.figure(figsize=(10,5))

plt.plot(time, virial)

plt.xlabel("Time(kyr)")
plt.ylabel("KE / |Egrav|")
plt.xlim(0, max(time))

plt.yscale("log")

plt.grid()
plt.tight_layout()
png = os.path.join(output_dir, f"virial_index.png")
plt.savefig(png, dpi=300)

plt.show()

# ============================================
# エネルギー誤差
# ============================================

plt.figure(figsize=(10,5))

plt.plot(time, E_err)

plt.xlabel("Time(kyr)")
plt.xlim(0,max(time))
plt.ylabel("Relative Energy Error")

plt.grid()
plt.tight_layout()
png = os.path.join(output_dir, f"relative_energy_error.png")
plt.savefig(png, dpi=300)

plt.show()


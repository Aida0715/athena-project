# -*- coding: utf-8 -*-
# Converted from 動画結合.ipynb

# %% cell 1
#個別の2つの動画を結合し、横並びで再生する新しい動画を作成するコード

import os
import subprocess

# ============================
# 入力動画
# ============================
video1 = "/home/aian/athena-project/results/jeans_3d-test26/xy_vector_field/xy_vector_field.mp4"
video2 = "/home/aian/athena-project/results/jeans_3d-test26/xy_density_maps/xy_density_maps_movie.mp4"

# ============================
# 出力先
# ============================
outdir = "/home/aian/athena-project/results/jeans_3d-test26/videos"
os.makedirs(outdir, exist_ok=True)

outfile = os.path.join(outdir, "〇〇.mp4")

# ============================
# ffmpeg コマンド
# ============================
cmd = [
    "ffmpeg",
    "-y",
    "-i", video1,
    "-i", video2,
    "-filter_complex", "hstack=shortest=1",
    outfile
]

# ============================
# 実行
# ============================
subprocess.run(cmd, check=True)

print("Saved to:", outfile)


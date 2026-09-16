# -*- coding: utf-8 -*-
# Converted from 動画.ipynb

# %% cell 1
#画像を紙芝居形式で動画化して保存する

import matplotlib.animation as animation
from PIL import Image
import natsort
import matplotlib.pyplot as plt
import os
import natsort

# ============================================================
# 設定
# ============================================================
# 画像ディレクトリ
img_dir = os.path.expanduser("~/athena-project/results/〇〇/〇〇")
# 出力動画ファイル名
output_video = os.path.join(img_dir, "〇〇_movie.mp4")
# フレーム間隔（秒）
interval = 1000  # 100ms → 10fps

# ============================================================
# 画像一覧を取得
# ============================================================
img_files = natsort.natsorted([os.path.join(img_dir, f) for f in os.listdir(img_dir) if f.endswith(".png")])
if len(img_files) == 0:
    raise FileNotFoundError("画像が見つかりません。パスを確認してください。")

# ============================================================
# Figure 作成
# ============================================================
fig, ax = plt.subplots(figsize=(6,6))
ax.axis('off')  # 軸は非表示
im = ax.imshow(Image.open(img_files[0]))

# ============================================================
# 更新関数
# ============================================================
def update(frame):
    im.set_array(Image.open(img_files[frame]))
    return [im]

# ============================================================
# アニメーション作成
# ============================================================
ani = animation.FuncAnimation(fig, update, frames=len(img_files), interval=interval, blit=True)

# ============================================================
# 動画保存
# ============================================================
ani.save(output_video, writer="ffmpeg", dpi=200)
print(f"✅ 動画を保存しました: {output_video}")


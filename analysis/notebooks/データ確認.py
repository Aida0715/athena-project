# -*- coding: utf-8 -*-
# Converted from データ確認.ipynb

# %% cell 1
# VTK格納データの変数名を検索
import pyvista as pv
import os
from pathlib import Path

vtk_dir = Path(
    os.path.expanduser(
        "~/athena-project/results/Toyouchi-test27"
    )
).resolve()

sample_file = sorted(vtk_dir.glob("*out2*.vtk"))[0]
sample_grid = pv.read(sample_file)

print("cell_data :", list(sample_grid.cell_data.keys()))
print("point_data:", list(sample_grid.point_data.keys()))
print("field_data:", list(sample_grid.field_data.keys()))


#!/usr/bin/env python3
"""
Convert bitmap to perfect grayscale, map values to 1–1000 (permeability),
and display using a colorful colormap with pump and drain locations.
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

input_bitmap = "./permeability.bmp"
output_grayscale = "permeability_grayscale.bmp"
output_scaled = "permeability_scaled.npy"

img = Image.open(input_bitmap).convert("L")
gray = np.array(img, dtype=np.float32)

min_val = gray.min()
max_val = gray.max()
permeability = 1 + (gray - min_val) * (999 / (max_val - min_val))

sources = [(0.4, 0.85), (0.9, 0.2)]
sinks   = [(0.9, 0.80), (0.25, 0.25)]

H, W = gray.shape

def to_pixel_coords(pt):
    x, y = pt
    return (x * W, (1 - y) * H)

sources_px = [to_pixel_coords(p) for p in sources]
sinks_px   = [to_pixel_coords(p) for p in sinks]

plt.figure(figsize=(8, 6))
plt.imshow(permeability, cmap="plasma")
plt.colorbar(label="Permeability")
plt.title("Permeability Map with pumps and drains")
plt.axis("off")

# Draw sources (injection wells)
for i, (x, y) in enumerate(sources_px):
    if i == 0:
        plt.scatter(x, y, s=120, facecolors='none', edgecolors='white',
                    linewidths=4, label="Source (Pump)")
    else:
        plt.scatter(x, y, s=120, facecolors='none', edgecolors='white',
                    linewidths=4)

# Draw sinks (extraction wells)
for i, (x, y) in enumerate(sinks_px):
    if i == 0:
        plt.scatter(x, y, s=160, marker='x', color='white',
                    linewidths=4, label="Sink (Drain)")
    else:
        plt.scatter(x, y, s=160, marker='x', color='white', linewidths=4)

# Add legend (in upper left corner)
plt.legend(loc="upper left", facecolor="black", labelcolor="white")

# Save result
plt.savefig("permeability_with_pumps.png", dpi=200, bbox_inches="tight")

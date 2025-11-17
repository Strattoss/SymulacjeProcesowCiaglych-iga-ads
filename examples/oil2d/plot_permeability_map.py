#!/usr/bin/env python3
"""
Convert bitmap to perfect grayscale, map values to 1–1000 (permeability),
and display using a colorful colormap.
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

input_bitmap = "./permeability.bmp"
output_grayscale = "permeability_grayscale.bmp"
output_scaled = "permeability_scaled.npy"

img = Image.open(input_bitmap).convert("L")
# img.save(output_grayscale)

gray = np.array(img, dtype=np.float32)

min_val = gray.min()
max_val = gray.max()

permeability = 1 + (gray - min_val) * (999 / (max_val - min_val))

# np.save(output_scaled, permeability)

plt.figure(figsize=(8, 6))
plt.imshow(permeability, cmap="plasma")
plt.colorbar(label="Permeability map")
plt.title("Permeability Map")
plt.axis("off")
plt.savefig("permeability_plot.png")

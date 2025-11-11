import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os

# Parametry
size = 20
timesteps = range(0, 9990, 10)

os.makedirs("wind_maps", exist_ok=True)

def wind_field(X, Y, t):
    speed = 40.0

    # Angle smoothly transitions from 0 (west) to π/2 (south)
    angle = (np.pi / 2.0) * (1.0 - np.exp(-t / 1500.0))

    # Compute vector components — broadcast over the arrays
    Bx = -speed * np.cos(angle) * np.ones_like(X)
    By = -speed * np.sin(angle) * np.ones_like(Y)

    return Bx, By

# Siatka
x = np.linspace(0, 100, size)
y = np.linspace(0, 100, size)
X, Y = np.meshgrid(x, y)

for t in timesteps:
    png_name = f"wind_maps/wind_field_t{t:05}.png"

    Bx, By = wind_field(X, Y, t)
    magnitude = np.sqrt(Bx**2 + By**2)

    plt.figure(figsize=(6, 5))
    Q = plt.quiver(X, Y, Bx, By, magnitude,
                   cmap="coolwarm", scale=150, angles="xy")
    plt.colorbar(Q, label="Wind magnitude")
    plt.title(
        f"Wind Field at t={t}\n"
        f"mean(Bx)={Bx.mean():.2f}, mean(By)={By.mean():.2f}, mean(|W|)={magnitude.mean():.2f}"
    )
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.tight_layout()
    plt.savefig(png_name, dpi=200)
    plt.close()
    print(f"Generated {png_name}")

print("✅ Wygenerowano mapy wiatru i wizualizacje w folderze 'wind_maps/'")

# Generate gif with all the wind plots
import glob

print("Generating gif...")

# Collect all PNGs in sorted order
png_files = sorted(glob.glob("wind_maps/wind_field_t*.png"))

frames = [Image.open(png) for png in png_files]

if frames:
    frames[0].save(
        "wind_maps/wind_field_animation.gif",
        save_all=True,
        append_images=frames[1:],
        duration=10,
        loop=0 # infinite loop
    )
    print("GIF saved as 'wind_maps/wind_field_animation.gif'")
else:
    print("⚠️ No PNG files found in wind_maps/")
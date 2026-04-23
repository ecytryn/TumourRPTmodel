import numpy as np
import matplotlib.pyplot as plt

# Load both CSV files
oxygen_pure = np.loadtxt('scripts/oxygen_field_day47.csv', delimiter=',')
oxygen_with_tumor = np.loadtxt('scripts/day_047_oxygen.csv', delimiter=',')

# Check if identical
print(f"Arrays identical: {np.allclose(oxygen_pure, oxygen_with_tumor)}")
print(f"Max difference: {np.max(np.abs(oxygen_pure - oxygen_with_tumor))}")

# Visualize both with identical normalization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
vmin, vmax = 0, 15000

axes[0].imshow(oxygen_pure, vmin=vmin, vmax=vmax, cmap='coolwarm')
axes[0].set_title('Pure Oxygen Field')

axes[1].imshow(oxygen_with_tumor, vmin=vmin, vmax=vmax, cmap='coolwarm')
axes[1].set_title('Oxygen (with tumor overlay)')

axes[2].imshow(oxygen_pure - oxygen_with_tumor, cmap='RdBu_r')
axes[2].set_title('Difference')
plt.colorbar(axes[2].images[0])

plt.show()
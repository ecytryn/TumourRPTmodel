import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

plt.ion()
fig, ax = plt.subplots()

# -------------------------------
# PARAMETERS
# -------------------------------
DOMAIN_SIZE = 4000.0  # microns (square domain)
TARGET_DENSITY = 420.0  # vessels/mm² — change this to get desired density

# N is derived from target density and domain area
DOMAIN_AREA_MM2 = (DOMAIN_SIZE / 1000.0) ** 2  # 16 mm²
N = int(round(TARGET_DENSITY * DOMAIN_AREA_MM2))
print(f"Target density: {TARGET_DENSITY} vessels/mm²  →  N = {N} vessels")

R_REPEL = 60.0  # repulsive radius (microns) — controls spacing uniformity
# doesn't need tuning for density; just affects how uniform
# the pattern is. 30 works well for most densities.
STEP_SIZE = 0.3
ITERATIONS = 200
PLOT_EVERY = 20

# Convergence threshold — stop early if max displacement is small
CONVERGE_THRESHOLD = 0.001  # microns

# -------------------------------
# INITIAL RANDOM POSITIONS
# -------------------------------
pts = np.random.rand(N, 2) * DOMAIN_SIZE


def reflect_boundaries(pts, domain_size):
    """Reflect points back into domain (no-flux boundary condition)."""
    # Reflect off left/bottom
    pts = np.abs(pts)
    # Reflect off right/top
    over = pts > domain_size
    pts[over] = 2 * domain_size - pts[over]
    # Clip any residual floating point issues
    pts = np.clip(pts, 0, domain_size)
    return pts


def periodic_boundaries(pts, domain_size):
    """Wrap points back into domain (periodic boundary condition)."""
    pts = pts % domain_size
    return pts


def save_as_grid(pts_um, target_density, xDim=400, yDim=400, cell_size_um=10):
    grid = np.zeros((xDim, yDim), dtype=int)
    xs_um = pts_um[:, 0]
    ys_um = pts_um[:, 1]
    mask = (
        (xs_um >= 0)
        & (xs_um < xDim * cell_size_um)
        & (ys_um >= 0)
        & (ys_um < yDim * cell_size_um)
    )
    xs = (xs_um[mask] / cell_size_um).astype(int)
    ys = (ys_um[mask] / cell_size_um).astype(int)
    grid[xs, ys] = 1
    actual_density = np.sum(grid) / DOMAIN_AREA_MM2
    filename = (
        f"scripts/GenerateVessels/Capillaries_Density{int(round(actual_density))}.csv"
    )
    np.savetxt(filename, grid, delimiter=",", fmt="%d")
    print(f"Saved {filename}")
    print(f"Kept {np.sum(grid)} vessels out of {N} (some may overlap in grid)")
    print(f"Actual vessel density: {actual_density:.1f} vessels/mm²")


# -------------------------------
# RELAXATION LOOP
# -------------------------------
for it in range(ITERATIONS):
    tree = cKDTree(pts, boxsize=DOMAIN_SIZE)
    pairs = tree.query_ball_tree(tree, r=R_REPEL)

    disp = np.zeros_like(pts)

    for i, neigh in enumerate(pairs):
        for j in neigh:
            if i >= j:
                continue
            rij = pts[j] - pts[i]
            # Minimum image convention for periodic boundaries
            rij = rij - DOMAIN_SIZE * np.round(rij / DOMAIN_SIZE)
            dist = np.linalg.norm(rij)
            if dist == 0:
                rij = np.random.randn(2)
                dist = np.linalg.norm(rij)
            direction = rij / dist
            overlap = R_REPEL - dist
            if overlap <= 0:
                continue
            force = overlap * direction
            disp[i] -= force
            disp[j] += force

    pts += STEP_SIZE * disp

    # Periodic boundary
    pts = periodic_boundaries(pts, DOMAIN_SIZE)

    # Check convergence
    max_disp = np.max(np.linalg.norm(STEP_SIZE * disp, axis=1))

    if it % PLOT_EVERY == 0:
        ax.clear()
        ax.scatter(pts[:, 0], pts[:, 1], s=2)
        ax.set_xlim(0, DOMAIN_SIZE)
        ax.set_ylim(0, DOMAIN_SIZE)
        ax.set_aspect("equal")
        ax.set_title(
            f"Iter {it}  max_disp={max_disp:.3f} µm  "
            f"N={N}  target={TARGET_DENSITY} v/mm²"
        )
        plt.pause(0.001)

#    if it > 500 and max_disp < CONVERGE_THRESHOLD:
#        print(f"Converged at iteration {it} (max_disp={max_disp:.4f} µm)")
#        break

# -------------------------------
# SAVE
# -------------------------------
save_as_grid(pts, TARGET_DENSITY)

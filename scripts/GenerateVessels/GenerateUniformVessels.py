import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

plt.ion()

fig, ax = plt.subplots()

# -------------------------------
# PARAMETERS
# -------------------------------
N = 6000
DOMAIN_SIZE = 4000.0      # microns (square domain)
#R_REPEL = 20              # repulsive radius --> 624.5 vessels/mm^2 
#R_REPEL = 40.0            # repulsive radius --> 616.1 vessels/mm^2 
#R_REPEL = 50.0            # repulsive radius --> 604.8 vessels/mm^2
#R_REPEL = 60.0            # repulsive radius --> 591.9 vessels/mm^2 
#R_REPEL = 80.0            # repulsive radius --> 566.3 vessels/mm^2 
R_REPEL = 50.0            # repulsive radius -->  vessels/mm^2 
STEP_SIZE = 0.1           # smaller = slower but more stable relaxation
ITERATIONS = 1000          # number of relaxation iterations
PLOT_EVERY = 20           # show progress

# -------------------------------
# INITIAL RANDOM POSITIONS
# -------------------------------
pts = np.random.rand(N, 2) * DOMAIN_SIZE

def save_as_grid(pts_um, filename,
                 xDim=400, yDim=400, cell_size_um=10):
    """
    Convert vessel coordinates in microns into a 2D grid (integers)
    and save as CSV (same format as original script).

    pts_um : N×2 array of positions in microns.
    filename : CSV output path.

    Each grid cell is 10 × 10 µm.
    """
    # Create empty grid
    grid = np.zeros((xDim, yDim), dtype=int)

    # Extract x,y
    xs_um = pts_um[:, 0]
    ys_um = pts_um[:, 1]

    # Keep only points inside the domain
    mask = (
        (xs_um >= 0) & (xs_um < xDim * cell_size_um) &
        (ys_um >= 0) & (ys_um < yDim * cell_size_um)
    )

    xs_um = xs_um[mask]
    ys_um = ys_um[mask]
        
    # Convert positions to grid index
    xs = (xs_um / cell_size_um).astype(int)
    ys = (ys_um / cell_size_um).astype(int)

    # Mark vessel locations
    grid[xs, ys] = 1

    # Save CSV
    np.savetxt(filename, grid, delimiter=",", fmt="%d")

    print(f"Saved {filename}")
    print(f"Kept {len(xs)} vessels out of {len(pts_um)}")
    print(f"Vessel density: {np.sum(grid)/16:.1f} vessels/mm²")


# -------------------------------
# RELAXATION LOOP
# -------------------------------
for it in range(ITERATIONS):
    tree = cKDTree(pts)
    
    # For every point, get neighbors within R_REPEL
    pairs = tree.query_ball_tree(tree, r=R_REPEL)

    disp = np.zeros_like(pts)

    for i, neigh in enumerate(pairs):
        for j in neigh:
            if i >= j:  
                continue  # avoid double-counting

            rij = pts[j] - pts[i]
            dist = np.linalg.norm(rij)

            if dist == 0:
                continue

            # Repulsive force direction
            direction = rij / dist
            
            # Repulsive magnitude decreases with distance
            overlap = R_REPEL - dist
            force = overlap * direction

            # Equal and opposite forces
            disp[i] -= force
            disp[j] += force

    # Apply displacement with step-size control
    pts += STEP_SIZE * disp

    # Keep points inside domain
#    pts = np.clip(pts, 0, DOMAIN_SIZE)

    if it % PLOT_EVERY == 0:
        ax.clear()
        ax.scatter(pts[:,0], pts[:,1], s=3)
        ax.set_title(f"Iteration {it}")
#        plt.xlim(0, DOMAIN_SIZE)
#        plt.ylim(0, DOMAIN_SIZE)
        plt.gca().set_aspect('equal')
        plt.pause(0.001)

# -------------------------------
# SAVE FINAL CONFIGURATION
# -------------------------------
#np.savetxt("vessels_relaxed.txt", pts, fmt="%.3f")
#print("Saved vessel positions → vessels_relaxed.txt")

save_as_grid(pts, "scripts/GenerateVessels/relaxed_vessels.csv")


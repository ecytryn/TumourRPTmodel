import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# Cell type indices (from code)
# -----------------------------
HEALTHY_CELL_TYPE   = 0
NORMAL_TCELL_TYPE   = 1
HYPO_TCELL_TYPE     = 2
NECRO_TCELL_TYPE    = 3
APOP_TCELL_TYPE     = 4
VESSEL_TYPE         = 5

# -----------------------------
# Which parameter set to visualize
# -----------------------------
INTERVAL  = 20
SKEW      = 10
REPLICATE = 1
SUFFIX    = ""   # if you use one elsewhere

# -----------------------------
# Data loading helpers
# -----------------------------
def load_run_file(interval, skew, replicate, filename, suffix=""):
    sweep_dir = f"results/IntervalSkewSweep{suffix}"
    timestamp = f"2026-02-26_10-56-10"
    run_dir = f"{sweep_dir}/IntervalSkewSweep_{timestamp}/interval_{interval}_skew_{int(skew)}_rep_{replicate}"
    file_path = f"{run_dir}/{filename}"

    if not Path(file_path).exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    data = np.loadtxt(file_path, delimiter=',')
    print(f"Loaded {filename}: shape = {data.shape}")
    return data

def load_populations(interval, skew, replicate, suffix=""):
    return load_run_file(interval, skew, replicate, "populations.csv", suffix)


def load_doses(interval, skew, replicate, suffix=""):
    return load_run_file(interval, skew, replicate, "doses.csv", suffix)

def save_png(interval, skew, replicate):
    # Save figure
    sweep_dir = f"results/IntervalSkewSweep{suffix}"
    timestamp = f"2026-02-26_10-56-10"
    run_dir = f"{sweep_dir}/IntervalSkewSweep_{timestamp}/interval_{interval}_skew_{int(skew)}_rep_{replicate}"
    output_file = f"{run_dir}/pop_dose_i{interval}_s{int(skew)}_r{replicate}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {output_file}")

# -----------------------------
# Load data
# -----------------------------
populations = load_populations(INTERVAL, SKEW, REPLICATE, SUFFIX)
doses       = load_doses(INTERVAL, SKEW, REPLICATE, SUFFIX)

# Time axis (hours → days)
t_hours = np.arange(populations.shape[0])
t_days  = t_hours / 24.0

# -----------------------------
# Tumour population definition
# -----------------------------
# Easy to modify which cell types contribute
tumour_cell_types = [
    NORMAL_TCELL_TYPE,
    HYPO_TCELL_TYPE,
    # NECRO_TCELL_TYPE,  # uncomment if desired
    # APOP_TCELL_TYPE,   # uncomment if desired
]

tumour_population = populations[:, tumour_cell_types].sum(axis=1)

# -----------------------------
# Plot
# -----------------------------
fig, ax1 = plt.subplots(figsize=(8, 5))

# Tumour population curve
ax1.plot(
    t_days,
    tumour_population,
    color="C0",
    lw=2,
    label="Tumour cells (normoxic + hypoxic)"
)

ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Tumour cell count")
ax1.tick_params(axis="y", labelcolor="C0")

# Secondary axis for dose
ax2 = ax1.twinx()

ax2.fill_between(
    t_days,
    doses,
    step="pre",
    alpha=0.3,
    color="C1",
    label="Delivered dose"
)

ax2.set_ylabel("Dose")
ax2.tick_params(axis="y", labelcolor="C1")

# -----------------------------
# Combined legend
# -----------------------------
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()

ax1.legend(
    lines_1 + lines_2,
    labels_1 + labels_2,
    loc="upper right"
)

plt.title(
    f"Interval={INTERVAL}, Skew={SKEW}, Replicate={REPLICATE}"
)

plt.tight_layout()

save_png(INTERVAL, SKEW, REPLICATE)

plt.show()

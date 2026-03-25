#!/usr/bin/env python3
"""
Plot captive RL (N_b^hot + N_ic^hot) vs time from DR sweep simulations
at a fixed receptor density, across all injected amounts.
Overlays analytical QSS approximation to validate the steady-state assumption.

Usage:
    python pk_saturation_plot.py
    python pk_saturation_plot.py 2026-03-19_10-00-00
"""

import glob
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update(
    {
        "figure.figsize": (3.35, 2.4),
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 6,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# =============================================================================
# PK parameters (must match SimParams.java)
# =============================================================================
K_ON = 1.5e-3  # L nmol^-1 min^-1
K_OFF = 1.2e-2  # min^-1
K_INT = 1.0e-3  # min^-1
K_REL = 2.0e-4  # min^-1
LAMBDA_BIO = 1.6e-4  # min^-1
LAMBDA_DEC = 7.1e-5  # min^-1
V_CEN = 0.5  # L
HOT_FRAC = 0.1

# Derived
K_M = V_CEN * (K_OFF + K_INT) / K_ON  # nmol
# K_M in mol for comparison with N_cen which is in mol
K_M_mol = K_M * 1e-9

MIN_PER_DAY = 1440.0

# Target receptor density
TARGET_RECEP = 6e-19  # mol/cell  — adjust if needed
RECEP_TOL = 0.05e-19  # match within ±0.05e-19

INJECTION_DAY = 5


def find_sweep_dir(timestamp=None):
    base = "results/DoseReceptorSweep/DoseReceptorSweep"
    if timestamp:
        d = f"{base}_{timestamp}"
        if not Path(d).exists():
            print(f"ERROR: {d} not found")
            sys.exit(1)
        return d
    dirs = sorted(glob.glob(f"{base}_*"))
    if not dirs:
        print("ERROR: No DoseReceptorSweep directories found.")
        sys.exit(1)
    return dirs[-1]


def load_captive_rl(run_dir):
    """Load N_b_hot + N_ic_hot from a single run directory. Returns array in nmol."""
    pk_file = Path(run_dir) / "pkStateVariables.csv"
    if not pk_file.exists():
        return None
    try:
        df = pd.read_csv(pk_file)
        captive = (df["N_b_hot"] + df["N_ic_hot"]) * 1e9  # mol -> nmol
        return captive.values
    except Exception as e:
        print(f"  Warning: could not read {pk_file}: {e}")
        return None


def analytical_captive(t_days, A_nmol, R_T_nmol):
    """
    QSS analytical approximation for N_captive^H(t).
    t_days: array of time in days (t=0 = injection)
    A_nmol: injected amount in nmol
    R_T_nmol: total receptors in nmol
    """
    t_min = t_days * MIN_PER_DAY
    prefactor = (K_INT / (K_REL + LAMBDA_DEC) + 1.0) * R_T_nmol
    N_cen_hot = HOT_FRAC * A_nmol * np.exp(-(LAMBDA_BIO + LAMBDA_DEC) * t_min)
    N_cen_tot = A_nmol * np.exp(-LAMBDA_BIO * t_min)
    return prefactor * N_cen_hot / (N_cen_tot + K_M)


def main():
    timestamp = sys.argv[1] if len(sys.argv) > 1 else None
    sweep_dir = find_sweep_dir(timestamp)
    print(f"Sweep dir: {sweep_dir}")
    print(f"K_M = {K_M:.1f} nmol")

    # Find all run directories matching target receptor density
    run_dirs = sorted(Path(sweep_dir).glob("dose_*_recep_*_rep_*"))
    if not run_dirs:
        print("ERROR: No run directories found.")
        sys.exit(1)

    # Parse directory names, filter by receptor density
    pattern = re.compile(r"dose_(\d+\.?\d*)_recep_([\d.e+\-]+)_rep_(\d+)")
    runs = []
    for d in run_dirs:
        m = pattern.match(d.name)
        if not m:
            continue
        dose = float(m.group(1))
        recep = float(m.group(2))
        rep = int(m.group(3))
        if abs(recep - TARGET_RECEP) <= RECEP_TOL:
            runs.append((dose, recep, rep, d))

    if not runs:
        print(f"ERROR: No runs found matching recep ≈ {TARGET_RECEP:.2e}")
        sys.exit(1)

    # Group by dose, average replicates
    doses = sorted(set(r[0] for r in runs))
    print(f"Found {len(doses)} dose levels: {doses}")

    # Get R_T from first run's sweep_summary
    summary = pd.read_csv(Path(sweep_dir) / "sweep_summary.csv")
    # We need R_T — infer from cell count and receptor density
    # R_T = n_cells * recep_per_cell (in nmol)
    # Get initial cell count from first run's populations.csv
    first_run_dir = runs[0][3]
    pop_file = first_run_dir / "populations.csv"
    try:
        pops = np.loadtxt(pop_file, delimiter=",", skiprows=1)
        inj_row = INJECTION_DAY * 24
        row = pops[inj_row]
        n_cells = int(row[1] + row[2] + row[3] + row[4])
        R_T_nmol = n_cells * TARGET_RECEP * 1e9  # mol -> nmol
        print(f"Cell count at injection: {n_cells}, R_T = {R_T_nmol:.4e} nmol")
    except Exception as e:
        print(f"Warning: could not read cell count: {e}. Using fallback R_T.")
        R_T_nmol = 28 * TARGET_RECEP * 1e9  # fallback

    # Colour ramp: pale to saturated blue
    blues = [
        "#C6DBEF",
        "#9ECAE1",
        "#6BAED6",
        "#4292C6",
        "#2171B5",
        "#08519C",
        "#08306B",
    ]
    # Interpolate to number of doses
    from matplotlib.colors import LinearSegmentedColormap

    blue_cmap = LinearSegmentedColormap.from_list("blues", blues)
    colours = [blue_cmap(i / (len(doses) - 1)) for i in range(len(doses))]

    fig, ax = plt.subplots()

    for idx, dose in enumerate(doses):
        # Collect all replicates for this dose
        rep_data = []
        for d_val, r_val, rep, run_dir in runs:
            if d_val == dose:
                arr = load_captive_rl(run_dir)
                if arr is not None:
                    rep_data.append(arr)

        if not rep_data:
            continue

        # Align lengths and average
        min_len = min(len(a) for a in rep_data)
        mean_captive = np.mean([a[:min_len] for a in rep_data], axis=0)

        # Time axis: rows are hourly, offset so injection = t=0
        inj_hour = INJECTION_DAY * 24
        t_hours = np.arange(min_len) - inj_hour
        t_days_full = t_hours / 24.0

        # Only plot from injection onwards
        mask = t_days_full >= 0
        t_plot = t_days_full[mask]
        captive_plot = mean_captive[mask]

        colour = colours[idx]
        label = f"{dose:.0f} nmol"

        ax.plot(
            t_plot, captive_plot, color=colour, linewidth=1.2, label=label, zorder=3
        )

        # Analytical overlay
        captive_analytical = analytical_captive(
            t_plot, dose * HOT_FRAC / HOT_FRAC, R_T_nmol
        )
        # Note: A_nmol in analytical is the TOTAL injected (hot+cold combined in N_cen)
        # but only hot fraction matters for N_cen_hot. The formula already uses
        # 0.1*A in numerator and A in denominator so pass full dose.
        captive_analytical = analytical_captive(t_plot, dose, R_T_nmol)
        ax.plot(
            t_plot,
            captive_analytical,
            color=colour,
            linewidth=1.0,
            linestyle="--",
            zorder=2,
        )

    # K_M annotation
    ax.axvline(x=0, color="gray", linewidth=0.6, linestyle=":", zorder=1)

    ax.set_xlabel("Time since injection (days)")
    ax.set_ylabel("$N_\mathrm{captive}^H$ (nmol)")
    # ax.set_yscale("log")

    # Legend: dose levels only (solid lines), plus one entry for solid/dashed meaning
    handles, labels = ax.get_legend_handles_labels()
    # Add proxy artists for solid/dashed explanation
    from matplotlib.lines import Line2D

    proxy_num = Line2D(
        [0], [0], color="gray", linewidth=1.2, linestyle="-", label="Simulation"
    )
    proxy_ana = Line2D(
        [0], [0], color="gray", linewidth=1.0, linestyle="--", label="Analytical"
    )
    ax.legend(
        handles=handles + [proxy_num, proxy_ana], ncol=2, fontsize=6, framealpha=0.9
    )

    ax.grid(True, alpha=0.3, linewidth=0.5)
    plt.tight_layout()

    out_pdf = Path(sweep_dir) / "pk_saturation.pdf"
    out_png = Path(sweep_dir) / "pk_saturation.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved to: {out_pdf}")
    plt.show()


if __name__ == "__main__":
    main()

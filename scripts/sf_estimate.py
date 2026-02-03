#!/usr/bin/env python3
"""
Estimate empirical survival fraction from populations.csv

Flow balance on killable cells (normoxic + hypoxic) each day:

    killable[t] = killable[t-1] + births[t] - new_apoptotic[t] - new_necrotic[t]

Rearranging:
    new_apoptotic[t] = killable[t-1] - killable[t] + births[t] - new_necrotic[t]

We don't observe births directly, but we can recover them:
    births[t] = killable[t] - killable[t-1] + new_apoptotic[t] + new_necrotic[t]

And new_apoptotic can be recovered from the apoptotic column:
    new_apoptotic[t] = apoptotic[t] - apoptotic[t-1] + removed_apoptotic[t]

But we don't observe removed_apoptotic either. However, we can sidestep both
unknowns. The total cell count (all types) satisfies:

    total[t] = total[t-1] + births[t] - removed_apoptotic[t]

(necrotic cells are never removed, apoptotic cells are removed stochastically)

So:
    births[t] - removed_apoptotic[t] = total[t] - total[t-1]

Substituting back into the killable flow balance:

    new_apoptotic[t] = killable[t-1] - killable[t] - new_necrotic[t] + births[t]

And births[t] = (total[t] - total[t-1]) + removed_apoptotic[t]

This still has removed_apoptotic as unknown. Cleanest approach: just use the
apoptotic column directly.

    new_apoptotic[t] = apoptotic[t] - apoptotic[t-1] + removed_apoptotic[t]

Since mean removal time is 2 days, removed_apoptotic[t] ≈ apoptotic[t-2] * removal_rate.
But this is approximate.

ACTUAL CLEANEST: Use total population to back out births, then use killable
balance to get new_apoptotic exactly.

    births[t] = total[t] - total[t-1] + removed_apoptotic[t]

Hmm, still coupled. Let's just do it the direct way:

Given: HEALTHY(0), NORMAL(1), HYPOXIC(2), NECROTIC(3), APOPTOTIC(4), VESSEL(5)

    killable[t]     = NORMAL[t] + HYPOXIC[t]
    new_necrotic[t] = NECROTIC[t] - NECROTIC[t-1]          (exact: necrotic never removed)
    
    total_living[t] = NORMAL[t] + HYPOXIC[t] + NECROTIC[t] + APOPTOTIC[t]
    
    births[t] - removed_apoptotic[t] = total_living[t] - total_living[t-1]
    
    From killable balance:
        new_apoptotic[t] = killable[t-1] - killable[t] - new_necrotic[t] + births[t]
    
    We need births[t] alone. Since removed_apoptotic is stochastic with mean
    removal time 2 days, we estimate:
        removed_apoptotic[t] ≈ apoptotic[t-1] / 2    (approx: half removed each day)
        births[t] ≈ (total_living[t] - total_living[t-1]) + apoptotic[t-1] / 2

    Then:
        new_apoptotic[t] = killable[t-1] - killable[t] - new_necrotic[t] + births[t]
    
    Empirical SF:
        SF[t] = 1 - new_apoptotic[t] / killable[t-1]

    Note: SF here is the average over all killable cells at time t-1. It is
    influenced by the mix of cell ages (birth cohorts) present. Newly born
    cells have SF≈1, old cells may have SF<<1. This gives a population-weighted
    average.

Usage:
    python estimate_sf.py path/to/sim_dir/
    python estimate_sf.py path/to/sim_dir/ --output sf_plot.png
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Match styling from other scripts
mpl.rcParams.update({
    "figure.figsize": (8, 10),
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# populations.csv columns (no header):
# 0: HEALTHY, 1: NORMAL, 2: HYPOXIC, 3: NECROTIC, 4: APOPTOTIC, 5: VESSEL
IDX_NORMAL   = 1
IDX_HYPOXIC  = 2
IDX_NECROTIC = 3
IDX_APOPTOTIC= 4

# pkStateVariables.csv columns:
# NcenHot(0), NcenCold(1), NvHot(2), NvCold(3), NecHot(4), NecCold(5),
# NbHot(6), NbCold(7), NicHot(8), NicCold(9), Ablob(10)
IDX_NB_HOT  = 6
IDX_NIC_HOT = 8


def load_data(sim_dir):
    """Load populations, PK state variables, and SF data if available"""
    pop_file = sim_dir / 'populations.csv'
    pk_file  = sim_dir / 'pkStateVariables.csv'
    sf_file  = sim_dir / 'sfData.csv'

    if not pop_file.exists():
        print(f"ERROR: {pop_file} not found")
        sys.exit(1)

    pop = np.loadtxt(pop_file, delimiter=',')
    pk  = np.loadtxt(pk_file, delimiter=',') if pk_file.exists() else None

    # sfData.csv columns: day, SF_cohort0_norm, SF_cohort0_hypo, SF_prevDay_norm, SF_prevDay_hypo
    sf_model = np.loadtxt(sf_file, delimiter=',') if sf_file.exists() else None
    if sf_model is not None:
        print(f"Loaded sfData.csv ({len(sf_model)} days)")
    else:
        print("No sfData.csv found - will show empirical SF only")

    # Truncate to same length if needed
    if pk is not None:
        n = min(len(pop), len(pk))
        pop = pop[:n]
        pk  = pk[:n]

    return pop, pk, sf_model


def estimate_sf(pop):
    """
    Estimate empirical SF from population flows.

    Returns arrays indexed by day (length = n_days - 1, since we need t and t-1).
    """
    n = len(pop)
    days = np.arange(n)

    killable      = pop[:, IDX_NORMAL] + pop[:, IDX_HYPOXIC]
    necrotic      = pop[:, IDX_NECROTIC]
    apoptotic     = pop[:, IDX_APOPTOTIC]
    total_living  = killable + necrotic + apoptotic

    # Derived quantities (from day 1 onward)
    new_necrotic = np.diff(necrotic)   # exact: necrotic never removed

    # Estimate removed_apoptotic: mean removal time 2 days, so each day
    # roughly half the apoptotic population is removed. Use previous day's count.
    removed_apoptotic = apoptotic[:-1] / 2.0

    # births = change in total_living + removed_apoptotic
    births = np.diff(total_living) + removed_apoptotic

    # new_apoptotic from killable flow balance:
    #   killable[t] = killable[t-1] + births[t] - new_apoptotic[t] - new_necrotic[t]
    new_apoptotic = killable[:-1] - killable[1:] - new_necrotic + births

    # Clamp negative values (can happen due to removal estimate noise)
    new_apoptotic = np.maximum(new_apoptotic, 0.0)

    # Empirical SF
    with np.errstate(divide='ignore', invalid='ignore'):
        sf = 1.0 - new_apoptotic / killable[:-1]
        sf = np.clip(sf, 0.0, 1.0)

    return {
        'days':            days[1:],   # day index for each computed value
        'sf':              sf,
        'killable':        killable,
        'new_apoptotic':   new_apoptotic,
        'new_necrotic':    new_necrotic,
        'births':          births,
    }


def plot_results(pop, pk, sf_data, sf_model=None, output_file=None):
    """Plot SF estimate alongside populations and captive RL"""

    days_full = np.arange(len(pop)) / 24.0   # populations.csv is logged hourly
    days_sf   = sf_data['days'] / 24.0        # same hourly indices

    killable      = sf_data['killable']
    new_apoptotic = sf_data['new_apoptotic']
    sf            = sf_data['sf']

    # Panel count: base 4 (populations, captive RL, new apoptotic, SF)
    # Plus 2 more if we have sfData (D and A/G_num)
    n_panels = 4 if pk is not None else 3
    if sf_model is not None:
        n_panels += 2
    fig, axes = plt.subplots(n_panels, 1, figsize=(8, 2.5*n_panels), sharex=True)

    # --- Panel 1: Populations ---
    ax = axes[0]
    ax.plot(days_full, pop[:, IDX_NORMAL],   label='Normoxic', color='blue')
    ax.plot(days_full, pop[:, IDX_HYPOXIC],  label='Hypoxic',  color='orange')
    ax.plot(days_full, pop[:, IDX_NECROTIC], label='Necrotic', color='gray', linestyle='--')
    ax.plot(days_full, pop[:, IDX_APOPTOTIC],label='Apoptotic',color='red',  linestyle='--')
    ax.set_ylabel('Cell count')
    ax.set_title('Cell populations')
    ax.legend(loc='upper left', fontsize=7)
    ax.set_yscale('log')
    ax.set_ylim(bottom=0.5)

    # --- Panel 2: Captive RL (if PK data available) ---
    if pk is not None:
        ax = axes[1]
        captive_rl = (pk[:, IDX_NB_HOT] + pk[:, IDX_NIC_HOT]) * 1e9  # convert to nmol
        ax.plot(days_full, captive_rl, color='darkred')
        ax.set_ylabel('Captive RL (nmol)')
        ax.set_title('Bound + Internalized radioligand (hot)')
        ax.set_yscale('log')
        panel_offset = 2
    else:
        panel_offset = 1

    # --- Panel 3: New apoptotic cells per day ---
    ax = axes[panel_offset]
    ax.plot(days_sf, new_apoptotic, color='red')
    ax.set_ylabel('New apoptotic / day')
    ax.set_title('Radiation-induced apoptosis rate')
#    ax.set_ylim(bottom=0)
    ax.set_yscale('log')

    # --- Panel 4: SF - empirical estimate and model values on same axes ---
    ax = axes[panel_offset + 1]
    ax.plot(days_sf, sf, color='darkgreen', label='Empirical (population flow)')

    if sf_model is not None:
        # sfData columns:
        #   day(0), SF_c0_norm(1), SF_c0_hypo(2), SF_prev_norm(3), SF_prev_hypo(4),
        #   D_c0(5), A_c0(6), Gnum_c0(7), D_prev(8), A_prev(9), Gnum_prev(10)
        sf_days      = sf_model[:, 0]
        sf_c0_norm   = sf_model[:, 1]
        sf_c0_hypo   = sf_model[:, 2]
        sf_prev_norm = sf_model[:, 3]
        sf_prev_hypo = sf_model[:, 4]

        ax.plot(sf_days, sf_c0_norm,   color='blue',   linestyle='--', label='Day-0 cohort (normoxic)')
        ax.plot(sf_days, sf_c0_hypo,   color='orange', linestyle='--', label='Day-0 cohort (hypoxic)')
        ax.plot(sf_days, sf_prev_norm, color='blue',   linestyle=':',  label='Prev-day cohort (normoxic)')
        ax.plot(sf_days, sf_prev_hypo, color='orange', linestyle=':',  label='Prev-day cohort (hypoxic)')

    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8)
    ax.set_ylabel('Survival fraction')
    ax.set_xlabel('Day')
    ax.set_title('SF: model cohorts vs empirical estimate')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='lower left', fontsize=7)

    # --- Panels 5 & 6: D, A, G_num (only if sfData available) ---
    if sf_model is not None:
        # Panel 5: Accumulated dose D
        ax = axes[panel_offset + 2]
        ax.plot(sf_days, sf_model[:, 5], color='blue',   linestyle='--', label='Day-0 cohort')
        ax.plot(sf_days, sf_model[:, 8], '.', color='blue',   linestyle='',  label='Prev-day cohort')
        ax.set_ylabel('D (Gy)')
        ax.set_title('Accumulated dose per cohort')
        ax.legend(loc='upper left', fontsize=7)
        ax.set_yscale('log')

        # Panel 6: A and G_num
        ax = axes[panel_offset + 3]
        ax.plot(sf_days, sf_model[:, 9], 'o', color='green',    linestyle='',  label='A (prev-day)')
        ax.plot(sf_days, sf_model[:, 6], '.', color='red',    linestyle='', label='A (day-0)')
        ax.plot(sf_days, sf_model[:, 7], color='purple', linestyle='--', label='G_num (day-0)')
        ax.plot(sf_days, sf_model[:,10], '.', color='purple', linestyle='',  label='G_num (prev-day)')
        ax.set_ylabel('ODE state')
        ax.set_xlabel('Day')
        ax.set_title('Damage accumulation: A and G_num')
        ax.legend(loc='upper left', fontsize=7)
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-20)

    fig.tight_layout()

    if output_file:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    else:
        plt.show()

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Estimate empirical SF from populations.csv',
        epilog="""
Examples:
  python estimate_sf.py results/debug_runs/Debug_I36_S0_20260130/
  python estimate_sf.py results/debug_runs/Debug_I36_S0_20260130/ --output sf.png
        """
    )
    parser.add_argument('sim_dir', help='Path to simulation output directory')
    parser.add_argument('--output', '-o', help='Output image file (default: show interactively)')

    args = parser.parse_args()
    sim_dir = Path(args.sim_dir)

    if not sim_dir.is_dir():
        print(f"ERROR: {sim_dir} is not a directory")
        sys.exit(1)

    pop, pk, sf_model = load_data(sim_dir)
    sf_data = estimate_sf(pop)

    # Default output: save to current working directory
    output_to_source_forder  = sim_dir / 'sf_estimate.png'
    
    output_file = args.output if args.output else output_to_source_forder
    plot_results(pop, pk, sf_data, sf_model=sf_model, output_file=output_file)


if __name__ == '__main__':
    main()

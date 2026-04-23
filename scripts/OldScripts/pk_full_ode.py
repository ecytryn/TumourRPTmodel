#!/usr/bin/env python3
"""
Full PK ODE model - integrates all 10 ODEs without QSS approximation

IMPORTANT: This is a STIFF ODE system due to multiple timescales:
- Fast: Vascular/EC equilibration (~milliseconds to seconds)
- Medium: Binding/unbinding (~tens of seconds)  
- Slow: Internalization, release, clearance, decay (~hours to days)

Uses scipy's solve_ivp with 'Radau' method (implicit solver for stiff systems).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# PARAMETERS - MUST MATCH SimParams.java EXACTLY
# =============================================================================

# Rate constants (SI units: 1/s or m³/(mol·s))
lambda_bio = 1.6e-4 / 60.0      # s^-1 (biological clearance)
lambda_decay = 7.14e-5 / 60.0   # s^-1 (Lu-177 decay)
k_on = 0.046 * 1e6 / 60.0       # m³/(mol·s) (binding rate)
k_off = 0.368 / 60.0            # s^-1 (unbinding rate)
k_int = 0.001 / 60.0            # s^-1 (internalization)
k_rel = 2e-4 / 60.0             # s^-1 (release from cells)

# Flow and permeability (from SimParams.java - per vessel, scaled up)
# NOTE: These are for a SINGLE vessel in SimParams, we need to scale by number of vessels
VESSEL_FLOW = 5e-9 / 1000.0     # m³/s per vessel (5 nL/s)
VESSEL_PS = (6.0/5.0) * VESSEL_FLOW  # m³/s per vessel

# =============================================================================
# TUMOR GEOMETRY - INPUT FROM JAVA SIMULATION
# =============================================================================

# **MANUALLY SET THESE FROM JAVA CONSOLE OUTPUT**
tumour_cells_2D = 3285          # From "Tumor cells (2D): " output
numberOfVessels = 365           # Typical value for small tumor

# Geometric constants (from SimParams.java)
CELL_LENGTH = 1e-5              # m (10 μm)
INTERSTITIAL_FRACTION = 0.4     # Fraction of tumor that is extracellular
AVOGADRO = 6.022e23             # molecules/mol
RECEPTORS_PER_CELL = 4e5        # receptors/cell

# =============================================================================
# COMPUTE VOLUMES (matching SimParams.computeTumorVolume)
# =============================================================================

# Tumor geometry (cylindrical extrusion)
radius_cells = np.sqrt(tumour_cells_2D / np.pi)  # radius in cell lengths
radius_m = radius_cells * CELL_LENGTH             # m
height_m = 2.0 * radius_m                         # m (cylinder height)

# Tumor volume
tumourVolume_m3 = np.pi * radius_m**2 * height_m  # m³

# Extracellular volume
V_ec = INTERSTITIAL_FRACTION * tumourVolume_m3    # m³

# Vascular volume
vesselVolume = (CELL_LENGTH**2) * height_m        # m³ per vessel
V_v = numberOfVessels * vesselVolume              # m³

# Central compartment (fixed)
V_cen = 0.458e-3                                  # m³ (0.458 L)

# Scale flow rates by number of vessels near tumor
F_T = numberOfVessels * VESSEL_FLOW               # m³/s (total flow)
PS = numberOfVessels * VESSEL_PS                  # m³/s (total permeability)

# =============================================================================
# RECEPTOR CALCULATIONS
# =============================================================================

# Total cells in 3D (with cylindrical extrusion)
height_cells = height_m / CELL_LENGTH
totalCells3D = tumour_cells_2D * height_cells

# Total receptors (in moles)
R_total = totalCells3D * RECEPTORS_PER_CELL / AVOGADRO  # mol

# =============================================================================
# INITIAL CONDITIONS
# =============================================================================

dose_mol = 100e-9               # mol (100 nmol injection)
hot_fraction = 0.1              # 10% radioactive

# Initial state vector [10 variables, all amounts in mol]
# Order: N_cen_H, N_cen_C, N_v_H, N_v_C, N_ec_H, N_ec_C, N_b_H, N_b_C, N_ic_H, N_ic_C
y0 = np.zeros(10)
y0[0] = hot_fraction * dose_mol      # N_cen_hot (injected)
y0[1] = (1 - hot_fraction) * dose_mol # N_cen_cold (injected)
# All other compartments start at zero

# =============================================================================
# ODE SYSTEM - FULL 10 EQUATIONS
# =============================================================================

def pk_odes(t, y):
    """
    Full PK ODE system without QSS approximation
    
    State vector y:
    [0] N_cen_H   - Central hot (mol)
    [1] N_cen_C   - Central cold (mol)
    [2] N_v_H     - Vascular hot (mol)
    [3] N_v_C     - Vascular cold (mol)
    [4] N_ec_H    - Extracellular hot (mol)
    [5] N_ec_C    - Extracellular cold (mol)
    [6] N_b_H     - Bound hot (mol)
    [7] N_b_C     - Bound cold (mol)
    [8] N_ic_H    - Intracellular hot (mol)
    [9] N_ic_C    - Intracellular cold (mol)
    
    Note: solve_ivp requires (t, y) order, not (y, t) like odeint
    """
    
    # Unpack state
    N_cen_H, N_cen_C = y[0], y[1]
    N_v_H, N_v_C = y[2], y[3]
    N_ec_H, N_ec_C = y[4], y[5]
    N_b_H, N_b_C = y[6], y[7]
    N_ic_H, N_ic_C = y[8], y[9]
    
    # Concentrations (needed for rate laws)
    C_cen_H = N_cen_H / V_cen
    C_cen_C = N_cen_C / V_cen
    C_v_H = N_v_H / V_v
    C_v_C = N_v_C / V_v
    C_ec_H = N_ec_H / V_ec
    C_ec_C = N_ec_C / V_ec
    
    # Free receptors (total - occupied)
    N_R_free = R_total - (N_b_H + N_b_C)
    C_R_free = N_R_free / V_ec
    
    # Ensure non-negative concentrations (numerical safety)
    C_R_free = max(0, C_R_free)
    
    # Initialize derivatives
    dy = np.zeros(10)
    
    # =========================================================================
    # HOT COMPARTMENTS (radioactive)
    # =========================================================================
    
    # Central hot
    dy[0] = (
        - lambda_bio * N_cen_H                  # Clearance
        - F_T * (C_cen_H - C_v_H)               # Flow: mol/s = (m³/s) × (mol/m³ difference)
        - lambda_decay * N_cen_H                # Decay
    )
    
    # Vascular hot
    dy[2] = (
        + F_T * (C_cen_H - C_v_H)               # Flow from central
        - PS * (C_v_H - C_ec_H)                 # Permeability OUT (corrected sign)
        - lambda_decay * N_v_H                  # Decay
    )
    
    # Extracellular hot
    dy[4] = (
        + PS * (C_v_H - C_ec_H)                 # Permeability IN (corrected sign)
        - k_on * C_ec_H * C_R_free * V_ec       # Binding (rate × volume)
        + k_off * N_b_H                         # Unbinding
        - lambda_decay * N_ec_H                 # Decay
    )
    
    # Bound hot
    dy[6] = (
        + k_on * C_ec_H * C_R_free * V_ec       # Binding
        - k_off * N_b_H                         # Unbinding
        - k_int * N_b_H                         # Internalization
        - lambda_decay * N_b_H                  # Decay
    )
    
    # Intracellular hot
    dy[8] = (
        + k_int * N_b_H                         # Internalization
        - k_rel * N_ic_H                        # Release
        - lambda_decay * N_ic_H                 # Decay
    )
    
    # =========================================================================
    # COLD COMPARTMENTS (decayed)
    # =========================================================================
    
    # Central cold
    dy[1] = (
        - lambda_bio * N_cen_C                  # Clearance
        - F_T * (C_cen_C - C_v_C)               # Flow to/from vascular
        + lambda_decay * N_cen_H                # Decay FROM hot
    )
    
    # Vascular cold
    dy[3] = (
        + F_T * (C_cen_C - C_v_C)               # Flow from central
        - PS * (C_v_C - C_ec_C)                 # Permeability OUT
        + lambda_decay * N_v_H                  # Decay FROM hot
    )
    
    # Extracellular cold
    dy[5] = (
        + PS * (C_v_C - C_ec_C)                 # Permeability IN
        - k_on * C_ec_C * C_R_free * V_ec       # Binding
        + k_off * N_b_C                         # Unbinding
        + lambda_decay * N_ec_H                 # Decay FROM hot
    )
    
    # Bound cold
    dy[7] = (
        + k_on * C_ec_C * C_R_free * V_ec       # Binding
        - k_off * N_b_C                         # Unbinding
        - k_int * N_b_C                         # Internalization
        + lambda_decay * N_b_H                  # Decay FROM hot
    )
    
    # Intracellular cold
    dy[9] = (
        + k_int * N_b_C                         # Internalization
        - k_rel * N_ic_C                        # Release
        + lambda_decay * N_ic_H                 # Decay FROM hot
    )
    
    return dy

# =============================================================================
# TIME INTEGRATION WITH STIFF SOLVER
# =============================================================================

t_max_days = 30                              # 30 days post-injection
t_max_sec = t_max_days * 24 * 3600           # Convert to seconds

# Time span for integration
t_span = (0, t_max_sec)

# Evaluation points (hourly)
t_eval = np.arange(0, t_max_sec, 3600)       # Every hour

print("\n" + "="*70)
print("FULL ODE MODEL - INTEGRATING WITH STIFF SOLVER")
print("="*70)
print(f"\nIntegrating over {t_max_days} days...")
print(f"Using Radau method (implicit, good for stiff systems)")
print(f"Evaluation points: {len(t_eval)} (hourly)")

# Integrate with stiff solver
# Radau: Implicit Runge-Kutta method of order 5 (good for stiff ODEs)
sol = solve_ivp(
    pk_odes, 
    t_span, 
    y0, 
    t_eval=t_eval,
    method='Radau',       # Stiff solver
    rtol=1e-8,            # Relative tolerance
    atol=1e-10            # Absolute tolerance (mol scale)
)

if not sol.success:
    print(f"WARNING: Integration failed - {sol.message}")
else:
    print("Integration successful!")

# Extract time and solution
t_sec = sol.t
t_hours = t_sec / 3600
t_days = t_sec / 86400
solution = sol.y.T  # Transpose to match (time, variables) shape

# =============================================================================
# EXTRACT RESULTS
# =============================================================================

# Extract time series
N_cen_H = solution[:, 0]
N_cen_C = solution[:, 1]
N_v_H = solution[:, 2]
N_v_C = solution[:, 3]
N_ec_H = solution[:, 4]
N_ec_C = solution[:, 5]
N_b_H = solution[:, 6]
N_b_C = solution[:, 7]
N_ic_H = solution[:, 8]
N_ic_C = solution[:, 9]

# Totals
N_cen_total = N_cen_H + N_cen_C
N_ic_total = N_ic_H + N_ic_C
N_tumor_hot = N_v_H + N_ec_H + N_b_H + N_ic_H

# Save results to dictionary
results = {
    't_hours': t_hours,
    't_days': t_days,
    't_sec': t_sec,
    'N_cen_hot': N_cen_H,          # For comparison
    'N_ic_hot': N_ic_H,            # For comparison
    'N_cen_total': N_cen_total,
    'N_ic_total': N_ic_total,
    'N_v_hot': N_v_H,
    'N_ec_hot': N_ec_H,
    'N_b_hot': N_b_H,
    'N_tumor_hot': N_tumor_hot,
}

# =============================================================================
# PRINT SUMMARY
# =============================================================================

print("\n=== FULL ODE MODEL SUMMARY ===")
print(f"\nParameters:")
print(f"  λ_bio = {lambda_bio:.4e} s⁻¹")
print(f"  λ_decay = {lambda_decay:.4e} s⁻¹")
print(f"  k_on = {k_on:.4e} m³/(mol·s)")
print(f"  k_off = {k_off:.4e} s⁻¹")
print(f"  k_int = {k_int:.4e} s⁻¹")
print(f"  k_rel = {k_rel:.4e} s⁻¹")
print(f"  F_T = {F_T:.4e} m³/s ({numberOfVessels} vessels)")
print(f"  PS = {PS:.4e} m³/s ({numberOfVessels} vessels)")

print(f"\nGeometry:")
print(f"  Tumor cells (2D): {tumour_cells_2D}")
print(f"  Tumor cells (3D): {totalCells3D:.0f}")
print(f"  V_ec = {V_ec:.6e} m³")
print(f"  V_v = {V_v:.6e} m³")
print(f"  V_cen = {V_cen:.6e} m³")
print(f"  R_total = {R_total:.4e} mol")

print(f"\nPeak values:")
print(f"  N_cen_hot (max) = {N_cen_H.max():.4e} mol at t = {t_days[N_cen_H.argmax()]:.1f} days")
print(f"  N_ic_hot (max) = {N_ic_H.max():.4e} mol at t = {t_days[N_ic_H.argmax()]:.1f} days")
print(f"  N_b_hot (max) = {N_b_H.max():.4e} mol at t = {t_days[N_b_H.argmax()]:.1f} days")

print(f"\nAt t = 1 day:")
idx_1day = np.argmin(np.abs(t_days - 1.0))
print(f"  N_cen_hot = {N_cen_H[idx_1day]:.4e} mol")
print(f"  N_v_hot = {N_v_H[idx_1day]:.4e} mol")
print(f"  N_ec_hot = {N_ec_H[idx_1day]:.4e} mol")
print(f"  N_b_hot = {N_b_H[idx_1day]:.4e} mol")
print(f"  N_ic_hot = {N_ic_H[idx_1day]:.4e} mol")

print(f"\nAt t = 7 days:")
idx_7day = np.argmin(np.abs(t_days - 7.0))
print(f"  N_cen_hot = {N_cen_H[idx_7day]:.4e} mol")
print(f"  N_v_hot = {N_v_H[idx_7day]:.4e} mol")
print(f"  N_ec_hot = {N_ec_H[idx_7day]:.4e} mol")
print(f"  N_b_hot = {N_b_H[idx_7day]:.4e} mol")
print(f"  N_ic_hot = {N_ic_H[idx_7day]:.4e} mol")

print("\n=== COMPLETE ===\n")

# Note: Plotting handled by pk_comparison_plot.py

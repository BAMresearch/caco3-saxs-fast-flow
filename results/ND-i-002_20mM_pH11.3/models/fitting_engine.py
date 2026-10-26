#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Fitting Engine and Optimization Workflows for ND-i-009.
Handles:
- Parameter object initialization tailored for Merged, Short Distance, and Long Distance
- Residuals and weighted objective functions
- Low-q oscillation detection for conditional large spheres
- Signal quality screening for difficult initial frames
- Multiprocessing chunking with all CPUs - 1
"""

import os
import multiprocessing
import numpy as np
import pandas as pd
from lmfit import Parameters, minimize, fit_report
from scipy.signal import savgol_filter

from .scattering_models import (
    calc_model_intensity,
    RHO2_CALIBRATED
)


def get_available_cpus():
    """
    Returns max_workers = os.cpu_count() - 1, keeping 1 CPU free for responsiveness.
    """
    total = os.cpu_count() or 4
    return max(1, total - 1)


def create_merged_params():
    """
    Creates Parameters for the merged dataset (0.002 to 2.5 nm^-1).
    All components (primary fractal, large spheres, OZ fluctuations) active.
    """
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)
    p.add('phi', value=0.006, min=1e-5, max=0.05, vary=True)
    p.add('Rmean', value=8.5, min=3.0, max=20.0, vary=True)
    p.add('stdDev', value=3.3, min=0.5, max=10.0, vary=True)
    p.add('D', value=1.92, min=1.1, max=2.999, vary=True)
    p.add('xi', value=3000.0, min=100.0, max=30000.0, vary=True)
    p.add('bkg', value=0.0, min=0.0, max=5.0, vary=False)

    # Large sphere population (active in merged)
    p.add('mean_radius_large', value=3200.0, min=500.0, max=8000.0, vary=True)
    p.add('std_radius_large', value=150.0, min=10.0, max=500.0, vary=True)
    p.add('total_volume_fraction_large', value=0.0007, min=1e-7, max=0.01, vary=True)

    # OZ fluctuations (active in merged)
    p.add('A_OZ', value=0.069, min=0.0, max=10.0, vary=True)
    p.add('xi_OZ', value=6.3, min=0.5, max=50.0, vary=True)

    return p


def create_short_distance_params(baseline_params=None):
    """
    Creates Parameters for short distance (3 m geometry, q: 0.022 - 2.5 nm^-1).
    - Large spheres EXCLUDED (outside window, phi_large=0, vary=False)
    - Correlation length xi FIXED to baseline value (outside window)
    - Primary particles (Rmean, stdDev, phi), fractal D, and OZ fluctuations ACTIVE
    """
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)

    phi_val = 0.006 if baseline_params is None else baseline_params.get('phi', 0.006)
    Rmean_val = 8.5 if baseline_params is None else baseline_params.get('Rmean', 8.5)
    stdDev_val = 3.3 if baseline_params is None else baseline_params.get('stdDev', 3.3)
    D_val = 1.92 if baseline_params is None else baseline_params.get('D', 1.92)
    xi_val = 3000.0 if baseline_params is None else baseline_params.get('xi', 3000.0)
    A_OZ_val = 0.069 if baseline_params is None else baseline_params.get('A_OZ', 0.069)
    xi_OZ_val = 6.3 if baseline_params is None else baseline_params.get('xi_OZ', 6.3)
    bkg_val = 0.0 if baseline_params is None else baseline_params.get('bkg', 0.0)

    p.add('phi', value=phi_val, min=1e-6, max=0.05, vary=True)
    p.add('Rmean', value=Rmean_val, min=2.0, max=25.0, vary=True)
    p.add('stdDev', value=stdDev_val, min=0.2, max=12.0, vary=True)
    p.add('D', value=D_val, min=1.1, max=2.999, vary=True)
    p.add('xi', value=xi_val, vary=False)
    p.add('bkg', value=bkg_val, min=0.0, max=5.0, vary=False)

    # Large spheres NOT applicable for short distance
    p.add('mean_radius_large', value=3200.0, vary=False)
    p.add('std_radius_large', value=150.0, vary=False)
    p.add('total_volume_fraction_large', value=0.0, vary=False)

    # High-q OZ fluctuations are active
    p.add('A_OZ', value=A_OZ_val, min=0.0, max=10.0, vary=True)
    p.add('xi_OZ', value=xi_OZ_val, min=0.2, max=50.0, vary=True)

    return p




def create_short_distance_viscoelastic_params(baseline_params=None):
    """
    Creates Parameters for Short Distance with Viscoelastic Network Fluctuations (Fisher-Burford).
    - Large spheres EXCLUDED (phi_large=0)
    - xi FIXED to merged baseline value
    - Primary particles (phi, Rmean, stdDev, D) ACTIVE
    - Viscoelastic Fluctuations (A_fluc, xi_fluc, mu_fluc) ACTIVE
    - Flat background ACTIVE
    """
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)

    phi_val = 0.003 if baseline_params is None else baseline_params.get('phi', 0.003)
    Rmean_val = 8.5 if baseline_params is None else baseline_params.get('Rmean', 8.5)
    stdDev_val = 3.0 if baseline_params is None else baseline_params.get('stdDev', 3.0)
    D_val = 1.85 if baseline_params is None else baseline_params.get('D', 1.85)
    xi_val = 3000.0 if baseline_params is None else baseline_params.get('xi', 3000.0)
    bkg_val = 0.015 if baseline_params is None else baseline_params.get('bkg', 0.015)

    p.add('phi', value=phi_val, min=1e-6, max=0.05, vary=True)
    p.add('Rmean', value=Rmean_val, min=2.0, max=25.0, vary=True)
    p.add('stdDev', value=stdDev_val, min=0.2, max=12.0, vary=True)
    p.add('D', value=D_val, min=1.1, max=2.999, vary=True)
    p.add('xi', value=xi_val, vary=False)
    p.add('bkg', value=bkg_val, min=0.0, max=0.2, vary=True)

    # Large spheres parked at 0
    p.add('mean_radius_large', value=3200.0, vary=False)
    p.add('std_radius_large', value=150.0, vary=False)
    p.add('total_volume_fraction_large', value=0.0, vary=False)

    # Viscoelastic Fluctuation Parameters (Fisher-Burford)
    A_fluc_val = 0.55 if baseline_params is None else baseline_params.get('A_fluc', 0.55)
    xi_fluc_val = 4.6 if baseline_params is None else baseline_params.get('xi_fluc', 4.6)
    mu_fluc_val = 4.0 if baseline_params is None else baseline_params.get('mu_fluc', 4.0)

    p.add('A_fluc', value=A_fluc_val, min=0.0, max=50.0, vary=True)
    p.add('xi_fluc', value=xi_fluc_val, min=0.5, max=25.0, vary=True)
    p.add('mu_fluc', value=mu_fluc_val, vary=False)  # Locked to Porod boundary decay mu=4.0

    return p


def create_long_distance_params(baseline_params=None, include_large_spheres=False):
    """
    Creates Parameters for long distance (31 m geometry, q: 0.002 - 0.22 nm^-1).
    - Primary particle dimensions Rmean and stdDev are FIXED (knee outside window)
    - OZ fluctuations EXCLUDED (A_OZ=0, xi_OZ=0, vary=False)
    - Aggregate scale phi, fractal D, and correlation cutoff xi are ACTIVE
    - Large spheres CONDITIONAL (active only if include_large_spheres=True)
    """
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)

    phi_val = 0.006 if baseline_params is None else baseline_params.get('phi', 0.006)
    Rmean_val = 8.5 if baseline_params is None else baseline_params.get('Rmean', 8.5)
    stdDev_val = 3.3 if baseline_params is None else baseline_params.get('stdDev', 3.3)
    D_val = 1.92 if baseline_params is None else baseline_params.get('D', 1.92)
    xi_val = 3000.0 if baseline_params is None else baseline_params.get('xi', 3000.0)
    bkg_val = 0.0 if baseline_params is None else baseline_params.get('bkg', 0.0)

    # Rmean and stdDev transferred from short distance or merged baseline
    p.add('Rmean', value=Rmean_val, vary=False)
    p.add('stdDev', value=stdDev_val, vary=False)

    # Aggregate network parameters are refined
    p.add('phi', value=phi_val, min=1e-6, max=0.08, vary=True)
    p.add('D', value=D_val, min=1.1, max=2.999, vary=True)
    p.add('xi', value=xi_val, min=50.0, max=50000.0, vary=True)
    p.add('bkg', value=bkg_val, min=0.0, max=5.0, vary=False)

    # High-q OZ fluctuations are outside long distance window
    p.add('A_OZ', value=0.0, vary=False)
    p.add('xi_OZ', value=0.0, vary=False)

    # Large spheres conditional configuration
    if include_large_spheres:
        R_large_val = 3200.0 if baseline_params is None else baseline_params.get('mean_radius_large', 3200.0)
        std_large_val = 150.0 if baseline_params is None else baseline_params.get('std_radius_large', 150.0)
        phi_large_val = 0.0007 if baseline_params is None else baseline_params.get('total_volume_fraction_large', 0.0007)

        p.add('mean_radius_large', value=R_large_val, min=1000.0, max=8000.0, vary=True)
        p.add('std_radius_large', value=std_large_val, min=10.0, max=500.0, vary=True)
        p.add('total_volume_fraction_large', value=phi_large_val, min=1e-8, max=0.005, vary=True)
    else:
        p.add('mean_radius_large', value=3200.0, vary=False)
        p.add('std_radius_large', value=150.0, vary=False)
        p.add('total_volume_fraction_large', value=0.0, vary=False)

    return p


def residuals_func(params, x, y_meas, dy_meas, include_large_spheres=True, include_oz=True, low_q_weight_factor=1.0, integration_steps=100):
    """
    Weighted residual calculation: (y_meas - y_model) / dy_meas.
    Applies optional low-q weighting factor when fine-tuning sphere oscillations.
    """
    y_calc = calc_model_intensity(params, x, include_large_spheres=include_large_spheres, include_oz=include_oz, integration_steps=integration_steps)
    weights = 1.0 / np.maximum(dy_meas, 1e-4)

    if low_q_weight_factor > 1.0:
        low_q_mask = (x < 0.006)
        weights[low_q_mask] *= low_q_weight_factor

    diff = (y_meas - y_calc) * weights
    # Guard against NaN/Inf
    return np.nan_to_num(diff, nan=1e4, posinf=1e4, neginf=-1e4)


def check_frame_signal_quality(q, y_meas, dy_meas, q_cutoff=0.05):
    """
    Assesses signal-to-noise ratio and data integrity for a single frame.
    In long distance SAXS, signal is concentrated in q < 0.05 nm^-1.
    Detects induction period frames dominated by noise, low amplitude, or negative intensities.
    """
    mask = (q <= q_cutoff)
    q_sub = q[mask]
    y_sub = y_meas[mask]
    dy_sub = dy_meas[mask]

    # Check for negative values in the innermost low-q region (q < 0.003 nm^-1)
    innermost_mask = (q_sub < 0.003)
    if np.any(y_sub[innermost_mask] <= 0):
        return False, 0.0, "Negative intensity observed in innermost low-q points (q < 0.003 nm^-1)"

    positive_mask = (y_sub > 0)
    if np.sum(positive_mask) < 20:
        return False, 0.0, "Insufficient positive points in low-q"

    # In reaction frames, maximum intensity in low-q exceeds 1e5 cm^-1
    # Induction/pre-nucleation frames have max intensity an order of magnitude lower (< 5e4)
    max_intensity = np.max(y_sub)
    if max_intensity < 5e4:
        return False, float(np.median(y_sub / np.maximum(dy_sub, 1e-4))), f"Pre-nucleation low intensity (max I = {max_intensity:.2e} < 5e4)"

    snr = np.median(y_sub[positive_mask] / np.maximum(dy_sub[positive_mask], 1e-4))
    if snr < 2.5:
        return False, float(snr), f"Low low-q SNR ({snr:.2f})"

    return True, float(snr), "Good signal"


def detect_low_q_oscillations(q, y_meas, dy_meas, q_cutoff=0.006):
    """
    Detects whether low-q region (q < 0.006 nm^-1) exhibits clear oscillations
    from micron-sized sphere populations (R ~ 3.2 um) vs monotonic power-law.
    Spherical form factor of R ~ 3.2 um produces a characteristic Kratky (I*q^2)
    peak at q ~ 0.0031 nm^-1 and minimum at q ~ 0.0038 nm^-1.
    """
    mask = (q >= 0.002) & (q <= q_cutoff) & (y_meas > 0)
    if np.sum(mask) < 15:
        return False, 0.0, {"reason": "Insufficient valid points in low-q"}

    q_sub = q[mask]
    y_sub = y_meas[mask]

    # Evaluate Kratky transform I * q^2 to remove the dominant 1/q^2 fractal background
    kratky = y_sub * (q_sub**2)
    
    # Locate characteristic spherical peak and trough
    idx_min1 = np.argmin(np.abs(q_sub - 0.0026))
    idx_peak = np.argmin(np.abs(q_sub - 0.0031))
    idx_trough = np.argmin(np.abs(q_sub - 0.0038))

    v_min1 = kratky[idx_min1]
    v_peak = kratky[idx_peak]
    v_trough = kratky[idx_trough]

    # Clear oscillation requires peak significantly exceeding subsequent trough and preceding knee
    ratio = v_peak / max(1e-6, v_trough)
    has_osc = (ratio >= 1.35) and (v_peak >= 0.95 * v_min1)
    
    details = {
        "v_min1": float(v_min1),
        "v_peak": float(v_peak),
        "v_trough": float(v_trough),
        "peak_to_trough_ratio": float(ratio)
    }

    return has_osc, float(ratio), details


def fit_single_frame(q, y, dy, init_params, q_range=None, include_large_spheres=True, 
                     include_oz=True, low_q_weight_factor=1.0, method='least_squares'):
    """
    Optimizes parameters for a single scattering curve using lmfit.minimize.
    """
    if q_range is not None:
        sel = (q >= q_range[0]) & (q <= q_range[1]) & (y > -1e3)
    else:
        sel = (y > -1e3)

    q_fit = q[sel]
    y_fit = y[sel]
    dy_fit = dy[sel]

    # Pre-clean non-positive points if near 0
    clean_mask = (dy_fit > 0)
    q_fit = q_fit[clean_mask]
    y_fit = y_fit[clean_mask]
    dy_fit = dy_fit[clean_mask]

    out = minimize(
        residuals_func,
        init_params,
        args=(q_fit, y_fit, dy_fit, include_large_spheres, include_oz, low_q_weight_factor),
        method=method,
        max_nfev=600
    )

    return out

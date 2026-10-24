#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Fit Merged Dataset for ND-i-009 as Physical Starting Point.

Establishes global baseline parameters across q in [0.002, 2.5] nm^-1:
- Primary particles: Rmean, stdDev, phi
- Mass fractal network: D, xi
- Large sphere population: R_large, std_large, phi_large
- High-q OZ concentration fluctuations: A_OZ, xi_OZ
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Ensure models module can be imported
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from models.scattering_models import (
    calc_model_intensity,
    decompose_model_components,
    RHO2_CALIBRATED
)
from models.fitting_engine import (
    create_merged_params,
    residuals_func
)
from lmfit import minimize, fit_report


def fit_merged_dataset():
    data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'dataNDi009_fitting_results_merged.csv'))
    output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'output', 'merged'))
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading merged dataset from: {data_path}")
    df = pd.read_csv(data_path)
    q = df['q'].values
    y = df['0'].values
    dy = df['error_0'].values
    dy = np.maximum(dy, 1e-4)

    # Filter strictly positive valid points
    mask = (q >= 0.002) & (q <= 2.5) & (y > 0)
    q_data = q[mask]
    y_data = y[mask]
    dy_data = dy[mask]

    print(f"Fitting merged curve: {len(q_data)} points from q = {q_data.min():.4e} to {q_data.max():.4e} nm^-1")

    # Initial parameter object
    params = create_merged_params()

    # --- 4-Stage Hierarchical Fit ---
    print("\n--- Running Stage 1: Primary Building Block & Fractal Fit (q in [0.008, 0.4]) ---")
    s1_mask = (q_data >= 0.008) & (q_data <= 0.4)
    params['phi'].vary = True
    params['Rmean'].vary = True
    params['stdDev'].vary = True
    params['D'].vary = True
    params['xi'].vary = False
    params['mean_radius_large'].vary = False
    params['std_radius_large'].vary = False
    params['total_volume_fraction_large'].vary = False
    params['A_OZ'].vary = False
    params['xi_OZ'].vary = False

    out1 = minimize(residuals_func, params, args=(q_data[s1_mask], y_data[s1_mask], dy_data[s1_mask], False, False), method='least_squares')
    params = out1.params.copy()

    print("\n--- Running Stage 2: Large Spheres Fit at Low-q (q < 0.006) ---")
    s2_mask = (q_data < 0.006)
    params['phi'].vary = False
    params['Rmean'].vary = False
    params['stdDev'].vary = False
    params['D'].vary = False
    params['xi'].vary = True
    params['mean_radius_large'].vary = True
    params['std_radius_large'].vary = True
    params['total_volume_fraction_large'].vary = True
    params['A_OZ'].vary = False
    params['xi_OZ'].vary = False

    out2 = minimize(residuals_func, params, args=(q_data[s2_mask], y_data[s2_mask], dy_data[s2_mask], True, False, 3.0), method='least_squares')
    params = out2.params.copy()

    print("\n--- Running Stage 3: Ornstein-Zernike Excess Fit at High-q (q > 0.3) ---")
    s3_mask = (q_data > 0.3)
    params['mean_radius_large'].vary = False
    params['std_radius_large'].vary = False
    params['total_volume_fraction_large'].vary = False
    params['xi'].vary = False
    params['A_OZ'].vary = True
    params['xi_OZ'].vary = True

    out3 = minimize(residuals_func, params, args=(q_data[s3_mask], y_data[s3_mask], dy_data[s3_mask], False, True), method='least_squares')
    params = out3.params.copy()

    print("\n--- Running Stage 4: Full Multi-Component Global Refinement ---")
    params['phi'].vary = True
    params['Rmean'].vary = True
    params['stdDev'].vary = True
    params['D'].vary = True
    params['xi'].vary = True
    params['mean_radius_large'].vary = True
    params['std_radius_large'].vary = True
    params['total_volume_fraction_large'].vary = True
    params['A_OZ'].vary = True
    params['xi_OZ'].vary = True

    out_final = minimize(residuals_func, params, args=(q_data, y_data, dy_data, True, True, 2.0), method='least_squares', max_nfev=800)
    final_params = out_final.params

    red_chi2 = out_final.redchi
    print("\n" + "=" * 60)
    print(f"Merged Baseline Fit Complete: Reduced Chi-Square = {red_chi2:.4f}")
    print("=" * 60)

    # Save parameters to dict and files
    param_dict = {}
    for p_name in final_params:
        p_val = float(final_params[p_name].value)
        p_err = float(final_params[p_name].stderr) if final_params[p_name].stderr is not None else 0.0
        param_dict[p_name] = p_val
        param_dict[f"std_err_{p_name}"] = p_err
        print(f"  {p_name:<28}: {p_val:12.6f} ± {p_err:10.6f}")

    param_dict['redchi'] = float(red_chi2)
    param_dict['aic'] = float(out_final.aic)

    # Save to both output/composite/ (standardized) and output/merged/ (legacy compatibility)
    composite_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'output', 'composite'))
    os.makedirs(composite_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    for target_dir in [composite_dir, output_dir]:
        json_path = os.path.join(target_dir, 'baseline_parameters.json')
        with open(json_path, 'w') as f:
            json.dump(param_dict, f, indent=4)
        print(f"Saved baseline parameters to: {json_path}")

        csv_path = os.path.join(target_dir, 'baseline_parameters.csv')
        pd.DataFrame([param_dict]).to_csv(csv_path, index=False)
        print(f"Saved baseline parameters to: {csv_path}")

    # --- Plot Baseline Fit & Decomposition ---
    print("\nGenerating baseline fit visualization...")
    q_dense = np.logspace(np.log10(q_data.min()), np.log10(q_data.max()), 600)
    I_tot, I_prim, I_large, I_oz, I_bkg = decompose_model_components(final_params, q_dense)

    fig = plt.figure(figsize=(11, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.1)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    # Top panel: Data & Model Decompositions
    ax_fit.loglog(q_data, y_data, 'o', color='#2b2d42', markersize=2, alpha=0.4, label='ND-i-009 Composite Data')
    ax_fit.loglog(q_dense, I_tot, 'r-', linewidth=2.0, label='Full Composite Fit')
    ax_fit.loglog(q_dense, I_prim, 'g--', linewidth=1.5, label=f'Primary Fractal (R={final_params["Rmean"].value:.2f} nm, D={final_params["D"].value:.2f})')
    ax_fit.loglog(q_dense, I_large, 'b-.', linewidth=1.5, label=f'Large Droplets (R={final_params["mean_radius_large"].value:.0f} nm)')
    ax_fit.loglog(q_dense, I_oz, 'm:', linewidth=1.5, label='OZ Fluctuations')

    ax_fit.set_ylabel("Scattering Intensity $I(q)$ [cm$^{-1}$]", fontsize=12)
    ax_fit.set_title(f"ND-i-009 Full-Range Composite Scattering Curve Fit (Model Condition, $\\chi^2_{{\\mathrm{{red}}}} = {red_chi2:.3f}$)", fontsize=13)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_fit.legend(loc='upper right', frameon=True, fontsize=10)
    ax_fit.set_ylim(bottom=1e-1, top=1e7)

    # Bottom panel: Residuals
    y_calc_data = calc_model_intensity(final_params, q_data)
    residuals = (y_data - y_calc_data) / dy_data
    ax_res.semilogx(q_data, residuals, 'o', color='#457b9d', markersize=2, alpha=0.5)
    ax_res.axhline(0, color='red', linestyle='--', linewidth=1)
    ax_res.set_xlabel("Scattering Vector $q$ [nm$^{-1}$]", fontsize=12)
    ax_res.set_ylabel("Residuals [$\\sigma$]", fontsize=11)
    ax_res.set_ylim(-5, 5)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)

    # Save to both composite and merged directories
    plt.savefig(os.path.join(composite_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(composite_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'merged_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved baseline fit plots to: {composite_dir} and {output_dir}")

    return param_dict


if __name__ == '__main__':
    fit_merged_dataset()

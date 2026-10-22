#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Scientific Evaluation: Alternative Particle Form Factors for ND-i-005 Short-Distance SAXS
-----------------------------------------------------------------------------------------
Investigates whether non-spherical primary building blocks (platelets/disks, cylinders/rods,
ellipsoids) can fit the short-distance data without a viscoelastic fluctuation component.

Tested Models:
1. Spheres (Unconstrained & Core-Conserved) + Fractal S(q) [Baseline]
2. Platelets / Disks (Thickness T, Radius R_disk) + Fractal S(q)
3. Cylinders / Rods (Length L, Radius R_cyl) + Fractal S(q)
4. Oblate / Prolate Ellipsoids + Fractal S(q)
5. Viscoelastic Spherical Model (Reference)
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
from scipy.special import j1, gamma
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
VISCO_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')

EVAL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(EVAL_DIR, 'output')
FITS_DIR = os.path.join(OUT_DIR, 'fits')

# Calibrated contrast factor
RHO2_SPHERE = 26.884225
DELTA_RHO_SQ = RHO2_SPHERE / (16.0 * (np.pi**2) / 9.0)  # = 1.5322 nm^-6

# 40-point Gauss-Legendre quadrature for orientational averaging
X_NODES, W_NODES = np.polynomial.legendre.leggauss(40)
X_NODES = 0.5 * (X_NODES + 1.0)
W_NODES = 0.5 * W_NODES


def cylinder_disk_form_factor(q, R, L):
    """
    Orientation-averaged form factor P(q) for cylinder/disk of radius R and length/thickness L.
    P(0) = 1.0.
    """
    q_mat = np.asarray(q)[:, None]
    x_mat = X_NODES[None, :]
    
    alpha = q_mat * R * np.sqrt(np.maximum(0.0, 1.0 - x_mat**2))
    beta = 0.5 * q_mat * L * x_mat
    
    term1 = np.ones_like(alpha)
    mask1 = alpha > 1e-6
    term1[mask1] = 2.0 * j1(alpha[mask1]) / alpha[mask1]
    
    term2 = np.ones_like(beta)
    mask2 = beta > 1e-6
    term2[mask2] = np.sin(beta[mask2]) / beta[mask2]
    
    integrand = (term1 * term2)**2
    return np.sum(integrand * W_NODES[None, :], axis=1)


def sphere_form_factor(q, R):
    """Normalized sphere form factor."""
    qr = np.asarray(q)[:, None] * np.asarray(R)[None, :]
    val = np.ones_like(qr)
    mask = qr > 1e-5
    qr_m = qr[mask]
    val[mask] = (3.0 * (np.sin(qr_m) - qr_m * np.cos(qr_m)) / (qr_m**3))**2
    return val


def mass_fractal_sq(q, R_eff, Dm, xi):
    """Teixeira mass fractal structure factor."""
    if xi <= 0 or Dm <= 1.0:
        return np.ones_like(q)
    qr = q * R_eff
    with np.errstate(divide='ignore', invalid='ignore'):
        term1 = Dm * gamma(Dm - 1.0) / (qr**Dm)
        if np.isinf(xi) or xi >= 1e6:
            term2 = np.sin((Dm - 1.0) * np.pi / 2.0)
            Sq = 1.0 + term1 * term2
        else:
            q_xi = q * xi
            term2 = np.sin((Dm - 1.0) * np.arctan(q_xi)) / ((1.0 + q_xi**(-2))**((Dm - 1.0) / 2.0))
            Sq = 1.0 + term1 * term2
    Sq = np.where(np.isnan(Sq) | np.isinf(Sq), 1.0, Sq)
    return np.maximum(Sq, 0.0)


def calc_platelet_model(params, q, xi_fixed=1193.58):
    """Integrated intensity for polydisperse platelets with mass-fractal structure factor."""
    phi = params['phi'].value
    R_disk = params['R_disk'].value
    sigma_R = params['sigma_R'].value
    T = params['thickness'].value
    Dm = params['D'].value
    bkg = params['bkg'].value
    
    if phi <= 0 or R_disk <= 0 or T <= 0:
        return np.full_like(q, bkg)
        
    mu = np.log(R_disk**2 / np.sqrt(sigma_R**2 + R_disk**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_disk**2)))
    
    r_vals = np.linspace(max(0.5, R_disk - 2.5 * sigma_R), R_disk + 3.0 * sigma_R, 15)
    dr = r_vals[1] - r_vals[0]
    pdf = (1.0 / (r_vals * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(r_vals) - mu)**2) / (2 * sigma**2))
    
    v_vals = np.pi * (r_vals**2) * T
    v_mean = np.sum(pdf * v_vals * dr)
    if v_mean <= 0:
        return np.full_like(q, bkg)
    n_density = phi / v_mean
    
    # Effective radius for fractal structure factor S(q)
    R_eff = (3.0 * v_mean / (4.0 * np.pi))**(1.0 / 3.0)
    Sq = mass_fractal_sq(q, R_eff, Dm, xi_fixed)
    
    I_platelet = np.zeros_like(q, dtype=float)
    for r, v, p in zip(r_vals, v_vals, pdf):
        P = cylinder_disk_form_factor(q, r, T)
        I_platelet += (v**2) * P * p * dr
        
    return DELTA_RHO_SQ * n_density * I_platelet * Sq + bkg


def calc_cylinder_model(params, q, xi_fixed=1193.58):
    """Integrated intensity for polydisperse cylinders with mass-fractal structure factor."""
    phi = params['phi'].value
    R_cyl = params['R_cyl'].value
    sigma_R = params['sigma_R'].value
    L = params['length'].value
    Dm = params['D'].value
    bkg = params['bkg'].value
    
    if phi <= 0 or R_cyl <= 0 or L <= 0:
        return np.full_like(q, bkg)
        
    mu = np.log(R_cyl**2 / np.sqrt(sigma_R**2 + R_cyl**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_cyl**2)))
    
    r_vals = np.linspace(max(0.2, R_cyl - 2.5 * sigma_R), R_cyl + 3.0 * sigma_R, 15)
    dr = r_vals[1] - r_vals[0]
    pdf = (1.0 / (r_vals * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(r_vals) - mu)**2) / (2 * sigma**2))
    
    v_vals = np.pi * (r_vals**2) * L
    v_mean = np.sum(pdf * v_vals * dr)
    if v_mean <= 0:
        return np.full_like(q, bkg)
    n_density = phi / v_mean
    
    R_eff = (3.0 * v_mean / (4.0 * np.pi))**(1.0 / 3.0)
    Sq = mass_fractal_sq(q, R_eff, Dm, xi_fixed)
    
    I_cyl = np.zeros_like(q, dtype=float)
    for r, v, p in zip(r_vals, v_vals, pdf):
        P = cylinder_disk_form_factor(q, r, L)
        I_cyl += (v**2) * P * p * dr
        
    return DELTA_RHO_SQ * n_density * I_cyl * Sq + bkg


def calc_sphere_fractal_model(params, q, xi_fixed=1193.58):
    """Integrated intensity for polydisperse spheres with mass-fractal structure factor."""
    phi = params['phi'].value
    R_mean = params['Rmean'].value
    sigma_R = params['stdDev'].value
    Dm = params['D'].value
    bkg = params['bkg'].value
    
    if phi <= 0 or R_mean <= 0:
        return np.full_like(q, bkg)
        
    mu = np.log(R_mean**2 / np.sqrt(sigma_R**2 + R_mean**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_mean**2)))
    
    r_vals = np.linspace(max(0.1, R_mean - 2.8 * sigma_R), R_mean + 3.5 * sigma_R, 25)
    dr = r_vals[1] - r_vals[0]
    pdf = (1.0 / (r_vals * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(r_vals) - mu)**2) / (2 * sigma**2))
    
    v_vals = (4.0 / 3.0) * np.pi * (r_vals**3)
    v_mean = np.sum(pdf * v_vals * dr)
    n_density = phi / v_mean
    
    Sq = mass_fractal_sq(q, R_mean, Dm, xi_fixed)
    P = sphere_form_factor(q, r_vals)
    integrand = (r_vals[None, :]**6) * P * pdf[None, :] * dr
    return RHO2_SPHERE * n_density * np.sum(integrand, axis=1) * Sq + bkg


def calc_viscoelastic_reference(params, q, xi_fixed=1193.58):
    """Reference viscoelastic model (spheres + OZ-like dynamic fluctuations)."""
    phi = params['phi'].value
    R_mean = params['Rmean'].value
    sigma_R = params['stdDev'].value
    Dm = params['D'].value
    bkg = params['bkg'].value
    A_fluc = params['A_fluc'].value
    xi_fluc = params['xi_fluc'].value
    
    # Spherical fractal
    I_prim = calc_sphere_fractal_model(params, q, xi_fixed) - bkg
    # Fluctuation component
    I_fluc = A_fluc / ((1.0 + (q * xi_fluc)**2)**(Dm / 2.0))
    return I_prim + I_fluc + bkg


def run_full_study():
    print("=" * 80)
    print("  SCIENTIFIC INVESTIGATION: ALTERNATIVE PARTICLE FORM FACTORS (ND-i-005)")
    print("  Testing Platelets (Disks) and Cylinders (Rods) vs Viscoelastic Model")
    print("=" * 80)
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(FITS_DIR, exist_ok=True)
    
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)
    
    with open(BASELINE_JSON, 'r') as f:
        bp = json.load(f)
    xi_net = float(bp['xi'])
    
    df_visco = pd.read_csv(VISCO_CSV)
    
    sample_frames = [0, 5, 10, 20, 30, 40, 49]
    records = []
    
    for f_idx in sample_frames:
        t_sec = float(df_meta.loc[df_meta['frame_index'] == f_idx, 'time_seconds'].iloc[0])
        y = df_data[str(f_idx)].values
        dy = np.maximum(df_data['error_' + str(f_idx)].values, 1e-4)
        
        mask = (q >= 0.035) & (q <= 2.2) & (y > 0)
        q_v = q[mask]
        y_v = y[mask]
        dy_v = dy[mask]
        
        # 1. Viscoelastic Reference
        row_v = df_visco[df_visco['frame_index'] == f_idx].iloc[0]
        p_visco = Parameters()
        p_visco.add('phi', value=float(row_v['phi']))
        p_visco.add('Rmean', value=float(row_v['Rmean']))
        p_visco.add('stdDev', value=float(row_v['stdDev']))
        p_visco.add('D', value=float(row_v['D']))
        p_visco.add('bkg', value=float(row_v['bkg']))
        p_visco.add('A_fluc', value=float(row_v['A_fluc']))
        p_visco.add('xi_fluc', value=float(row_v['xi_fluc']))
        calc_v = calc_viscoelastic_reference(p_visco, q_v, xi_net)
        res_v = (y_v - calc_v) / dy_v
        chi2_visco = float(row_v['chi2_red'])
        
        # 2. Polydisperse Spheres (Unconstrained)
        p_sph = Parameters()
        p_sph.add('phi', value=float(row_v['phi']), min=1e-7, max=0.01)
        p_sph.add('Rmean', value=2.0, min=0.5, max=25.0)
        p_sph.add('stdDev', value=1.2, min=0.1, max=10.0)
        p_sph.add('D', value=float(row_v['D']), min=1.5, max=2.8)
        p_sph.add('bkg', value=0.001, min=1e-5, max=0.02)
        out_sph = minimize(lambda p: (y_v - calc_sphere_fractal_model(p, q_v, xi_net)) / dy_v, p_sph, method='least_squares')
        calc_sph = calc_sphere_fractal_model(out_sph.params, q_v, xi_net)
        res_sph = (y_v - calc_sph) / dy_v
        
        # 3. Platelets / Disks (Thickness T, Lateral Radius R_disk)
        p_dsk = Parameters()
        p_dsk.add('phi', value=float(row_v['phi']), min=1e-7, max=0.01)
        p_dsk.add('R_disk', value=12.0, min=2.0, max=60.0)
        p_dsk.add('sigma_R', value=3.0, min=0.5, max=20.0)
        p_dsk.add('thickness', value=2.0, min=0.4, max=10.0)
        p_dsk.add('D', value=float(row_v['D']), min=1.5, max=2.8)
        p_dsk.add('bkg', value=0.0005, min=1e-5, max=0.02)
        out_dsk = minimize(lambda p: (y_v - calc_platelet_model(p, q_v, xi_net)) / dy_v, p_dsk, method='least_squares')
        calc_dsk = calc_platelet_model(out_dsk.params, q_v, xi_net)
        res_dsk = (y_v - calc_dsk) / dy_v
        
        # 4. Cylinders / Rods (Length L, Radius R_cyl)
        p_cyl = Parameters()
        p_cyl.add('phi', value=float(row_v['phi']), min=1e-7, max=0.01)
        p_cyl.add('R_cyl', value=2.0, min=0.5, max=12.0)
        p_cyl.add('sigma_R', value=0.6, min=0.1, max=5.0)
        p_cyl.add('length', value=30.0, min=5.0, max=150.0)
        p_cyl.add('D', value=float(row_v['D']), min=1.5, max=2.8)
        p_cyl.add('bkg', value=0.0005, min=1e-5, max=0.02)
        out_cyl = minimize(lambda p: (y_v - calc_cylinder_model(p, q_v, xi_net)) / dy_v, p_cyl, method='least_squares')
        calc_cyl = calc_cylinder_model(out_cyl.params, q_v, xi_net)
        res_cyl = (y_v - calc_cyl) / dy_v
        
        # Compute regional residual statistics
        def get_stat(res):
            l = np.mean(res[q_v < 0.1]**2)
            m = np.mean(res[(q_v >= 0.1) & (q_v < 0.5)]**2)
            h = np.mean(res[q_v >= 0.5]**2)
            mx = np.max(np.abs(res))
            return l, m, h, mx
            
        vl, vm, vh, vmx = get_stat(res_v)
        sl, sm, sh, smx = get_stat(res_sph)
        dl, dm, dh, dmx = get_stat(res_dsk)
        cl, cm, ch, cmx = get_stat(res_cyl)
        
        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'chi2_visco': chi2_visco,
            'chi2_spheres': out_sph.redchi,
            'chi2_platelets': out_dsk.redchi,
            'chi2_cylinders': out_cyl.redchi,
            'sph_R': out_sph.params['Rmean'].value,
            'sph_std': out_sph.params['stdDev'].value,
            'sph_bkg': out_sph.params['bkg'].value,
            'sph_low_res2': sl, 'sph_mid_res2': sm, 'sph_high_res2': sh, 'sph_max_res': smx,
            'dsk_R': out_dsk.params['R_disk'].value,
            'dsk_T': out_dsk.params['thickness'].value,
            'dsk_bkg': out_dsk.params['bkg'].value,
            'dsk_low_res2': dl, 'dsk_mid_res2': dm, 'dsk_high_res2': dh, 'dsk_max_res': dmx,
            'cyl_R': out_cyl.params['R_cyl'].value,
            'cyl_L': out_cyl.params['length'].value,
            'cyl_bkg': out_cyl.params['bkg'].value,
            'cyl_low_res2': cl, 'cyl_mid_res2': cm, 'cyl_high_res2': ch, 'cyl_max_res': cmx,
        }
        records.append(rec)
        
        print(f"\nFrame {f_idx:02d} (t={t_sec:.3f}s):")
        print(f"  Viscoelastic Ref: chi2 = {chi2_visco:.3f}")
        print(f"  Spheres:          chi2 = {out_sph.redchi:.3f} | R = {out_sph.params['Rmean'].value:.2f} nm, bkg = {out_sph.params['bkg'].value:.5f} | low_res = {sl:.2f}, high_res = {sh:.2f}")
        print(f"  Platelets:        chi2 = {out_dsk.redchi:.3f} | R_dsk = {out_dsk.params['R_disk'].value:.2f} nm, T = {out_dsk.params['thickness'].value:.2f} nm, bkg = {out_dsk.params['bkg'].value:.5f} | low_res = {dl:.2f}, high_res = {dh:.2f}")
        print(f"  Cylinders:        chi2 = {out_cyl.redchi:.3f} | R_cyl = {out_cyl.params['R_cyl'].value:.2f} nm, L = {out_cyl.params['length'].value:.2f} nm, bkg = {out_cyl.params['bkg'].value:.5f} | low_res = {cl:.2f}, high_res = {ch:.2f}")

        # Generate individual frame comparative plot
        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8.5, 7.5), gridspec_kw={'height_ratios': [2.2, 1.0], 'hspace': 0.12}, sharex=True)
        ax_top.set_title(f"ND-i-005 Frame {f_idx:02d} (t = {t_sec:.3f} s): Particle Shape Form Factor Comparison", fontsize=11, fontweight='bold')
        
        ax_top.errorbar(q_v, y_v, yerr=dy_v, fmt='.', color='#94a3b8', markersize=3.5, alpha=0.6, label='Observed Data', zorder=1)
        ax_top.plot(q_v, calc_v, color='#059669', linewidth=2.2, label=f'Viscoelastic Reference ($\\chi^2 = {chi2_visco:.3f}$)', zorder=5)
        ax_top.plot(q_v, calc_sph, color='#d97706', linewidth=1.8, linestyle='--', label=f'Polydisperse Spheres ($\\chi^2 = {out_sph.redchi:.3f}$, $R={out_sph.params["Rmean"].value:.1f}$ nm)', zorder=4)
        ax_top.plot(q_v, calc_dsk, color='#dc2626', linewidth=1.8, linestyle='-.', label=f'Platelets/Disks ($\\chi^2 = {out_dsk.redchi:.3f}$, $R={out_dsk.params["R_disk"].value:.1f}, T={out_dsk.params["thickness"].value:.1f}$ nm)', zorder=3)
        ax_top.plot(q_v, calc_cyl, color='#2563eb', linewidth=1.8, linestyle=':', label=f'Cylinders/Rods ($\\chi^2 = {out_cyl.redchi:.3f}$, $R={out_cyl.params["R_cyl"].value:.1f}, L={out_cyl.params["length"].value:.1f}$ nm)', zorder=2)
        
        ax_top.set_xscale('log')
        ax_top.set_yscale('log')
        ax_top.set_ylabel(r"$I(q)$ [cm$^{-1}$]", fontsize=10, fontweight='bold')
        ax_top.grid(True, which='both', linestyle='--', alpha=0.4)
        ax_top.legend(loc='lower left', fontsize=8.5, framealpha=0.92)
        
        # Residuals
        ax_bot.axhline(0, color='black', linewidth=1.0)
        ax_bot.plot(q_v, res_v, color='#059669', linewidth=1.2, label='Viscoelastic')
        ax_bot.plot(q_v, res_sph, color='#d97706', linewidth=1.2, linestyle='--', label='Spheres')
        ax_bot.plot(q_v, res_dsk, color='#dc2626', linewidth=1.2, linestyle='-.', label='Platelets')
        ax_bot.plot(q_v, res_cyl, color='#2563eb', linewidth=1.2, linestyle=':', label='Cylinders')
        
        ax_bot.set_xscale('log')
        ax_bot.set_ylabel(r"$(I_{\mathrm{obs}} - I_{\mathrm{calc}}) / \sigma$", fontsize=9, fontweight='bold')
        ax_bot.set_xlabel(r"$q$ [nm$^{-1}$]", fontsize=10, fontweight='bold')
        ax_bot.set_ylim(-6.5, 6.5)
        ax_bot.grid(True, which='both', linestyle='--', alpha=0.4)
        ax_bot.legend(loc='upper right', fontsize=8.0, ncol=2)
        
        plt.tight_layout()
        plot_path = os.path.join(FITS_DIR, f"shape_comparison_frame_{f_idx:02d}.png")
        plt.savefig(plot_path, dpi=200)
        plt.close(fig)
        
    df_results = pd.DataFrame(records)
    summary_path = os.path.join(OUT_DIR, 'alternative_form_factors_summary.csv')
    df_results.to_csv(summary_path, index=False)
    print(f"\nSaved summary to: {summary_path}")
    
    # Master comparison plot across sample frames
    generate_master_figure(df_results, q_v, y_v, dy_v, calc_v, calc_sph, calc_dsk, calc_cyl, res_v, res_sph, res_dsk, res_cyl)


def generate_master_figure(df_results, q_v, y_v, dy_v, calc_v, calc_sph, calc_dsk, calc_cyl, res_v, res_sph, res_dsk, res_cyl):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    fig.suptitle("ND-i-005 Form Factor Evaluation: Platelets & Cylinders vs Viscoelastic Phase Separation", fontsize=13, fontweight='bold', y=0.98)
    
    # Panel A: Frame 00 Model Fits
    ax1 = axes[0, 0]
    ax1.errorbar(q_v, y_v, yerr=dy_v, fmt='.', color='#94a3b8', markersize=3, alpha=0.5, label='Frame 00 Data')
    ax1.plot(q_v, calc_v, color='#059669', linewidth=2.0, label='Viscoelastic Model')
    ax1.plot(q_v, calc_sph, color='#d97706', linewidth=1.6, linestyle='--', label='Polydisperse Spheres')
    ax1.plot(q_v, calc_dsk, color='#dc2626', linewidth=1.6, linestyle='-.', label='Platelets / Disks')
    ax1.plot(q_v, calc_cyl, color='#2563eb', linewidth=1.6, linestyle=':', label='Cylinders / Rods')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r"$I(q)$ [cm$^{-1}$]", fontsize=10, fontweight='bold')
    ax1.set_xlabel(r"$q$ [nm$^{-1}$]", fontsize=10, fontweight='bold')
    ax1.set_title(r"(A) Frame 00 Fit Profiles across Shapes ($q = 0.035 - 2.2$ nm$^{-1}$)", fontsize=11, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.4)
    ax1.legend(loc='lower left', fontsize=8.5)
    
    # Panel B: Frame 00 Residuals
    ax2 = axes[0, 1]
    ax2.axhline(0, color='black', linewidth=1.0)
    ax2.plot(q_v, res_v, color='#059669', linewidth=1.2, label='Viscoelastic')
    ax2.plot(q_v, res_sph, color='#d97706', linewidth=1.2, linestyle='--', label='Spheres')
    ax2.plot(q_v, res_dsk, color='#dc2626', linewidth=1.2, linestyle='-.', label='Platelets')
    ax2.plot(q_v, res_cyl, color='#2563eb', linewidth=1.2, linestyle=':', label='Cylinders')
    ax2.set_xscale('log')
    ax2.set_ylabel(r"$(I_{\mathrm{obs}} - I_{\mathrm{calc}}) / \sigma$", fontsize=10, fontweight='bold')
    ax2.set_xlabel(r"$q$ [nm$^{-1}$]", fontsize=10, fontweight='bold')
    ax2.set_ylim(-6.0, 6.0)
    ax2.set_title(r"(B) Frame 00 Weighted Residuals (Systematic Wave Test)", fontsize=11, fontweight='bold')
    ax2.grid(True, which='both', linestyle='--', alpha=0.4)
    ax2.legend(loc='upper right', fontsize=8.5, ncol=2)
    
    # Panel C: Chi2 Comparison Across Time
    ax3 = axes[1, 0]
    t = df_results['time_seconds']
    ax3.plot(t, df_results['chi2_visco'], 'o-', color='#059669', linewidth=2.0, label='Viscoelastic Model')
    ax3.plot(t, df_results['chi2_spheres'], 's--', color='#d97706', linewidth=1.6, label='Polydisperse Spheres')
    ax3.plot(t, df_results['chi2_platelets'], '^-.', color='#dc2626', linewidth=1.6, label='Platelets / Disks')
    ax3.plot(t, df_results['chi2_cylinders'], 'd:', color='#2563eb', linewidth=1.6, label='Cylinders / Rods')
    ax3.set_xscale('log')
    ax3.set_ylabel(r"Reduced $\chi^2$", fontsize=10, fontweight='bold')
    ax3.set_xlabel("Reaction Time $t$ [s]", fontsize=10, fontweight='bold')
    ax3.set_title(r"(C) Goodness-of-Fit $\chi^2_{\mathrm{red}}$ vs Time", fontsize=11, fontweight='bold')
    ax3.grid(True, which='both', linestyle='--', alpha=0.4)
    ax3.legend(loc='upper right', fontsize=8.5)
    
    # Panel D: Low-q Residual Misfit
    ax4 = axes[1, 1]
    ax4.plot(t, np.sqrt(df_results['sph_low_res2']), 's--', color='#d97706', linewidth=1.6, label='Spheres Low-q RMS Error')
    ax4.plot(t, np.sqrt(df_results['dsk_low_res2']), '^-.', color='#dc2626', linewidth=1.6, label='Platelets Low-q RMS Error')
    ax4.plot(t, np.sqrt(df_results['cyl_low_res2']), 'd:', color='#2563eb', linewidth=1.6, label='Cylinders Low-q RMS Error')
    ax4.axhline(1.0, color='#059669', linestyle='-', linewidth=2.0, label=r'Ideal Statistical Noise ($1\sigma$)')
    ax4.set_xscale('log')
    ax4.set_ylabel(r"Low-$q$ RMS Residual [$\sigma$]", fontsize=10, fontweight='bold')
    ax4.set_xlabel("Reaction Time $t$ [s]", fontsize=10, fontweight='bold')
    ax4.set_title(r"(D) Low-$q$ Primary Knee Distortion ($\langle \Delta^2 \rangle_{q < 0.1}^{1/2}$)", fontsize=11, fontweight='bold')
    ax4.grid(True, which='both', linestyle='--', alpha=0.4)
    ax4.legend(loc='upper right', fontsize=8.5)
    
    plt.tight_layout()
    master_path = os.path.join(OUT_DIR, 'alternative_form_factors_comparison.png')
    plt.savefig(master_path, dpi=200)
    plt.close(fig)
    print(f"Master comparative figure saved to: {master_path}")


if __name__ == '__main__':
    run_full_study()

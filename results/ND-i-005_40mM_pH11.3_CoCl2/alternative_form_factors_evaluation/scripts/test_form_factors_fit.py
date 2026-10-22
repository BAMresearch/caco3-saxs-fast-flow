#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script to evaluate alternative form factors (Platelets/Disks, Cylinders/Rods, Ellipsoids)
against Polydisperse Spheres and the Viscoelastic model on ND-i-005 short distance data.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.special import j1
from lmfit import minimize, Parameters

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
RHO2 = 2.0e-5

# 40-point Gauss-Legendre quadrature for orientation averaging
X_NODES, W_NODES = np.polynomial.legendre.leggauss(40)
X_NODES = 0.5 * (X_NODES + 1.0)
W_NODES = 0.5 * W_NODES

def cylinder_P(q, R, L):
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

def calc_cylinder_disk_model(params, q, is_disk=False):
    phi = params['phi'].value
    R_mean = params['R_mean'].value
    sigma_R = params['sigma_R'].value
    L = params['thickness'].value if is_disk else params['length'].value
    bkg = params['bkg'].value
    
    if phi <= 0 or R_mean <= 0 or L <= 0:
        return np.full_like(q, bkg)
        
    mu = np.log(R_mean**2 / np.sqrt(sigma_R**2 + R_mean**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_mean**2)))
    
    r_min = max(0.2, R_mean - 2.5 * sigma_R)
    r_max = R_mean + 3.0 * sigma_R
    r_vals = np.linspace(r_min, r_max, 15)
    dr = r_vals[1] - r_vals[0]
    
    pdf = (1.0 / (r_vals * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(r_vals) - mu)**2) / (2 * sigma**2))
    v_vals = np.pi * (r_vals**2) * L
    v_mean = np.sum(pdf * v_vals * dr)
    
    if v_mean <= 0:
        return np.full_like(q, bkg)
        
    n_density = phi / v_mean
    I_tot = np.zeros_like(q, dtype=float)
    for r, v, p in zip(r_vals, v_vals, pdf):
        P = cylinder_P(q, r, L)
        I_tot += (v**2) * P * p * dr
        
    return RHO2 * n_density * I_tot + bkg

def sphere_P(q, r):
    qr = np.asarray(q)[:, None] * np.asarray(r)[None, :]
    val = np.ones_like(qr)
    mask = qr > 1e-5
    qr_m = qr[mask]
    val[mask] = (3.0 * (np.sin(qr_m) - qr_m * np.cos(qr_m)) / (qr_m**3))**2
    return val

def calc_sphere_model(params, q):
    phi = params['phi'].value
    R_mean = params['Rmean'].value
    sigma_R = params['stdDev'].value
    bkg = params['bkg'].value
    
    if phi <= 0 or R_mean <= 0:
        return np.full_like(q, bkg)
        
    mu = np.log(R_mean**2 / np.sqrt(sigma_R**2 + R_mean**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_mean**2)))
    
    r_min = max(0.1, R_mean - 2.8 * sigma_R)
    r_max = R_mean + 3.5 * sigma_R
    r_vals = np.linspace(r_min, r_max, 25)
    dr = r_vals[1] - r_vals[0]
    
    pdf = (1.0 / (r_vals * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(r_vals) - mu)**2) / (2 * sigma**2))
    v_vals = (4.0 / 3.0) * np.pi * (r_vals**3)
    v_mean = np.sum(pdf * v_vals * dr)
    
    n_density = phi / v_mean
    P = sphere_P(q, r_vals)
    integrand = (v_vals[None, :]**2) * P * pdf[None, :] * dr
    return RHO2 * n_density * np.sum(integrand, axis=1) + bkg

def run_comparative_fit_frame(frame_idx=0):
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    y = df_data[str(frame_idx)].values
    dy = np.maximum(df_data['error_' + str(frame_idx)].values, 1e-4)
    
    mask = (q >= 0.035) & (q <= 2.2) & (y > 0)
    q_f = q[mask]
    y_f = y[mask]
    dy_f = dy[mask]
    
    print(f"\n================ FRAME {frame_idx:02d} FITTING COMPARISON ================")
    
    # 1. Polydisperse Spheres (Unconstrained)
    p_sph = Parameters()
    p_sph.add('phi', value=2e-4, min=1e-6, max=0.01)
    p_sph.add('Rmean', value=2.0, min=0.5, max=25.0)
    p_sph.add('stdDev', value=1.0, min=0.1, max=10.0)
    p_sph.add('bkg', value=5e-4, min=0.0, max=0.01)
    
    def res_sph(p):
        return (y_f - calc_sphere_model(p, q_f)) / dy_f
        
    out_sph = minimize(res_sph, p_sph, method='least_squares')
    chi2_sph = out_sph.redchi
    print(f"1. Polydisperse Spheres (Unconstrained): chi2_red = {chi2_sph:.3f} | R = {out_sph.params['Rmean'].value:.2f} nm, sigma = {out_sph.params['stdDev'].value:.2f} nm, bkg = {out_sph.params['bkg'].value:.5f}")

    # 2. Platelets / Disks (Thin circular disks: R_disk >> Thickness T)
    p_disk = Parameters()
    p_disk.add('phi', value=2e-4, min=1e-6, max=0.01)
    p_disk.add('R_mean', value=10.0, min=2.0, max=80.0)  # Disk radius
    p_disk.add('sigma_R', value=3.0, min=0.5, max=30.0)
    p_disk.add('thickness', value=1.5, min=0.3, max=8.0) # Disk thickness T
    p_disk.add('bkg', value=5e-4, min=0.0, max=0.01)
    
    def res_disk(p):
        return (y_f - calc_cylinder_disk_model(p, q_f, is_disk=True)) / dy_f
        
    out_disk = minimize(res_disk, p_disk, method='least_squares')
    chi2_disk = out_disk.redchi
    print(f"2. Platelets / Disks:                  chi2_red = {chi2_disk:.3f} | R_disk = {out_disk.params['R_mean'].value:.2f} nm, T = {out_disk.params['thickness'].value:.2f} nm, bkg = {out_disk.params['bkg'].value:.5f}")

    # 3. Cylinders / Rods (Elongated needles: Length L >> Radius R)
    p_cyl = Parameters()
    p_cyl.add('phi', value=2e-4, min=1e-6, max=0.01)
    p_cyl.add('R_mean', value=2.0, min=0.5, max=15.0)   # Cross-sectional radius
    p_cyl.add('sigma_R', value=0.6, min=0.1, max=6.0)
    p_cyl.add('length', value=25.0, min=5.0, max=150.0) # Length L
    p_cyl.add('bkg', value=5e-4, min=0.0, max=0.01)
    
    def res_cyl(p):
        return (y_f - calc_cylinder_disk_model(p, q_f, is_disk=False)) / dy_f
        
    out_cyl = minimize(res_cyl, p_cyl, method='least_squares')
    chi2_cyl = out_cyl.redchi
    print(f"3. Cylinders / Rods:                   chi2_red = {chi2_cyl:.3f} | R_cyl = {out_cyl.params['R_mean'].value:.2f} nm, L = {out_cyl.params['length'].value:.2f} nm, bkg = {out_cyl.params['bkg'].value:.5f}")

    # Compute residual breakdown across q-bands
    def calc_residuals(fit_fn, params):
        calc = fit_fn(params, q_f)
        diff = (y_f - calc) / dy_f
        r_low = np.mean(diff[q_f < 0.1]**2)
        r_mid = np.mean(diff[(q_f >= 0.1) & (q_f < 0.5)]**2)
        r_high = np.mean(diff[q_f >= 0.5]**2)
        return r_low, r_mid, r_high, np.max(np.abs(diff))

    r_l_sph, r_m_sph, r_h_sph, max_sph = calc_residuals(calc_sphere_model, out_sph.params)
    r_l_dsk, r_m_dsk, r_h_dsk, max_dsk = calc_residuals(lambda p, q: calc_cylinder_disk_model(p, q, is_disk=True), out_disk.params)
    r_l_cyl, r_m_cyl, r_h_cyl, max_cyl = calc_residuals(lambda p, q: calc_cylinder_disk_model(p, q, is_disk=False), out_cyl.params)

    print("\nResidual Breakdown (<Delta^2> = mean((y - calc)^2 / dy^2)):")
    print(f"  Spheres:   low-q = {r_l_sph:.2f}, mid-q = {r_m_sph:.2f}, high-q = {r_h_sph:.2f}, max_abs = {max_sph:.2f} sigma")
    print(f"  Platelets: low-q = {r_l_dsk:.2f}, mid-q = {r_m_dsk:.2f}, high-q = {r_h_dsk:.2f}, max_abs = {max_dsk:.2f} sigma")
    print(f"  Cylinders: low-q = {r_l_cyl:.2f}, mid-q = {r_m_cyl:.2f}, high-q = {r_h_cyl:.2f}, max_abs = {max_cyl:.2f} sigma")

if __name__ == '__main__':
    run_comparative_fit_frame(0)
    run_comparative_fit_frame(10)
    run_comparative_fit_frame(30)
    run_comparative_fit_frame(49)

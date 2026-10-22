#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Quick exploratory test for alternative form factors on ND-i-005 short distance data:
- Cylinders (Rods / Needles)
- Platelets (Disks / Lamellae)
- Ellipsoids (Oblate / Prolate)
- Spheres (Baseline)
- Viscoelastic Model (Reference)
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.special import j1
from lmfit import minimize, Parameters

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_short_distance.csv')

# Precompute Gauss-Legendre quadrature nodes for orientation averaging
X_NODES, W_NODES = np.polynomial.legendre.leggauss(40)
X_NODES = 0.5 * (X_NODES + 1.0)
W_NODES = 0.5 * W_NODES

def cylinder_P(q, R, L):
    """Orientation-averaged form factor for a cylinder/disk of radius R and length/thickness L."""
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

def calc_polydisperse_cylinder_intensity(q, R_mean, sigma_R, L, phi, bkg, rho2=2.0e-5, n_points=15):
    """Polydisperse cylinders with log-normal radius R, constant length L."""
    if phi <= 0 or R_mean <= 0 or L <= 0:
        return np.full_like(q, bkg)
    
    # Log-normal distribution on R
    mu = np.log(R_mean**2 / np.sqrt(sigma_R**2 + R_mean**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_mean**2)))
    
    r_min = max(0.2, R_mean - 2.8 * sigma_R)
    r_max = R_mean + 3.5 * sigma_R
    r_vals = np.linspace(r_min, r_max, n_points)
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
        
    return rho2 * n_density * I_tot + bkg

def ellipsoid_P(q, R_eq, eps):
    """Orientation-averaged form factor of ellipsoid with semi-axes (R_eq, R_eq, eps*R_eq)."""
    q_mat = np.asarray(q)[:, None]
    x_mat = X_NODES[None, :]
    u = q_mat * R_eq * np.sqrt(1.0 + (eps**2 - 1.0) * (x_mat**2))
    
    phi_u = np.ones_like(u)
    mask = u > 1e-6
    phi_u[mask] = 3.0 * (np.sin(u[mask]) - u[mask] * np.cos(u[mask])) / (u[mask]**3)
    
    return np.sum((phi_u**2) * W_NODES[None, :], axis=1)

def calc_polydisperse_ellipsoid_intensity(q, R_mean, sigma_R, eps, phi, bkg, rho2=2.0e-5, n_points=15):
    """Polydisperse ellipsoids with log-normal equatorial radius R_eq, axis ratio eps."""
    if phi <= 0 or R_mean <= 0 or eps <= 0:
        return np.full_like(q, bkg)
        
    mu = np.log(R_mean**2 / np.sqrt(sigma_R**2 + R_mean**2))
    sigma = np.sqrt(np.log(1.0 + (sigma_R**2 / R_mean**2)))
    
    r_min = max(0.2, R_mean - 2.8 * sigma_R)
    r_max = R_mean + 3.5 * sigma_R
    r_vals = np.linspace(r_min, r_max, n_points)
    dr = r_vals[1] - r_vals[0]
    
    pdf = (1.0 / (r_vals * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(r_vals) - mu)**2) / (2 * sigma**2))
    v_vals = (4.0 / 3.0) * np.pi * (r_vals**3) * eps
    v_mean = np.sum(pdf * v_vals * dr)
    
    if v_mean <= 0:
        return np.full_like(q, bkg)
        
    n_density = phi / v_mean
    
    I_tot = np.zeros_like(q, dtype=float)
    for r, v, p in zip(r_vals, v_vals, pdf):
        P = ellipsoid_P(q, r, eps)
        I_tot += (v**2) * P * p * dr
        
    return rho2 * n_density * I_tot + bkg

print("Test functions compiled successfully.")

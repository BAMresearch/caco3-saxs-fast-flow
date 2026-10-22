#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Core SAXS Scattering Models for ND-i-009.
Includes:
- Polydisperse primary spheres with Teixeira mass-fractal structure factor
- Polydisperse large spheres (secondary precipitation population)
- Ornstein-Zernike (OZ) concentration fluctuations (high-q excess)

Optimized for numerical stability and high-performance multiprocessing.
"""

import numpy as np
from scipy.special import gamma, gammaln

# Calibrated contrast factor for amorphous calcium carbonate in water
# (SLD_ACC - SLD_water)^2 = 5.185^2 = 26.884225
RHO2_CALIBRATED = 26.884225


def sphere_form_factor(q, radius):
    """
    Normalized form factor P(q, R) for a hard sphere of radius R.
    P(0) = 1.0. Vectorized across q and radius dimensions.
    """
    qr = q * radius
    # Handle q*r approaching 0 using Taylor series limit
    with np.errstate(divide='ignore', invalid='ignore'):
        P = (3.0 * (np.sin(qr) - qr * np.cos(qr)) / (qr**3))**2
    return np.where(qr < 1e-5, 1.0, P)


def mass_fractal_structure_factor(q, radius, Dm, xi):
    """
    Teixeira mass-fractal inter-particle structure factor S(q).
    Accounts for cluster size correlation cutoff xi and fractal dimension Dm.
    If xi is infinite or >= 1e6 nm, evaluates the exact infinite-network asymptotic limit:
    S(q) = 1 + Dm * gamma(Dm - 1) * sin((Dm - 1)*pi/2) / (q * radius)^Dm
    """
    if xi <= 0 or Dm <= 1.0:
        return np.ones_like(q * radius)
        
    qr = q * radius
    
    with np.errstate(divide='ignore', invalid='ignore'):
        term1 = Dm * gamma(Dm - 1.0) / (qr**Dm)
        if np.isinf(xi) or xi >= 1e6:
            # Exact analytical limit for infinitely large percolated network (xi -> infinity)
            term2 = np.sin((Dm - 1.0) * np.pi / 2.0)
            Sq = 1.0 + term1 * term2
        else:
            q_xi = q * xi
            term2 = np.sin((Dm - 1.0) * np.arctan(q_xi)) / ((1.0 + q_xi**(-2))**((Dm - 1.0) / 2.0))
            Sq = 1.0 + term1 * term2
        
    # Replace non-finite or low-q singularity with smooth physical limit
    if np.isinf(xi) or xi >= 1e6:
        Sq = np.where(np.isnan(Sq) | np.isinf(Sq), 1.0, Sq)
    else:
        Sq = np.where(np.isnan(Sq) | np.isinf(Sq), 1.0 + Dm * gamma(Dm - 1.0) * (xi / radius)**Dm, Sq)
    return np.maximum(Sq, 0.0)


def schulz_flory_pdf(radius, mean_r, std_r):
    """
    Schulz-Flory (Schulz-Zimm / Gamma) probability density function for particle radius:
    f(r) = (1 / Gamma(z+1)) * ((z+1)/<R>)^(z+1) * r^z * exp(-(z+1)*r/<R>)
    where z = (<R> / std_r)^2 - 1.
    """
    mean_r = max(mean_r, 1e-5)
    std_r = max(std_r, 1e-5)
    z = (mean_r / std_r)**2 - 1.0
    if z <= -0.99:
        z = -0.99
    scaled_r = (z + 1.0) * radius / mean_r
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        log_pdf = (z + 1.0) * np.log(z + 1.0) - gammaln(z + 1.0) - np.log(mean_r) + z * np.log(radius / mean_r) - scaled_r
        pdf = np.exp(log_pdf)
    return np.where(radius <= 0, 0.0, np.nan_to_num(pdf, nan=0.0))


def log_normal_pdf(radius, mu, sigma):
    """
    Log-normal probability density function.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        pdf = (1.0 / (radius * sigma * np.sqrt(2.0 * np.pi))) * np.exp(-((np.log(radius) - mu)**2) / (2.0 * sigma**2))
    return np.where(radius <= 0, 0.0, np.nan_to_num(pdf, nan=0.0))


def convert_mean_std_to_mu_sigma(mean_r, std_r):
    """
    Converts arithmetic mean and standard deviation to log-normal distribution parameters (mu, sigma).
    """
    mean_r = max(mean_r, 1e-5)
    std_r = max(std_r, 1e-5)
    mu = np.log(mean_r**2 / np.sqrt(std_r**2 + mean_r**2))
    sigma = np.sqrt(np.log(1.0 + (std_r**2 / mean_r**2)))
    return mu, sigma


def calc_total_volume_per_unit(mu, sigma, mean_r, std_r, integration_steps=300):
    """
    Calculates unit volume of the polydisperse sphere population.
    Uses numerical integration matching Tomek's established grid representation.
    """
    r_max = mean_r + 10.0 * std_r
    r_grid = np.linspace(1e-5, r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]
    integrand = (4.0 / 3.0) * np.pi * (r_grid**3) * log_normal_pdf(r_grid, mu, sigma)
    vol = np.sum(integrand * dr)
    return max(vol, 1e-12)


def calc_primary_fractal_intensity(q, mean_radius, std_radius, phi, Dm, xi, rho2=RHO2_CALIBRATED, integration_steps=300):
    """
    Vectorized integrated scattering intensity for primary building blocks with fractal structure factor.
    I(q) = rho2 * int [ r^6 * n(r) * P(q, r) * S(q, r, Dm, xi) ] dr
    """
    if phi <= 0:
        return np.zeros_like(q)
        
    mu, sigma = convert_mean_std_to_mu_sigma(mean_radius, std_radius)
    vol_per_unit = calc_total_volume_per_unit(mu, sigma, mean_radius, std_radius, integration_steps)
    number_density = phi / vol_per_unit

    r_max = mean_radius + 10.0 * std_radius
    r_grid = np.linspace(1e-5, r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]

    # Reshape for broadcasting: r_grid along rows (axis 0), q along cols (axis 1)
    r_col = r_grid[:, np.newaxis]
    q_row = np.asarray(q)[np.newaxis, :]

    P_qr = sphere_form_factor(q_row, r_col)
    S_qr = mass_fractal_structure_factor(q_row, r_col, Dm, xi)
    dist = log_normal_pdf(r_col, mu, sigma) * number_density

    integrand = (r_col**6) * dist * P_qr * S_qr * dr
    np.nan_to_num(integrand, copy=False, nan=0.0)
    
    return rho2 * np.sum(integrand, axis=0)


def calc_large_sphere_intensity(q, mean_radius_large, std_radius_large, phi_large, rho2=RHO2_CALIBRATED, integration_steps=300):
    """
    Vectorized integrated scattering intensity for the secondary large sphere population.
    I_large(q) = rho2 * int [ r^6 * n(r) * P(q, r) ] dr
    """
    if phi_large <= 0 or mean_radius_large <= 0:
        return np.zeros_like(q)
        
    mu, sigma = convert_mean_std_to_mu_sigma(mean_radius_large, std_radius_large)
    vol_per_unit = calc_total_volume_per_unit(mu, sigma, mean_radius_large, std_radius_large, integration_steps)
    number_density = phi_large / vol_per_unit

    r_max = mean_radius_large + 10.0 * std_radius_large
    r_grid = np.linspace(1e-5, r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]

    r_col = r_grid[:, np.newaxis]
    q_row = np.asarray(q)[np.newaxis, :]

    P_qr = sphere_form_factor(q_row, r_col)
    dist = log_normal_pdf(r_col, mu, sigma) * number_density

    integrand = (r_col**6) * dist * P_qr * dr
    np.nan_to_num(integrand, copy=False, nan=0.0)
    
    return rho2 * np.sum(integrand, axis=0)


def calc_oz_intensity(q, A_OZ, xi_OZ):
    """
    Ornstein-Zernike concentration fluctuations: I_OZ(q) = A_OZ / (1 + (q * xi_OZ)^2).
    """
    if A_OZ <= 0:
        return np.zeros_like(q)
    return A_OZ / (1.0 + (q * xi_OZ)**2)


def calc_fb_intensity(q, A_fluc, xi_fluc, mu_fluc=4.0):
    """
    Generalized Fisher-Burford / Tanaka viscoelastic concentration fluctuations:
    I_fluc(q) = A_fluc / [ 1 + (2 / (3 * mu)) * (q * xi_fluc)^2 ]^(mu / 2).
    For mu=2: reduces to Ornstein-Zernike.
    For mu=4: Porod boundary decay limit.
    """
    if A_fluc <= 0:
        return np.zeros_like(q)
    return A_fluc / (1.0 + (2.0 / (3.0 * mu_fluc)) * (q * xi_fluc)**2)**(mu_fluc / 2.0)


def calc_pnc_intensity(q, R_pnc, std_pnc, phi_pnc, rho2=RHO2_CALIBRATED, integration_steps=120):
    """
    Vectorized integrated scattering intensity for sub-nanometer prenucleation clusters (DOLLOPs).
    I_pnc(q) = rho2 * int [ r^6 * n(r) * P(q, r) ] dr
    """
    if phi_pnc <= 0 or R_pnc <= 0:
        return np.zeros_like(q)
        
    mu, sigma = convert_mean_std_to_mu_sigma(R_pnc, std_pnc)
    vol_per_unit = calc_total_volume_per_unit(mu, sigma, R_pnc, std_pnc, integration_steps)
    number_density = phi_pnc / vol_per_unit

    r_min = max(1e-4, R_pnc - 3.5 * std_pnc)
    r_max = R_pnc + 6.0 * std_pnc
    r_grid = np.linspace(r_min, r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]

    r_col = r_grid[:, np.newaxis]
    q_row = np.asarray(q)[np.newaxis, :]

    P_qr = sphere_form_factor(q_row, r_col)
    dist = log_normal_pdf(r_col, mu, sigma) * number_density

    integrand = (r_col**6) * dist * P_qr * dr
    np.nan_to_num(integrand, copy=False, nan=0.0)
    
    return rho2 * np.sum(integrand, axis=0)


def calc_porod_intensity(q, A_porod, exp_porod=4.0):
    """
    Macroscopic Porod surface scattering component: I_porod(q) = A_porod / (q ** exp_porod).
    Captures asymptotic surface scattering from giant droplets/precipitates where q*R >> 1.
    """
    if A_porod <= 0:
        return np.zeros_like(q)
    return A_porod / (np.asarray(q) ** float(exp_porod))


def calc_model_intensity(params, q, include_large_spheres=True, include_oz=True, integration_steps=300):
    """
    Calculates total scattering intensity from lmfit Parameters dictionary.
    Decomposes into primary fractal, large spheres, Porod tail, high-q component (PNC, FB, or OZ), and flat background.
    """
    rho2 = params['rho2'].value if hasattr(params['rho2'], 'value') else params['rho2']
    phi = params['phi'].value if hasattr(params['phi'], 'value') else params['phi']
    Rmean = params['Rmean'].value if hasattr(params['Rmean'], 'value') else params['Rmean']
    stdDev = params['stdDev'].value if hasattr(params['stdDev'], 'value') else params['stdDev']
    D = params['D'].value if hasattr(params['D'], 'value') else params['D']
    xi = params['xi'].value if hasattr(params['xi'], 'value') else params['xi']
    bkg = params['bkg'].value if hasattr(params['bkg'], 'value') else params['bkg']

    I_total = calc_primary_fractal_intensity(q, Rmean, stdDev, phi, D, xi, rho2=rho2, integration_steps=integration_steps)

    if include_large_spheres and 'total_volume_fraction_large' in params:
        phi_large = params['total_volume_fraction_large'].value if hasattr(params['total_volume_fraction_large'], 'value') else params['total_volume_fraction_large']
        if phi_large > 0:
            R_large = params['mean_radius_large'].value if hasattr(params['mean_radius_large'], 'value') else params['mean_radius_large']
            std_large = params['std_radius_large'].value if hasattr(params['std_radius_large'], 'value') else params['std_radius_large']
            I_large = calc_large_sphere_intensity(q, R_large, std_large, phi_large, rho2=rho2, integration_steps=integration_steps)
            I_total += I_large

    # Macroscopic Porod surface scattering component (for giant objects / precipitates with Guinier knee < q_min)
    if 'A_porod' in params:
        A_porod = params['A_porod'].value if hasattr(params['A_porod'], 'value') else params['A_porod']
        if A_porod > 0:
            exp_p = params['exp_porod'].value if ('exp_porod' in params and hasattr(params['exp_porod'], 'value')) else params.get('exp_porod', 4.0)
            I_total += calc_porod_intensity(q, A_porod, exp_p)

    # High-Q option 1: PNCs (small polydisperse spheres)
    if 'phi_pnc' in params and 'R_pnc' in params:
        phi_pnc = params['phi_pnc'].value if hasattr(params['phi_pnc'], 'value') else params['phi_pnc']
        if phi_pnc > 0:
            R_pnc = params['R_pnc'].value if hasattr(params['R_pnc'], 'value') else params['R_pnc']
            std_pnc = params['std_pnc'].value if hasattr(params['std_pnc'], 'value') else params['std_pnc']
            I_pnc = calc_pnc_intensity(q, R_pnc, std_pnc, phi_pnc, rho2=rho2, integration_steps=120)
            I_total += I_pnc

    # High-Q option 2: Viscoelastic network fluctuations (Fisher-Burford)
    elif 'A_fluc' in params and 'xi_fluc' in params:
        A_fluc = params['A_fluc'].value if hasattr(params['A_fluc'], 'value') else params['A_fluc']
        if A_fluc > 0:
            xi_fluc = params['xi_fluc'].value if hasattr(params['xi_fluc'], 'value') else params['xi_fluc']
            mu_fluc = params['mu_fluc'].value if hasattr(params['mu_fluc'], 'value') else (params.get('mu_fluc', 4.0))
            I_fb = calc_fb_intensity(q, A_fluc, xi_fluc, mu_fluc=mu_fluc)
            I_total += I_fb

    # High-Q option 3: Ornstein-Zernike fluctuations
    elif include_oz and 'A_OZ' in params and 'xi_OZ' in params:
        A_OZ = params['A_OZ'].value if hasattr(params['A_OZ'], 'value') else params['A_OZ']
        xi_OZ = params['xi_OZ'].value if hasattr(params['xi_OZ'], 'value') else params['xi_OZ']
        if A_OZ > 0:
            I_oz = calc_oz_intensity(q, A_OZ, xi_OZ)
            I_total += I_oz

    return I_total + bkg


def decompose_model_components(params, q, integration_steps=300):
    """
    Returns individual component profiles: (total, primary_fractal, large_spheres, high_q_term, background)
    """
    rho2 = params['rho2'].value if hasattr(params['rho2'], 'value') else params['rho2']
    phi = params['phi'].value if hasattr(params['phi'], 'value') else params['phi']
    Rmean = params['Rmean'].value if hasattr(params['Rmean'], 'value') else params['Rmean']
    stdDev = params['stdDev'].value if hasattr(params['stdDev'], 'value') else params['stdDev']
    D = params['D'].value if hasattr(params['D'], 'value') else params['D']
    xi = params['xi'].value if hasattr(params['xi'], 'value') else params['xi']
    bkg = params['bkg'].value if hasattr(params['bkg'], 'value') else params['bkg']

    I_primary = calc_primary_fractal_intensity(q, Rmean, stdDev, phi, D, xi, rho2=rho2, integration_steps=integration_steps)

    I_large = np.zeros_like(q)
    if 'total_volume_fraction_large' in params:
        phi_large = params['total_volume_fraction_large'].value if hasattr(params['total_volume_fraction_large'], 'value') else params['total_volume_fraction_large']
        if phi_large > 0:
            R_large = params['mean_radius_large'].value if hasattr(params['mean_radius_large'], 'value') else params['mean_radius_large']
            std_large = params['std_radius_large'].value if hasattr(params['std_radius_large'], 'value') else params['std_radius_large']
            I_large = calc_large_sphere_intensity(q, R_large, std_large, phi_large, rho2=rho2, integration_steps=integration_steps)

    if 'A_porod' in params:
        A_porod = params['A_porod'].value if hasattr(params['A_porod'], 'value') else params['A_porod']
        if A_porod > 0:
            exp_p = params['exp_porod'].value if ('exp_porod' in params and hasattr(params['exp_porod'], 'value')) else params.get('exp_porod', 4.0)
            I_porod = calc_porod_intensity(q, A_porod, exp_p)
            I_large = I_large + I_porod

    I_high = np.zeros_like(q)
    if 'phi_pnc' in params and 'R_pnc' in params:
        phi_pnc = params['phi_pnc'].value if hasattr(params['phi_pnc'], 'value') else params['phi_pnc']
        if phi_pnc > 0:
            R_pnc = params['R_pnc'].value if hasattr(params['R_pnc'], 'value') else params['R_pnc']
            std_pnc = params['std_pnc'].value if hasattr(params['std_pnc'], 'value') else params['std_pnc']
            I_high = calc_pnc_intensity(q, R_pnc, std_pnc, phi_pnc, rho2=rho2, integration_steps=120)
    elif 'A_fluc' in params and 'xi_fluc' in params:
        A_fluc = params['A_fluc'].value if hasattr(params['A_fluc'], 'value') else params['A_fluc']
        if A_fluc > 0:
            xi_fluc = params['xi_fluc'].value if hasattr(params['xi_fluc'], 'value') else params['xi_fluc']
            mu_val = params.get('mu_fluc', 4.0)
            mu_fluc = mu_val.value if hasattr(mu_val, 'value') else mu_val
            I_high = calc_fb_intensity(q, A_fluc, xi_fluc, mu_fluc=mu_fluc)
    elif 'A_OZ' in params and 'xi_OZ' in params:
        A_OZ = params['A_OZ'].value if hasattr(params['A_OZ'], 'value') else params['A_OZ']
        xi_OZ = params['xi_OZ'].value if hasattr(params['xi_OZ'], 'value') else params['xi_OZ']
        if A_OZ > 0:
            I_high = calc_oz_intensity(q, A_OZ, xi_OZ)

    I_total = I_primary + I_large + I_high + bkg
    return I_total, I_primary, I_large, I_high, np.full_like(q, bkg)


def calc_sharp_sq_vectorized(q, r0, D, Rc_grid, weights, n_u=250):
    """
    Vectorized structure factor for sharp-cutoff mass fractal with Gaussian polydispersity on Rc:
    S(q, Rc) = 1 + D * (Rc/r0)^D * (1/(q*Rc)) * int_0^1 u^(D-2) * sin(q*Rc*u) du
    (Freltoft, Kjems & Sinha 1986; Sorensen 2001)
    """
    u = np.linspace(1e-5, 1.0, n_u)
    du = u[1] - u[0]
    u_pow = u**(D - 2.0)
    
    x = q[:, np.newaxis] * Rc_grid[np.newaxis, :]
    xu = x[:, :, np.newaxis] * u[np.newaxis, np.newaxis, :]
    sin_xu = np.sin(xu)
    
    integral = np.sum(sin_xu * u_pow[np.newaxis, np.newaxis, :] * du, axis=2)
    J = np.where(x == 0, 1.0 / D, integral / x)
    
    N_Rc = D * (Rc_grid / r0)**D
    Sq_matrix = 1.0 + N_Rc[np.newaxis, :] * J
    Sq_poly = np.sum(Sq_matrix * weights[np.newaxis, :], axis=1)
    return np.maximum(Sq_poly, 1.0)



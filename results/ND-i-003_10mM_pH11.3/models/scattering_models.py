#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Core SAXS Scattering Models for ND-i-003 (Low Supersaturation Regime: 10 mM CaCO3, pH 11.3).
-----------------------------------------------------------------------------------------
Physical Model Architecture:
1. Low Supersaturation Physics:
   - Resides below the spinodal boundary; barrierless spinodal demixing is suppressed.
   - No continuous bicontinuous mass fractal network forms (Dm ~ 1.84 is absent).
   - Precipitation occurs via classical nucleation & growth of discrete amorphous calcium carbonate (ACC) nanospheres.
2. Short-Distance SAXS (3.0 m, q in [0.022, 2.54] nm^-1):
   - Purely flat solvent background. Zero detectable primary building blocks or high-q PNC fluctuations.
3. Long-Distance SAXS (31.0 m, q in [0.002, 0.25] nm^-1):
   - Governed by the form factor P(q) of growing, polydisperse spherical ACC nanoparticles (lognormal distribution).
   - Features a flat Guinier plateau (q < 0.01 nm^-1), form factor knee (q ~ 0.015 - 0.025 nm^-1),
     and textbook Porod power-law decay (I ~ q^-4), proving compact 3D objects with smooth interfaces (D = 3).
"""

import numpy as np
from scipy.special import gammaln

# Calibrated contrast factor for amorphous calcium carbonate in water
# (SLD_ACC - SLD_water)^2 = 5.185^2 = 26.884225
RHO2_CALIBRATED = 26.884225


def sphere_form_factor(q, radius):
    """
    Normalized form factor P(q, R) for a hard sphere of radius R.
    P(0) = 1.0. Vectorized across q and radius dimensions.
    """
    qr = q * radius
    with np.errstate(divide='ignore', invalid='ignore'):
        P = (3.0 * (np.sin(qr) - qr * np.cos(qr)) / (qr**3))**2
    return np.where(qr < 1e-5, 1.0, P)


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
    Calculates unit volume of the polydisperse sphere population:
    <V> = int (4/3)*pi*r^3 * f(r) dr
    """
    r_max = mean_r + 10.0 * std_r
    r_grid = np.linspace(1e-5, r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]
    integrand = (4.0 / 3.0) * np.pi * (r_grid**3) * log_normal_pdf(r_grid, mu, sigma)
    vol = np.sum(integrand * dr)
    return max(vol, 1e-12)


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


def calc_sphere_intensity(q, Rmean, stdDev, phi, bkg=0.0, rho2=1.0, dist_type='lognormal', integration_steps=300):
    """
    Vectorized integrated scattering intensity for polydisperse discrete spheres:
    I(q) = rho2 * int [ r^6 * n(r) * P(q, r) ] dr + bkg
    Supports dist_type='lognormal' and dist_type='schulz' (or 'schulz_flory').
    """
    if phi <= 0 or Rmean <= 0:
        return np.full_like(q, max(0.0, bkg))

    r_max = Rmean + 8.0 * stdDev
    r_grid = np.linspace(1e-4, r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]

    dist_key = str(dist_type).lower()
    if 'schulz' in dist_key or 'flory' in dist_key:
        pdf = schulz_flory_pdf(r_grid, Rmean, stdDev)
    else:
        mu, sigma = convert_mean_std_to_mu_sigma(Rmean, stdDev)
        pdf = log_normal_pdf(r_grid, mu, sigma)

    vol_per_unit = np.sum((4.0 / 3.0) * np.pi * (r_grid**3) * pdf * dr)
    number_density = phi / max(vol_per_unit, 1e-12)

    # Reshape for 2D broadcasting: r along rows (axis 0), q along cols (axis 1)
    r_col = r_grid[:, np.newaxis]
    q_row = np.asarray(q)[np.newaxis, :]

    P_qr = sphere_form_factor(q_row, r_col)
    dist = (pdf * number_density)[:, np.newaxis]

    integrand = (r_col**6) * dist * P_qr * dr
    np.nan_to_num(integrand, copy=False, nan=0.0)

    return rho2 * np.sum(integrand, axis=0) + bkg


def calc_large_sphere_intensity(q, mean_radius_large, std_radius_large, phi_large, rho2=1.0, integration_steps=200):
    """
    Vectorized scattering intensity for large macroscopic droplets (radius > 1000 nm):
    I_large(q) = rho2 * int [ r^6 * n(r) * P(q, r) ] dr
    Accurately captures the forward scattering rollover at q < 0.005 nm^-1.
    """
    if phi_large <= 0 or mean_radius_large <= 0:
        return np.zeros_like(q)
    mu, sigma = convert_mean_std_to_mu_sigma(mean_radius_large, std_radius_large)
    vol_per_unit = calc_total_volume_per_unit(mu, sigma, mean_radius_large, std_radius_large, integration_steps)
    number_density = phi_large / max(vol_per_unit, 1e-12)

    r_max = mean_radius_large + 6.0 * std_radius_large
    r_grid = np.linspace(max(1.0, mean_radius_large - 4.0 * std_radius_large), r_max, integration_steps, endpoint=True)
    dr = r_grid[1] - r_grid[0]

    r_col = r_grid[:, np.newaxis]
    q_row = np.asarray(q)[np.newaxis, :]

    P_qr = sphere_form_factor(q_row, r_col)
    dist = (log_normal_pdf(r_col, mu, sigma) * number_density)

    integrand = (r_col**6) * dist * P_qr * dr
    np.nan_to_num(integrand, copy=False, nan=0.0)
    return rho2 * np.sum(integrand, axis=0)


def calc_porod_intensity(q, A_porod, exp_porod=4.0):
    """
    Macroscopic Porod surface scattering from large objects / interfaces:
    I_porod(q) = A_porod / (q ** exp_porod)
    """
    if A_porod <= 0:
        return np.zeros_like(q)
    return float(A_porod) / (np.asarray(q, dtype=float) ** float(exp_porod))


def calc_short_distance_model(params, q):
    """
    Scattering model for short-distance data in ND-i-003.
    Since no primary nodes or fractal mesh exist at 10 mM, this returns flat solvent background.
    """
    bkg = params['bkg'].value if 'bkg' in params else 0.0
    return np.full_like(q, bkg)


def calc_equivalent_porod_amplitude(R_large, phi_large):
    """
    Computes asymptotic Porod amplitude A_porod = 27 * phi_large / (8 * pi * R_large)
    corresponding to macroscopic droplets of radius R_large and volume fraction phi_large.
    """
    if R_large <= 0 or phi_large <= 0:
        return 0.0
    return float(27.0 * phi_large / (8.0 * np.pi * R_large))


def calc_long_distance_model(params, q, dist_type=None):
    """
    Scattering model for long-distance data in ND-i-003:
    Discrete polydisperse spherical ACC nanoparticles + macroscopic droplets (or Porod interface) + flat background.
    """
    def _get_val(key, default=0.0):
        if key not in params:
            return default
        val = params[key]
        return val.value if hasattr(val, 'value') else val

    Rmean = _get_val('Rmean', 0.0)
    stdDev = _get_val('stdDev', 0.0)
    phi = _get_val('phi', 0.0)
    bkg = _get_val('bkg', 0.0)
    rho2 = _get_val('rho2', 1.0)

    if dist_type is None:
        dist_type = _get_val('dist_type', 'lognormal')

    I_tot = calc_sphere_intensity(q, Rmean, stdDev, phi, bkg=0.0, rho2=rho2, dist_type=dist_type)

    # Macroscopic droplet form factor
    if 'mean_radius_large' in params and _get_val('total_volume_fraction_large', 0.0) > 0:
        R_large = _get_val('mean_radius_large', 0.0)
        std_large = _get_val('std_radius_large', 250.0)
        phi_large = _get_val('total_volume_fraction_large', 0.0)
        I_tot = I_tot + calc_large_sphere_intensity(q, R_large, std_large, phi_large, rho2=rho2)
    elif 'R_large' in params and _get_val('phi_large', 0.0) > 0:
        R_large = _get_val('R_large', 0.0)
        std_large = _get_val('std_large', 250.0)
        phi_large = _get_val('phi_large', 0.0)
        I_tot = I_tot + calc_large_sphere_intensity(q, R_large, std_large, phi_large, rho2=rho2)
    elif 'A_porod' in params and _get_val('A_porod', 0.0) > 0:
        A_p = _get_val('A_porod', 0.0)
        exp_p = _get_val('exp_porod', 4.0)
        I_tot = I_tot + calc_porod_intensity(q, A_p, exp_p)

    return I_tot + bkg


def calc_composite_model(params, q, dist_type=None):
    """
    Scattering model for full-range composite curve:
    - Primary Mature ACC Nanospheres: Rmean (or R1), stdDev (or std1), phi (or phi1)
    - Growing / Secondary ACC Nanospheres (if present): R2, stdDev2, phi2
    - Macroscopic Droplets (if present): mean_radius_large, std_radius_large, total_volume_fraction_large
    - Porod Interface (if present): A_porod, exp_porod
    - Incoherent Background: bkg
    """
    def _get_val(key, default=0.0):
        if key not in params:
            return default
        val = params[key]
        return val.value if hasattr(val, 'value') else val

    rho2 = _get_val('rho2', 1.0)
    bkg = _get_val('bkg', 0.0)

    # Primary / Mature ACC population
    R1 = _get_val('R1', _get_val('Rmean', 0.0))
    std1 = _get_val('std1', _get_val('stdDev', 0.0))
    phi1 = _get_val('phi1', _get_val('phi', 0.0))
    I_tot = calc_sphere_intensity(q, R1, std1, phi1, bkg=0.0, rho2=rho2, dist_type='lognormal')

    # Secondary / Growing ACC population
    if 'R2' in params and _get_val('phi2', 0.0) > 0:
        R2 = _get_val('R2', 0.0)
        std2 = _get_val('stdDev2', _get_val('std2', 0.0))
        phi2 = _get_val('phi2', 0.0)
        I_tot = I_tot + calc_sphere_intensity(q, R2, std2, phi2, bkg=0.0, rho2=rho2, dist_type='lognormal')

    # Macroscopic droplet population
    if 'mean_radius_large' in params and _get_val('total_volume_fraction_large', 0.0) > 0:
        R_large = _get_val('mean_radius_large', 0.0)
        std_large = _get_val('std_radius_large', 0.0)
        phi_large = _get_val('total_volume_fraction_large', 0.0)
        I_tot = I_tot + calc_large_sphere_intensity(q, R_large, std_large, phi_large, rho2=rho2)

    # Porod interface
    if 'A_porod' in params and _get_val('A_porod', 0.0) > 0:
        A_p = _get_val('A_porod', 0.0)
        exp_p = _get_val('exp_porod', 4.0)
        I_tot = I_tot + calc_porod_intensity(q, A_p, exp_p)

    return I_tot + bkg


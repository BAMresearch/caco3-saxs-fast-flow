#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Fitting Engine and Parameter Builders for ND-i-003.
---------------------------------------------------
Implements:
1. Standardized lmfit Parameters builders for short-distance, long-distance, and composite models.
2. Sequential continuity regularization ensuring smooth frame-to-frame drift while preventing unphysical jumping.
"""

from lmfit import Parameters, minimize
import numpy as np


def create_short_distance_parameters(bkg_init=1.0e-4):
    """ Builds parameter set for short-distance background model. """
    p = Parameters()
    p.add('bkg', value=bkg_init, min=-0.05, max=0.05, vary=True)
    return p


def create_induction_parameters(
    R_large_init=2600.0, std_large_init=250.0, phi_large_init=5.0e-5, bkg_init=0.01
):
    """ Builds parameter set for induction frames (macroscopic droplets + flat background). """
    p = Parameters()
    p.add('rho2', value=1.0, vary=False)
    p.add('mean_radius_large', value=R_large_init, min=1500.0, max=4000.0, vary=True)
    p.add('std_radius_large', value=std_large_init, min=50.0, max=600.0, vary=False)
    p.add('total_volume_fraction_large', value=phi_large_init, min=1.0e-8, max=1.0e-3, vary=True)
    p.add('bkg', value=bkg_init, min=-0.05, max=0.25, vary=True)
    return p


def create_long_distance_parameters(
    Rmean_init=120.0, stdDev_init=18.0, phi_init=4.2e-4,
    R_large_init=2600.0, std_large_init=250.0, phi_large_init=5.5e-5,
    bkg_init=0.01, dist_type='lognormal'
):
    """ Builds parameter set for long-distance polydisperse ACC sphere model + macroscopic droplets. """
    p = Parameters()
    p.add('rho2', value=1.0, vary=False)
    p.add('Rmean', value=Rmean_init, min=15.0, max=160.0, vary=True)
    p.add('stdDev', value=stdDev_init, min=1.0, max=25.0, vary=True)
    p.add('phi', value=phi_init, min=1.0e-12, max=0.01, vary=True)
    p.add('mean_radius_large', value=R_large_init, min=1500.0, max=4000.0, vary=True)
    p.add('std_radius_large', value=std_large_init, min=50.0, max=600.0, vary=False)
    p.add('total_volume_fraction_large', value=phi_large_init, min=1.0e-8, max=1.0e-3, vary=True)
    p.add('bkg', value=bkg_init, min=-0.01, max=0.25, vary=True)
    p.dist_type = dist_type
    return p


def create_composite_parameters(
    R1_init=160.0, std1_init=14.0, phi1_init=2.0e-4,
    R2_init=105.0, std2_init=22.0, phi2_init=4.0e-4,
    R_large_init=2700.0, std_large_init=300.0, phi_large_init=8.0e-5,
    bkg_init=1.0e-4
):
    """
    Builds parameter set for composite anchor curve:
    - Mature ACC Nanospheres: Rmean (R1), stdDev (std1), phi (phi1)
    - Growing ACC Nanospheres: R2, stdDev2, phi2
    - Macroscopic Droplets: mean_radius_large, std_radius_large, total_volume_fraction_large
    - Incoherent Background: bkg
    """
    p = Parameters()
    p.add('rho2', value=1.0, vary=False)

    # Mature ACC Nanospheres (narrow distribution generating distinct Bessel oscillations)
    p.add('Rmean', value=R1_init, min=130.0, max=180.0, vary=True)
    p.add('stdDev', value=std1_init, min=5.0, max=18.0, vary=True)
    p.add('phi', value=phi1_init, min=1.0e-7, max=1.0e-3, vary=True)

    # Growing / Secondary ACC Nanospheres (intermediate size filling oscillation trough)
    p.add('R2', value=R2_init, min=50.0, max=130.0, vary=True)
    p.add('stdDev2', value=std2_init, min=5.0, max=35.0, vary=True)
    p.add('phi2', value=phi2_init, min=1.0e-7, max=1.0e-3, vary=True)

    # Macroscopic Droplets (capturing ultra-low q forward scattering rollover)
    p.add('mean_radius_large', value=R_large_init, min=1200.0, max=5000.0, vary=True)
    p.add('std_radius_large', value=std_large_init, min=50.0, max=800.0, vary=True)
    p.add('total_volume_fraction_large', value=phi_large_init, min=1.0e-7, max=1.0e-3, vary=True)

    # Flat solvent background
    p.add('bkg', value=bkg_init, min=0.0, max=0.01, vary=True)
    return p


class SequentialContinuityFitter:
    """
    Handles frame-to-frame evolutionary regularization.
    """

    def __init__(self, model_func, reg_weights=None):
        self.model_func = model_func
        self.reg_weights = reg_weights or {'Rmean': 0.05, 'stdDev': 0.05, 'phi': 0.02}

    def residual(self, params, q, y, dy, prev_params=None):
        ymod = self.model_func(params, q)
        res = (y - ymod) / dy

        if prev_params is not None:
            reg_terms = []
            for name, weight in self.reg_weights.items():
                if name in params and name in prev_params and params[name].vary:
                    curr_val = params[name].value
                    prev_val = prev_params[name].value
                    scale = max(abs(prev_val), 1e-4)
                    reg_terms.append(weight * (curr_val - prev_val) / scale)
            if reg_terms:
                res = np.concatenate([res, np.array(reg_terms)])

        return res

    def fit_frame(self, params, q, y, dy, prev_params=None, method='least_squares'):
        out = minimize(self.residual, params, args=(q, y, dy, prev_params), method=method)
        return out

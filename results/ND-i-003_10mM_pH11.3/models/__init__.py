#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Models module for ND-i-003 unified SAXS analysis suite.
Discrete Amorphous Calcium Carbonate (ACC) Nanosphere Form Factor & Kinetics.
"""

from .scattering_models import (
    RHO2_CALIBRATED,
    sphere_form_factor,
    log_normal_pdf,
    convert_mean_std_to_mu_sigma,
    calc_total_volume_per_unit,
    calc_sphere_intensity,
    calc_large_sphere_intensity,
    calc_porod_intensity,
    calc_equivalent_porod_amplitude,
    calc_short_distance_model,
    calc_long_distance_model,
    calc_composite_model
)

from .fitting_engine import (
    SequentialContinuityFitter,
    create_short_distance_parameters,
    create_long_distance_parameters,
    create_induction_parameters,
    create_composite_parameters
)

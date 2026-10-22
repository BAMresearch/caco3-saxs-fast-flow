#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Standardised SAXS Plotting Palette and Technical Style Definitions
------------------------------------------------------------------
Provides invariant colour mappings and technical line styles for SAXS/WAXS
scattering models and component decompositions across all reaction series.

Standardised physical components:
1. Total Unified Fit: #dc2626 (crimson red solid line)
2. Primary Spinodal Fractal Network: #059669 (emerald green dashed line)
3. Discrete ACC Spheres (Mode 1): #0284c7 (ocean blue dash-dot line)
4. Discrete ACC Spheres (Mode 2): #0ea5e9 (sky blue dotted line)
5. Tanaka Viscoelastic PNC Fluctuation Mesh: #9333ea (violet purple dotted line)
6. Macroscopic Droplets / Porod Tail: #d97706 (amber gold dash-dot line)
7. Solvent Background: #64748b (slate grey dotted line)
8. Experimental Scattering Data: #334155 (slate charcoal scatter markers)
9. Residuals: #475569 (muted slate markers)
"""

SAXS_COLORS = {
    # Core physical scattering components
    'data': '#334155',             # Slate charcoal
    'total_fit': '#dc2626',        # Crimson red
    'fractal_network': '#059669',  # Emerald green
    'acc_spheres': '#0284c7',      # Ocean blue (primary mode)
    'acc_spheres_mode2': '#0ea5e9',# Sky blue (secondary mode)
    'pnc_fluctuations': '#9333ea', # Violet purple (Tanaka viscoelastic mesh / PNC fluctuations)
    'droplets_porod': '#d97706',   # Amber gold (macroscopic droplets / Porod tail)
    'background': '#64748b',       # Slate grey (incoherent solvent baseline)
    'residuals': '#475569',        # Muted slate

    # Concentration regime identifiers
    'regime_10mM': '#38bdf8',      # Light blue (sub-spinodal)
    'regime_15mM': '#2dd4bf',      # Teal (crossover boundary)
    'regime_20mM': '#4ade80',      # Light green (spinodal gelation)
    'regime_40mM': '#f87171',      # Coral red (benchmark spinodal)
    'regime_additive': '#c084fc',  # Light purple (additive perturbation branch)
}

SAXS_STYLES = {
    'data': {
        'color': SAXS_COLORS['data'],
        'marker': 'o',
        'markersize': 2.5,
        'alpha': 0.45,
        'linestyle': 'none',
        'label': 'Experimental Scattering Data'
    },
    'total_fit': {
        'color': SAXS_COLORS['total_fit'],
        'linestyle': '-',
        'linewidth': 2.2,
        'label': 'Total Unified Fit'
    },
    'fractal_network': {
        'color': SAXS_COLORS['fractal_network'],
        'linestyle': '--',
        'linewidth': 1.8,
        'label': 'Spinodal Fractal Network'
    },
    'acc_spheres': {
        'color': SAXS_COLORS['acc_spheres'],
        'linestyle': '-.',
        'linewidth': 1.8,
        'label': 'Discrete ACC Spheres (Mode 1)'
    },
    'acc_spheres_mode2': {
        'color': SAXS_COLORS['acc_spheres_mode2'],
        'linestyle': ':',
        'linewidth': 1.6,
        'label': 'Discrete ACC Spheres (Mode 2)'
    },
    'pnc_fluctuations': {
        'color': SAXS_COLORS['pnc_fluctuations'],
        'linestyle': ':',
        'linewidth': 1.8,
        'label': 'Tanaka PNC Fluctuation Mesh'
    },
    'droplets_porod': {
        'color': SAXS_COLORS['droplets_porod'],
        'linestyle': '-.',
        'linewidth': 1.8,
        'label': 'Macroscopic Droplets / Porod Tail'
    },
    'background': {
        'color': SAXS_COLORS['background'],
        'linestyle': ':',
        'linewidth': 1.2,
        'label': 'Solvent Background'
    }
}


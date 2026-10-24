#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 05: Build Unified Multi-Scale Interactive Dashboard for ND-i-003.
-----------------------------------------------------------------------
Generates a standalone, publication-grade HTML dashboard with 6 standardized tabs:
1. Overview & Composite Curve: Full composite scattering curve (q in [0.002, 2.54] nm^-1) + baseline parameters
2. Cross-Scale Kinetics: Multi-scale kinetic synthesis (3.0 m & 31.0 m)
3. Short Distance (3.0 m): Interactive frame reviewer (pure flat solvent background verification)
4. Long Distance (31.0 m): Interactive frame reviewer (classical ACC nanosphere nucleation & growth)
5. Parameter Tables: Full searchable & scrollable parameter tables (+/- 1 sigma)
6. Physical Formalism: Mathematical equations and models rendered with MathJax

Output:
  output/dashboard/index.html
"""

import os
import sys
import json
import pandas as pd

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)


def build_dashboard():
    dashboard_dir = os.path.join(PROJECT_DIR, 'output', 'dashboard')
    os.makedirs(dashboard_dir, exist_ok=True)

    # 1. Load Short Distance Parameters
    p_short_csv = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
    df_short = pd.read_csv(p_short_csv)
    short_telemetry = {}
    for idx, r in df_short.iterrows():
        f = int(r['frame'])
        t_sec = float(r.get('time_seconds', 0.0))
        bkg = float(r['bkg'])
        err_bkg = float(r.get('std_err_bkg', 0.0))
        bkg_str = f"{bkg*1e4:.2f} ± {err_bkg*1e4:.2f}" if err_bkg > 0 else f"{bkg*1e4:.2f}"

        short_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': round(float(r['redchi']), 3),
            'bkg': bkg_str,
            'mean_I': f"{float(r.get('mean_I', 0.0))*1e4:.2f}",
            'pos_pct': round(float(r.get('pos_pct', 50.0)), 1),
            'pts': int(r.get('pts', 1638))
        }

    # 2. Load Long Distance Parameters
    p_long_csv = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')
    df_long = pd.read_csv(p_long_csv)
    long_telemetry = {}
    for idx, r in df_long.iterrows():
        f = int(r['frame'])
        t_sec = float(r.get('time_seconds', 0.0))
        is_ind = bool(r.get('is_induction', False))

        r_mean = float(r['Rmean']) if pd.notnull(r.get('Rmean')) else 0.0
        err_r = float(r.get('std_err_Rmean', 0.0)) if pd.notnull(r.get('std_err_Rmean')) else 0.0
        r_str = f"{r_mean:.1f} ± {err_r:.1f}" if err_r > 0 else f"{r_mean:.1f}"

        std_r = float(r['stdDev']) if pd.notnull(r.get('stdDev')) else 0.0
        err_std = float(r.get('std_err_stdDev', 0.0)) if pd.notnull(r.get('std_err_stdDev')) else 0.0
        std_str = f"{std_r:.1f} ± {err_std:.1f}" if err_std > 0 else f"{std_r:.1f}"

        phi = float(r['phi']) if pd.notnull(r.get('phi')) else 0.0
        err_phi = float(r.get('std_err_phi', 0.0)) if pd.notnull(r.get('std_err_phi')) else 0.0
        phi_str = f"{phi:.2e} ± {err_phi:.2e}" if err_phi > 0 else f"{phi:.2e}"

        r_large = float(r['mean_radius_large']) if pd.notnull(r.get('mean_radius_large')) else 0.0
        err_rlarge = float(r.get('std_err_mean_radius_large', 0.0)) if pd.notnull(r.get('std_err_mean_radius_large')) else 0.0
        rlarge_str = f"{r_large:.0f} ± {err_rlarge:.0f} nm" if err_rlarge > 0 else f"{r_large:.0f} nm"

        phi_large = float(r['total_volume_fraction_large']) if pd.notnull(r.get('total_volume_fraction_large')) else 0.0
        err_philarge = float(r.get('std_err_total_volume_fraction_large', 0.0)) if pd.notnull(r.get('std_err_total_volume_fraction_large')) else 0.0
        philarge_str = f"{phi_large*1e5:.2f} ± {err_philarge*1e5:.2f} × 10⁻⁵" if err_philarge > 0 else f"{phi_large*1e5:.2f} × 10⁻⁵"
        droplet_str = f"{r_large:.0f} nm, {phi_large*1e5:.2f} × 10⁻⁵"

        a_porod = float(r['A_porod']) if pd.notnull(r.get('A_porod')) else 0.0
        err_aporod = float(r.get('std_err_A_porod', 0.0)) if pd.notnull(r.get('std_err_A_porod')) else 0.0
        aporod_str = f"{a_porod*1e8:.2f} ± {err_aporod*1e8:.2f}" if err_aporod > 0 else f"{a_porod*1e8:.2f}"

        bkg = float(r['bkg']) if pd.notnull(r.get('bkg')) else 0.0
        err_bkg = float(r.get('std_err_bkg', 0.0)) if pd.notnull(r.get('std_err_bkg')) else 0.0
        bkg_str = f"{bkg:.3f} ± {err_bkg:.3f}" if err_bkg > 0 else f"{bkg:.3f}"

        long_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': round(float(r['redchi']), 3) if pd.notnull(r.get('redchi')) else '-',
            'Rmean': r_str if not is_ind else 'N/A (Induction)',
            'stdDev': std_str if not is_ind else 'N/A (Induction)',
            'poly_pct': round(std_r / max(r_mean, 1e-3) * 100, 1) if not is_ind else ':',
            'phi': phi_str if not is_ind else ':',
            'droplet': droplet_str,
            'A_porod': f"{aporod_str} × 10⁻⁸",
            'bkg': bkg_str,
            'is_ind': is_ind
        }

    # 3. Load Baseline Parameters
    bp_json = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
    with open(bp_json, 'r') as f:
        bp = json.load(f)

    # 4. Construct Full HTML
    raw_template = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ND-i-003 Multi-Scale SAXS Kinetics: Classical ACC Nanosphere Nucleation &amp; Growth</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <script>
    window.MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        displayMath: [['$$', '$$'], ['\\[', '\\]']],
        processEscapes: true
      },
      options: {
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
      }
    };
    </script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <style>
        :root {
            --bg-primary: #0b0f19;
            --bg-card: #151d2f;
            --bg-card-hover: #1c263c;
            --border-color: rgba(255, 255, 255, 0.08);
            --border-hover: rgba(56, 189, 248, 0.3);
            --text-primary: #f1f5f9;
            --text-secondary: #94a3b8;
            --text-muted: #64748b;
            --accent-blue: #38bdf8;
            --accent-teal: #2dd4bf;
            --accent-purple: #c084fc;
            --accent-amber: #fbbf24;
            --accent-green: #34d399;
            --accent-red: #f43f5e;
            --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            --font-mono: 'JetBrains Mono', monospace;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            background-color: var(--bg-primary);
            color: var(--text-primary);
            font-family: var(--font-sans);
            line-height: 1.5;
            -webkit-font-smoothing: antialiased;
            padding: 1.5rem;
        }

        .container { max-width: 1440px; margin: 0 auto; }

        header {
            margin-bottom: 1.5rem;
            padding-bottom: 1.5rem;
            border-bottom: 1px solid var(--border-color);
        }
        .badge-bar { display: flex; flex-wrap: wrap; gap: 0.5rem; margin-bottom: 0.8rem; }
        .badge {
            font-size: 0.75rem;
            font-weight: 600;
            padding: 0.25rem 0.6rem;
            border-radius: 9999px;
            display: inline-flex;
            align-items: center;
            border: 1px solid transparent;
        }
        .badge-blue { background: rgba(56, 189, 248, 0.12); color: var(--accent-blue); border-color: rgba(56, 189, 248, 0.25); }
        .badge-teal { background: rgba(45, 212, 191, 0.12); color: var(--accent-teal); border-color: rgba(45, 212, 191, 0.25); }
        .badge-purple { background: rgba(192, 132, 252, 0.12); color: var(--accent-purple); border-color: rgba(192, 132, 252, 0.25); }
        .badge-amber { background: rgba(251, 191, 36, 0.12); color: var(--accent-amber); border-color: rgba(251, 191, 36, 0.25); }

        h1 { font-size: 1.6rem; font-weight: 700; color: #fff; letter-spacing: -0.02em; margin-bottom: 0.4rem; }
        p.subtitle { color: var(--text-secondary); font-size: 0.95rem; max-width: 1100px; }

        .nav-tabs {
            display: flex; gap: 0.5rem; margin-bottom: 1.5rem; border-bottom: 1px solid var(--border-color);
            padding-bottom: 0.75rem; overflow-x: auto;
        }
        .nav-tab {
            background: transparent; border: 1px solid transparent; color: var(--text-secondary);
            padding: 0.6rem 1.0rem; border-radius: 8px; font-weight: 600; font-size: 0.88rem;
            cursor: pointer; transition: all 0.2s ease; display: inline-flex; align-items: center; gap: 0.4rem;
            white-space: nowrap;
        }
        .nav-tab:hover { color: var(--text-primary); background: rgba(255, 255, 255, 0.04); }
        .nav-tab.active {
            color: #fff; background: var(--bg-card); border-color: var(--border-hover);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
        }

        .section-pane { display: none; }
        .section-pane.active { display: block; animation: fadeIn 0.25s ease-in-out; }
        @keyframes fadeIn { from { opacity: 0; transform: translateY(4px); } to { opacity: 1; transform: translateY(0); } }

        .card {
            background: var(--bg-card); border: 1px solid var(--border-color);
            border-radius: 12px; padding: 1.5rem; margin-bottom: 1.5rem;
        }
        .card-title { font-size: 1.15rem; font-weight: 700; color: #fff; margin-bottom: 0.8rem; }

        .reviewer-layout {
            display: grid; grid-template-columns: 1fr 340px; gap: 1.5rem; align-items: start;
        }
        @media (max-width: 950px) { .reviewer-layout { grid-template-columns: 1fr; } }

        .plot-container {
            background: #090d16; border: 1px solid var(--border-color);
            border-radius: 8px; padding: 0.5rem; display: flex; justify-content: center;
            align-items: center; min-height: 480px;
        }
        .plot-container img { width: 100%; height: auto; max-height: 520px; object-fit: contain; border-radius: 4px; }

        .controls-panel { display: flex; flex-direction: column; gap: 1.2rem; }
        .slider-box {
            background: rgba(11, 15, 25, 0.6); border: 1px solid var(--border-color);
            border-radius: 8px; padding: 1.2rem;
        }
        .slider-label {
            display: flex; justify-content: space-between; font-size: 0.85rem;
            font-weight: 600; margin-bottom: 0.6rem;
        }
        input[type=range] {
            width: 100%; height: 6px; background: #24344d; border-radius: 3px;
            outline: none; cursor: pointer; accent-color: var(--accent-blue);
        }
        .playback-bar { display: flex; gap: 0.5rem; margin-top: 0.8rem; }
        .btn {
            background: #1e293b; border: 1px solid var(--border-color); color: var(--text-primary);
            padding: 0.4rem 0.75rem; border-radius: 6px; font-size: 0.82rem; font-weight: 600;
            cursor: pointer; transition: all 0.2s ease; flex: 1; text-align: center;
        }
        .btn:hover { background: rgba(255, 255, 255, 0.08); border-color: var(--accent-blue); }
        .btn-play { background: rgba(56, 189, 248, 0.15); color: var(--accent-blue); border-color: rgba(56, 189, 248, 0.4); }
        .btn-pause { background: rgba(244, 63, 94, 0.15); color: var(--accent-red); border-color: rgba(244, 63, 94, 0.4); }

        .telemetry-card {
            background: rgba(11, 15, 25, 0.6); border: 1px solid var(--border-color);
            border-radius: 8px; padding: 1.0rem; font-family: 'JetBrains Mono', monospace; font-size: 0.82rem;
        }
        .tele-row {
            display: flex; justify-content: space-between; padding: 0.35rem 0;
            border-bottom: 1px dashed rgba(255, 255, 255, 0.06);
        }
        .tele-row:last-child { border-bottom: none; }
        .tele-label { color: var(--text-secondary); }
        .tele-val { color: var(--text-primary); font-weight: 600; }

        .hero-img { width: 100%; border-radius: 8px; border: 1px solid var(--border-color); margin-top: 0.8rem; }
        .equation-box {
            background: #090d16; border: 1px solid var(--border-color); border-left: 4px solid var(--accent-blue);
            border-radius: 6px; padding: 1.2rem; margin: 1.0rem 0; overflow-x: auto;
        }
        .data-table-container {
            max-height: 480px; overflow-y: auto; overflow-x: auto; border: 1px solid var(--border-color);
            border-radius: 8px; margin-top: 1rem;
        }
        .data-table {
            width: 100%; border-collapse: collapse; font-family: var(--font-mono); font-size: 0.82rem;
            text-align: right;
        }
        .data-table th, .data-table td {
            padding: 0.5rem 0.75rem; border-bottom: 1px solid rgba(255, 255, 255, 0.05);
        }
        .data-table th {
            position: sticky; top: 0; background: #1e293b; color: var(--accent-blue);
            font-weight: 600; z-index: 10;
        }
        .data-table tr:hover td { background: rgba(255, 255, 255, 0.03); }
        .telemetry-table { width: 100%; border-collapse: collapse; }
        .telemetry-table th, .telemetry-table td {
            padding: 0.6rem 0.8rem; border-bottom: 1px solid rgba(255, 255, 255, 0.06);
            text-align: right; font-family: var(--font-mono); font-size: 0.85rem;
        }
        .telemetry-table th:first-child, .telemetry-table td:first-child { text-align: left; }
    </style>
</head>
<body>
<div class="container">
    <header>
        <div class="badge-bar">
            <a href="../../../index.html" class="badge badge-blue" style="text-decoration: none; cursor: pointer; border-color: var(--accent-blue);">&larr; Master Dashboard</a>
            <span class="badge badge-blue">Reaction Series ND-i-003</span>
            <span class="badge badge-teal">Low Supersaturation: 10 mM CaCO3, pH 11.3</span>
            <span class="badge badge-purple">Distances: 3.0 m &amp; 31.0 m (ESRF ID02)</span>
            <span class="badge badge-amber">Classical ACC Nucleation &amp; Growth</span>
        </div>
        <h1>CaCO<sub>3</sub> Multi-Scale SAXS Kinetics: Discrete ACC Nanosphere Formation</h1>
        <p class="subtitle">
            Physical investigation of the low supersaturation regime (10 mM pure CaCO<sub>3</sub>, pH 11.3) below the spinodal decomposition threshold.
            At 3.0 m, the short-distance signal verifies flat solvent background without high-q prenucleation clusters or fractal gel nodes.
            At 31.0 m, the long-distance signal demonstrates the emergence and growth of discrete polydisperse amorphous calcium carbonate (ACC)
            nanospheres with textbook Porod decay ($I \propto q^{-4}$).
        </p>
    </header>

    <!-- Navigation Tabs -->
    <div class="nav-tabs">
        <button class="nav-tab active" onclick="switchSection('overview')">1. Overview &amp; Composite Curve</button>
        <button class="nav-tab" onclick="switchSection('trajectories')">2. Cross-Scale Kinetics</button>
        <button class="nav-tab" onclick="switchSection('short')">3. Short Distance (3.0 m)</button>
        <button class="nav-tab" onclick="switchSection('long')">4. Long Distance (31.0 m)</button>
        <button class="nav-tab" onclick="switchSection('tables')">5. Parameter Tables</button>
        <button class="nav-tab" onclick="switchSection('theory')">6. Physical Formalism</button>
    </div>

    <!-- ==================== TAB 1: OVERVIEW & COMPOSITE CURVE ==================== -->
    <div id="pane-overview" class="section-pane active">
        <div class="card">
            <div class="card-title">Full-Range Composite Baseline Anchor Curve ($q \in [0.002, 2.54]\text{ nm}^{-1}$)</div>
            <p style="color: var(--text-secondary); font-size: 0.9rem; margin-bottom: 1.0rem;">
                Combined composite scattering curve spanning over 3 decades in $q$ and 6 orders of magnitude in intensity, fitted with bimodal ACC nanospheres
                (narrow mature population preserving Bessel oscillations + growing population matching the intermediate shoulder) + macroscopic droplet form factor
                (reproducing ultra-low-$q$ forward scattering rollover) + flat solvent background.
            </p>
            <div style="display: flex; justify-content: center;">
                <img style="max-width: 950px; width: 100%; border-radius: 8px; border: 1px solid var(--border-color);" src="../composite/composite_baseline_fit.png" alt="Composite Baseline Fit">
            </div>
        </div>

        <div class="card">
            <div class="card-title">Converged Structural Baseline Parameters (&plusmn;1&sigma;)</div>
            <table class="telemetry-table" style="font-size: 0.9rem;">
                <thead>
                    <tr style="border-bottom: 2px solid var(--border-color);">
                        <th style="color: var(--accent-blue);">Component</th>
                        <th style="color: var(--accent-blue);">Parameter</th>
                        <th style="color: var(--accent-blue);">Symbol</th>
                        <th style="color: var(--accent-blue); text-align: right;">Value</th>
                        <th style="color: var(--accent-blue); text-align: right;">&plusmn;1&sigma; Error</th>
                        <th style="color: var(--accent-blue); text-align: left; padding-left: 1.5rem;">Physical Interpretation</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td rowspan="3" style="font-weight: 600; color: var(--accent-blue); border-right: 1px solid var(--border-color);">Mature ACC Nanospheres</td>
                        <td>Mean Radius</td><td>&lang;R<sub>1</sub>&rang;</td>
                        <td>__BP_RMEAN__ nm</td><td>&plusmn; __BP_ERR_RMEAN__ nm</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mature nanospheres generating distinct Bessel form factor oscillations</td>
                    </tr>
                    <tr>
                        <td>Dispersion Width</td><td>&sigma;<sub>1</sub></td>
                        <td>__BP_STDDEV__ nm</td><td>&plusmn; __BP_ERR_STDDEV__ nm</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Narrow log-normal dispersion (__BP_POLY_PCT__% polydispersity)</td>
                    </tr>
                    <tr>
                        <td>Volume Fraction</td><td>&phi;<sub>1</sub></td>
                        <td>__BP_PHI__ &times; 10<sup>-4</sup></td><td>&plusmn; __BP_ERR_PHI__ &times; 10<sup>-4</sup></td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of mature ACC nanosphere population</td>
                    </tr>
                    <tr>
                        <td rowspan="3" style="font-weight: 600; color: var(--accent-amber); border-right: 1px solid var(--border-color);">Growing ACC Nanospheres</td>
                        <td>Mean Radius</td><td>&lang;R<sub>2</sub>&rang;</td>
                        <td>__BP_R2__ nm</td><td>&plusmn; __BP_ERR_R2__ nm</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Intermediate growing population consistent with late time-resolved series</td>
                    </tr>
                    <tr>
                        <td>Dispersion Width</td><td>&sigma;<sub>2</sub></td>
                        <td>__BP_STDDEV2__ nm</td><td>&plusmn; __BP_ERR_STDDEV2__ nm</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Broader log-normal dispersion (__BP_POLY2_PCT__% polydispersity)</td>
                    </tr>
                    <tr>
                        <td>Volume Fraction</td><td>&phi;<sub>2</sub></td>
                        <td>__BP_PHI2__ &times; 10<sup>-4</sup></td><td>&plusmn; __BP_ERR_PHI2__ &times; 10<sup>-4</sup></td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction filling in mid-q shoulder under Bessel troughs</td>
                    </tr>
                    <tr>
                        <td rowspan="3" style="font-weight: 600; color: var(--accent-purple); border-right: 1px solid var(--border-color);">Macroscopic Droplets</td>
                        <td>Droplet Radius</td><td>&lang;R<sub>large</sub>&rang;</td>
                        <td>__BP_RLARGE__ nm</td><td>&plusmn; __BP_ERR_RLARGE__ nm</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Capillary-pinched macroscopic droplets (replaces pure Porod)</td>
                    </tr>
                    <tr>
                        <td>Dispersion Width</td><td>&sigma;<sub>large</sub></td>
                        <td>__BP_STDLARGE__ nm</td><td>&plusmn; __BP_ERR_STDLARGE__ nm</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Droplet dispersion width (__BP_POLYLARGE_PCT__% polydispersity)</td>
                    </tr>
                    <tr>
                        <td>Volume Fraction</td><td>&phi;<sub>large</sub></td>
                        <td>__BP_PHILARGE__ &times; 10<sup>-5</sup></td><td>&plusmn; __BP_ERR_PHILARGE__ &times; 10<sup>-5</sup></td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Captures forward scattering rollover at q &lt; 0.005 nm<sup>-1</sup></td>
                    </tr>
                    <tr>
                        <td colspan="2" style="font-weight: 600;">Solvent Background</td><td>I<sub>bkg</sub></td>
                        <td>__BP_BKG__ &times; 10<sup>-5</sup> cm<sup>-1</sup></td><td>&plusmn; __BP_ERR_BKG__ &times; 10<sup>-5</sup></td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Flat incoherent solvent background</td>
                    </tr>
                    <tr>
                        <td colspan="2" style="font-weight: 600;">Goodness of Fit</td><td>&chi;<sup>2</sup><sub>red</sub></td>
                        <td style="color: var(--accent-green); font-weight: 700;">__BP_REDCHI__</td><td>&mdash;</td>
                        <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Full composite curve reduced chi-squared</td>
                    </tr>
                </tbody>
            </table>
        </div>
    </div>

    <!-- ==================== TAB 2: CROSS-SCALE KINETICS ==================== -->
    <div id="pane-trajectories" class="section-pane">
        <div class="card">
            <h2 class="card-title">Cross-Scale Kinetic Trajectories Synthesis ($t = 0.060\text{ s} \to 10.88\text{ s}$)</h2>
            <p style="color: var(--text-secondary); margin-bottom: 1rem; font-size: 0.9rem;">
                Cross-scale synthesis combining 3.0 m (WAXS/MAXS) and 31.0 m (SAXS) measurements with shaded &plusmn;1&sigma; uncertainty envelopes.
            </p>
            <img class="hero-img" src="../unified_trajectories/unified_cross_scale_trajectories.png" alt="Unified Multi-Scale Trajectories">

            <h3 style="font-size: 1.05rem; color: var(--accent-blue); margin-top: 1.8rem;">Long-Distance ACC Nanosphere Kinetics (31.0 m)</h3>
            <img class="hero-img" src="../long_distance/trajectories.png" alt="Long Distance ACC Nanosphere Trajectories">

            <h3 style="font-size: 1.05rem; color: var(--accent-purple); margin-top: 1.8rem;">Short-Distance Flat Background Verification (3.0 m)</h3>
            <img class="hero-img" src="../short_distance/trajectories.png" alt="Short Distance Flat Background Verification">
        </div>
    </div>

    <!-- ==================== TAB 3: SHORT DISTANCE ==================== -->
    <div id="pane-short" class="section-pane">
        <div class="card">
            <div class="card-header">
                <div>
                    <h2 class="card-title">Short Distance (3.0 m) Frame-by-Frame Reviewer</h2>
                    <p style="color: var(--text-secondary); font-size: 0.85rem; margin-top: 0.2rem;">
                        Verification of flat solvent background across $q \in [0.022, 2.54]\text{ nm}^{-1}$ ($t = 0.060$ to $10.877\text{ s}$).
                    </p>
                </div>
            </div>

            <div class="reviewer-layout">
                <div class="plot-container">
                    <img id="short-plot-img" src="../short_distance/fits/fit_frame_049.png" alt="Short Distance Fit">
                </div>

                <div class="controls-panel">
                    <div class="slider-box">
                        <div class="slider-label">
                            <span id="short-frame-label">Frame 049 / 049</span>
                            <span id="short-time-label" style="color: var(--accent-blue);">t = 10.877 s</span>
                        </div>
                        <input type="range" id="short-slider" min="0" max="49" value="49" oninput="updateShortFrame(this.value)">
                        
                        <div class="playback-bar">
                            <button class="btn" onclick="stepShortFrame(-10)">&minus;10</button>
                            <button class="btn" onclick="stepShortFrame(-1)">&minus;1</button>
                            <button class="btn btn-play" id="short-play-btn" onclick="toggleShortPlay()">Play</button>
                            <button class="btn" onclick="stepShortFrame(1)">+1</button>
                            <button class="btn" onclick="stepShortFrame(10)">+10</button>
                        </div>
                    </div>

                    <div class="telemetry-card">
                        <div style="font-weight: 600; color: var(--text-primary); margin-bottom: 0.6rem; border-bottom: 1px solid var(--border-color); padding-bottom: 0.3rem;">
                            Frame Telemetry (<span id="tele-short-f">049</span>)
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Physical Time:</span>
                            <span class="tele-val" id="tele-short-t">10.877 s</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Reduced &chi;&sup2;:</span>
                            <span class="tele-val" id="tele-short-chi2">0.094</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Fitted Bkg I<sub>bkg</sub>:</span>
                            <span class="tele-val" id="tele-short-bkg">(0.88 ± 0.11) &times; 10<sup>-4</sup> cm<sup>-1</sup></span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Mean Subtracted I:</span>
                            <span class="tele-val" id="tele-short-meani">0.04 &times; 10<sup>-4</sup> cm<sup>-1</sup></span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Positive Points:</span>
                            <span class="tele-val" id="tele-short-pospct">57.9 % (pure noise)</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Physical State:</span>
                            <span class="tele-val" style="color: var(--accent-teal);">Flat Solvent (Water)</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- ==================== TAB 4: LONG DISTANCE ==================== -->
    <div id="pane-long" class="section-pane">
        <div class="card">
            <div class="card-header">
                <div>
                    <h2 class="card-title">Long Distance (31.0 m) Frame-by-Frame Reviewer</h2>
                    <p style="color: var(--text-secondary); font-size: 0.85rem; margin-top: 0.2rem;">
                        Discrete ACC nanosphere form factor fits across $q \in [0.002, 0.25]\text{ nm}^{-1}$ ($t = 0.060$ to $10.877\text{ s}$).
                    </p>
                </div>
            </div>

            <div class="reviewer-layout">
                <div class="plot-container">
                    <img id="long-plot-img" src="../long_distance/fits/fit_frame_049.png" alt="Long Distance Fit">
                </div>

                <div class="controls-panel">
                    <div class="slider-box">
                        <div class="slider-label">
                            <span id="long-frame-label">Frame 049 / 049</span>
                            <span id="long-time-label" style="color: var(--accent-blue);">t = 10.877 s</span>
                        </div>
                        <input type="range" id="long-slider" min="0" max="49" value="49" oninput="updateLongFrame(this.value)">
                        
                        <div class="playback-bar">
                            <button class="btn" onclick="stepLongFrame(-10)">&minus;10</button>
                            <button class="btn" onclick="stepLongFrame(-1)">&minus;1</button>
                            <button class="btn btn-play" id="long-play-btn" onclick="toggleLongPlay()">Play</button>
                            <button class="btn" onclick="stepLongFrame(1)">+1</button>
                            <button class="btn" onclick="stepLongFrame(10)">+10</button>
                        </div>
                    </div>

                    <div class="telemetry-card">
                        <div style="font-weight: 600; color: var(--text-primary); margin-bottom: 0.6rem; border-bottom: 1px solid var(--border-color); padding-bottom: 0.3rem;">
                            Frame Telemetry (<span id="tele-long-f">049</span>)
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Physical Time:</span>
                            <span class="tele-val" id="tele-long-t">10.877 s</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Reduced &chi;&sup2;:</span>
                            <span class="tele-val" id="tele-long-chi2">0.190</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Mean Radius &lang;R&rang;:</span>
                            <span class="tele-val" id="tele-long-rmean" style="color: var(--accent-blue);">101.8 ± 1.2 nm</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Dispersion &sigma;<sub>R</sub>:</span>
                            <span class="tele-val" id="tele-long-stdr">23.5 ± 0.4 nm</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Polydispersity:</span>
                            <span class="tele-val" id="tele-long-poly">23.1 %</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Volume Fraction &phi;:</span>
                            <span class="tele-val" id="tele-long-phi" style="color: var(--accent-green);">4.41e-04 ± 3.39e-06</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Droplets R<sub>large</sub>, &phi;<sub>large</sub>:</span>
                            <span class="tele-val" id="tele-long-droplet" style="color: var(--accent-amber);">2167 nm, 5.25 &times; 10<sup>-5</sup></span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Porod Equiv. A<sub>porod</sub>:</span>
                            <span class="tele-val" id="tele-long-aporod" style="color: var(--accent-amber);">2.76 ± 0.06 &times; 10<sup>-8</sup></span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Solvent Background:</span>
                            <span class="tele-val" id="tele-long-bkg">0.121 ± 0.051</span>
                        </div>
                        <div class="tele-row">
                            <span class="tele-label">Regime:</span>
                            <span class="tele-val" id="tele-long-regime" style="color: var(--accent-teal);">Nanosphere Growth + Macroscopic Droplets</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- ==================== TAB 5: PARAMETER TABLES ==================== -->
    <div id="pane-tables" class="section-pane">
        <div class="card">
            <div class="card-title">Long-Distance Fitted Parameters (31.0 m, 50 Frames)</div>
            <div class="data-table-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>Frame</th>
                            <th>Time [s]</th>
                            <th>&chi;<sup>2</sup><sub>red</sub></th>
                            <th>&lang;R<sub>ACC</sub>&rang; [nm]</th>
                            <th>&sigma;<sub>ACC</sub> [nm]</th>
                            <th>&phi;<sub>ACC</sub> &times; 10<sup>4</sup></th>
                            <th>R<sub>large</sub> [nm]</th>
                            <th>&phi;<sub>large</sub> &times; 10<sup>5</sup></th>
                            <th>Bkg</th>
                            <th>Regime</th>
                        </tr>
                    </thead>
                    <tbody>
__LONG_TABLE_ROWS__
                    </tbody>
                </table>
            </div>
        </div>

        <div class="card">
            <div class="card-title">Short-Distance Fitted Parameters (3.0 m, 50 Frames)</div>
            <div class="data-table-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>Frame</th>
                            <th>Time [s]</th>
                            <th>&chi;<sup>2</sup><sub>red</sub></th>
                            <th>I<sub>bkg</sub> &times; 10<sup>4</sup> [cm<sup>-1</sup>]</th>
                            <th>Mean I &times; 10<sup>4</sup></th>
                            <th>% Positive Points</th>
                            <th>Points</th>
                        </tr>
                    </thead>
                    <tbody>
__SHORT_TABLE_ROWS__
                    </tbody>
                </table>
            </div>
        </div>
    </div>

    <!-- ==================== TAB 6: PHYSICAL FORMALISM ==================== -->
    <div id="pane-theory" class="section-pane">
        <div class="card">
            <h2 class="card-title">Physical Model Formulations &amp; Mathematical Equations</h2>
            
            <h3 style="margin-top: 1.0rem; color: var(--accent-blue);">1. Low Supersaturation Scenario vs Spinodal Decomposition</h3>
            <p style="margin-top: 0.4rem; color: var(--text-secondary);">
                At 10 mM CaCO<sub>3</sub> (ND-i-003), the reaction solution resides below the spinodal boundary.
                Unlike high-supersaturation regimes (e.g. 40 mM in ND-i-001/ND-i-009) where barrierless phase separation produces a bicontinuous mass fractal gel network ($D_m \approx 1.84$),
                ND-i-003 exhibits classical homogeneous nucleation of discrete, non-percolating amorphous calcium carbonate (ACC) nanospheres.
            </p>

            <h3 style="margin-top: 1.2rem; color: var(--accent-teal);">2. Unified Multi-Component Scattering Model</h3>
            <div class="equation-box">
                $$I_{\rm tot}(q) = I_{\rm sphere}(q) + I_{\rm droplets}(q) + I_{\rm bkg}$$
                $$I_{\rm droplets}(q) = (\Delta\rho)^2 \int_0^\infty n_{\rm large}(r) V(r)^2 P(q, r)\, dr \xrightarrow{q \gg 1/R_{\rm large}} \frac{27 \phi_{\rm large}}{8\pi R_{\rm large} q^4} = \frac{A_{\rm porod}}{q^4}$$
                Macroscopic droplets ($R_{\rm large} \approx 2.1 - 2.8\ \mu\mathrm{m}$, $\phi_{\rm large} \approx (5.0 - 6.0) \times 10^{-5}$) represent mixing/shear-induced emulsion droplets identical to the phase observed in ND-i-001 and ND-i-009.
                This form factor captures the low-$q$ forward scattering rollover ($q < 0.005\text{ nm}^{-1}$) avoiding unphysical $1/q^4$ divergence, while smoothly decaying as $q^{-4}$ at higher $q$.
            </div>

            <h3 style="margin-top: 1.2rem; color: var(--accent-blue);">3. Polydisperse Sphere Form Factor</h3>
            <div class="equation-box">
                $$I_{\rm sphere}(q) = (\Delta\rho)^2 \int_0^\infty n(r) V(r)^2 P(q, r)\, dr$$
                $$P(q, r) = \left[ \frac{3 \left(\sin(qr) - qr \cos(qr)\right)}{(qr)^3} \right]^2, \quad V(r) = \frac{4}{3}\pi r^3$$
            </div>

            <h3 style="margin-top: 1.2rem; color: var(--accent-purple);">4. Lognormal Size Distribution</h3>
            <div class="equation-box">
                $$n(r) = \frac{\phi}{\langle V \rangle} \frac{1}{r \sigma_{\rm log} \sqrt{2\pi}} \exp\left( - \frac{(\ln r - \mu_{\rm log})^2}{2\sigma_{\rm log}^2} \right)$$
                $$\mu_{\rm log} = \ln\left(\frac{\langle R \rangle^2}{\sqrt{\sigma_R^2 + \langle R \rangle^2}}\right), \quad \sigma_{\rm log} = \sqrt{\ln\left(1 + \frac{\sigma_R^2}{\langle R \rangle^2}\right)}$$
            </div>

            <h3 style="margin-top: 1.2rem; color: var(--accent-green);">5. Porod Law for Compact 3D Objects with Sharp Interfaces</h3>
            <div class="equation-box">
                $$\lim_{q \gg 1/\langle R \rangle} I_{\rm sphere}(q) = 2\pi (\Delta\rho)^2 S_{\rm spec} q^{-4} \quad (D = 3)$$
                The observation of strict $q^{-4}$ decay across $q \in [0.03, 0.10]\text{ nm}^{-1}$ confirms that the particles are smooth, dense 3D spheres rather than diffuse mass fractals ($q^{-D_m}$ with $D_m < 3$).
            </div>

            <h3 style="margin-top: 1.2rem; color: var(--accent-amber);">6. Stoichiometric Precursor Mass Balance</h3>
            <div class="equation-box">
                $$\phi_{\rm max} = \frac{[\text{Ca}^{2+}]_0 \cdot M_{\text{CaCO}_3}}{\rho_{\text{ACC}}} = \frac{0.010\text{ mol/L} \times 100.09\text{ g/mol}}{2.0\text{ g/cm}^3 \times 1000\text{ cm}^3/\text{L}} \approx 5.0 \times 10^{-4}$$
                The fitted volume fraction $\phi(t)$ plateaus at $(4.41 \pm 0.03) \times 10^{-4}$, representing near complete conversion towards the stoichiometric theoretical maximum $\phi_{\rm max} \approx 5.0 \times 10^{-4}$.
            </div>
        </div>
    </div>
</div>

<script>
    const shortData = __SHORT_DATA__;
    const longData = __LONG_DATA__;

    let shortTimer = null;
    let longTimer = null;

    function switchSection(secId) {
        document.querySelectorAll('.nav-tab').forEach(t => t.classList.remove('active'));
        document.querySelectorAll('.section-pane').forEach(p => p.classList.remove('active'));

        const pane = document.getElementById('pane-' + secId);
        if (pane) pane.classList.add('active');

        const tabs = document.querySelectorAll('.nav-tab');
        if (secId === 'overview') tabs[0].classList.add('active');
        else if (secId === 'trajectories') tabs[1].classList.add('active');
        else if (secId === 'short') tabs[2].classList.add('active');
        else if (secId === 'long') tabs[3].classList.add('active');
        else if (secId === 'tables') tabs[4].classList.add('active');
        else if (secId === 'theory') tabs[5].classList.add('active');

        if (window.MathJax && window.MathJax.typesetPromise) {
            window.MathJax.typesetPromise();
        }
    }

    // --- Long Distance Reviewer Functions ---
    function updateLongFrame(val) {
        const f = parseInt(val);
        const fStr = String(f).padStart(3, '0');
        const img = document.getElementById('long-plot-img');
        if (img) img.src = `../long_distance/fits/fit_frame_${fStr}.png`;

        const data = longData[f] || {};
        document.getElementById('long-frame-label').textContent = `Frame ${fStr} / 049`;
        document.getElementById('long-time-label').textContent = `t = ${data.time || '0.000'} s`;
        document.getElementById('long-slider').value = f;

        document.getElementById('tele-long-f').textContent = fStr;
        document.getElementById('tele-long-t').textContent = `${data.time || '0.000'} s`;
        document.getElementById('tele-long-chi2').textContent = data.chi2 !== undefined ? data.chi2 : '-';
        document.getElementById('tele-long-rmean').textContent = data.Rmean ? (data.is_ind ? data.Rmean : `${data.Rmean} nm`) : '-';
        document.getElementById('tele-long-stdr').textContent = data.stdDev ? (data.is_ind ? data.stdDev : `${data.stdDev} nm`) : '-';
        document.getElementById('tele-long-poly').textContent = data.poly_pct ? (data.is_ind ? data.poly_pct : `${data.poly_pct} %`) : '-';
        document.getElementById('tele-long-phi').textContent = data.phi || '-';
        document.getElementById('tele-long-droplet').textContent = data.droplet || '-';
        document.getElementById('tele-long-aporod').textContent = data.A_porod || '-';
        document.getElementById('tele-long-bkg').textContent = data.bkg || '-';
        
        const regimeEl = document.getElementById('tele-long-regime');
        if (data.is_ind) {
            regimeEl.textContent = 'Induction (Macroscopic Droplets)';
            regimeEl.style.color = 'var(--text-secondary)';
        } else {
            regimeEl.textContent = 'Nanosphere Growth + Macroscopic Droplets';
            regimeEl.style.color = 'var(--accent-teal)';
        }
    }

    function stepLongFrame(delta) {
        let f = parseInt(document.getElementById('long-slider').value) + delta;
        if (f < 0) f = 0;
        if (f > 49) f = 49;
        updateLongFrame(f);
    }

    function toggleLongPlay() {
        const btn = document.getElementById('long-play-btn');
        if (longTimer) {
            clearInterval(longTimer);
            longTimer = null;
            btn.textContent = 'Play';
            btn.classList.remove('btn-pause');
            btn.classList.add('btn-play');
        } else {
            btn.textContent = 'Pause';
            btn.classList.remove('btn-play');
            btn.classList.add('btn-pause');
            longTimer = setInterval(() => {
                let f = parseInt(document.getElementById('long-slider').value) + 1;
                if (f > 49) f = 0;
                updateLongFrame(f);
            }, 180);
        }
    }

    // --- Short Distance Reviewer Functions ---
    function updateShortFrame(val) {
        const f = parseInt(val);
        const fStr = String(f).padStart(3, '0');
        const img = document.getElementById('short-plot-img');
        if (img) img.src = `../short_distance/fits/fit_frame_${fStr}.png`;

        const data = shortData[f] || {};
        document.getElementById('short-frame-label').textContent = `Frame ${fStr} / 049`;
        document.getElementById('short-time-label').textContent = `t = ${data.time || '0.000'} s`;
        document.getElementById('short-slider').value = f;

        document.getElementById('tele-short-f').textContent = fStr;
        document.getElementById('tele-short-t').textContent = `${data.time || '0.000'} s`;
        document.getElementById('tele-short-chi2').textContent = data.chi2 !== undefined ? data.chi2 : '-';
        document.getElementById('tele-short-bkg').textContent = data.bkg ? `${data.bkg} × 10⁻⁴ cm⁻¹` : '-';
        document.getElementById('tele-short-meani').textContent = data.mean_I ? `${data.mean_I} × 10⁻⁴ cm⁻¹` : '-';
        document.getElementById('tele-short-pospct').textContent = data.pos_pct ? `${data.pos_pct} % (noise)` : '-';
    }

    function stepShortFrame(delta) {
        let f = parseInt(document.getElementById('short-slider').value) + delta;
        if (f < 0) f = 0;
        if (f > 49) f = 49;
        updateShortFrame(f);
    }

    function toggleShortPlay() {
        const btn = document.getElementById('short-play-btn');
        if (shortTimer) {
            clearInterval(shortTimer);
            shortTimer = null;
            btn.textContent = 'Play';
            btn.classList.remove('btn-pause');
            btn.classList.add('btn-play');
        } else {
            btn.textContent = 'Pause';
            btn.classList.remove('btn-play');
            btn.classList.add('btn-pause');
            shortTimer = setInterval(() => {
                let f = parseInt(document.getElementById('short-slider').value) + 1;
                if (f > 49) f = 0;
                updateShortFrame(f);
            }, 180);
        }
    }

    // Initialize default states
    window.addEventListener('DOMContentLoaded', () => {
        updateLongFrame(49);
        updateShortFrame(49);
    });
</script>
</body>
</html>
"""

    # Generate Long Table Rows
    long_rows = []
    for _, r in df_long.iterrows():
        f = int(r['frame'])
        t_s = float(r.get('time_seconds', 0.0))
        chi2 = float(r['redchi']) if pd.notnull(r.get('redchi')) else 0.0
        is_ind = bool(r.get('is_induction', False))
        if is_ind:
            r_str = "&mdash;"
            std_str = "&mdash;"
            phi_str = "&mdash;"
            regime = "Induction"
        else:
            r_val = float(r['Rmean'])
            r_err = float(r.get('std_err_Rmean', 0.0))
            r_str = f"{r_val:.1f} &plusmn; {r_err:.1f}" if r_err > 0 else f"{r_val:.1f}"
            std_val = float(r['stdDev'])
            std_err = float(r.get('std_err_stdDev', 0.0))
            std_str = f"{std_val:.1f} &plusmn; {std_err:.1f}" if std_err > 0 else f"{std_val:.1f}"
            phi_val = float(r['phi']) * 1e4
            phi_err = float(r.get('std_err_phi', 0.0)) * 1e4
            phi_str = f"{phi_val:.2f} &plusmn; {phi_err:.2f}" if phi_err > 0 else f"{phi_val:.2f}"
            regime = "ACC Growth"
        r_l_val = float(r['mean_radius_large']) if pd.notnull(r.get('mean_radius_large')) else 0.0
        r_l_err = float(r.get('std_err_mean_radius_large', 0.0)) if pd.notnull(r.get('std_err_mean_radius_large')) else 0.0
        r_l_str = f"{r_l_val:.0f} &plusmn; {r_l_err:.0f}" if r_l_err > 0 else f"{r_l_val:.0f}"
        phi_l_val = float(r['total_volume_fraction_large']) * 1e5 if pd.notnull(r.get('total_volume_fraction_large')) else 0.0
        phi_l_err = float(r.get('std_err_total_volume_fraction_large', 0.0)) * 1e5 if pd.notnull(r.get('std_err_total_volume_fraction_large')) else 0.0
        phi_l_str = f"{phi_l_val:.2f} &plusmn; {phi_l_err:.2f}" if phi_l_err > 0 else f"{phi_l_val:.2f}"
        bkg_val = float(r['bkg'])
        bkg_err = float(r.get('std_err_bkg', 0.0))
        bkg_str = f"{bkg_val:.3f} &plusmn; {bkg_err:.3f}" if bkg_err > 0 else f"{bkg_val:.3f}"
        long_rows.append(f"<tr><td>{f:02d}</td><td>{t_s:.3f}</td><td>{chi2:.3f}</td><td>{r_str}</td><td>{std_str}</td><td>{phi_str}</td><td>{r_l_str}</td><td>{phi_l_str}</td><td>{bkg_str}</td><td style='text-align: left;'>{regime}</td></tr>")

    # Generate Short Table Rows
    short_rows = []
    for _, r in df_short.iterrows():
        f = int(r['frame'])
        t_s = float(r.get('time_seconds', 0.0))
        chi2 = float(r['redchi']) if pd.notnull(r.get('redchi')) else 0.0
        bkg_val = float(r['bkg']) * 1e4
        bkg_err = float(r.get('std_err_bkg', 0.0)) * 1e4
        bkg_str = f"{bkg_val:.2f} &plusmn; {bkg_err:.2f}" if bkg_err > 0 else f"{bkg_val:.2f}"
        mean_i = float(r.get('mean_I', 0.0)) * 1e4
        pos_pct = float(r.get('pos_pct', 50.0))
        pts = int(r.get('pts', 1638))
        short_rows.append(f"<tr><td>{f:02d}</td><td>{t_s:.3f}</td><td>{chi2:.3f}</td><td>{bkg_str}</td><td>{mean_i:.2f}</td><td>{pos_pct:.1f}%</td><td>{pts}</td></tr>")

    # Replace Template Placeholders
    poly_pct = (bp['stdDev'] / bp['Rmean']) * 100.0
    poly2_pct = (bp.get('stdDev2', 0.0) / bp.get('R2', 1.0)) * 100.0 if bp.get('R2', 0.0) > 0 else 0.0
    poly_large_pct = (bp.get('std_radius_large', 0.0) / bp.get('mean_radius_large', 1.0)) * 100.0 if bp.get('mean_radius_large', 0.0) > 0 else 0.0

    html_content = raw_template
    html_content = html_content.replace('__SHORT_DATA__', json.dumps(short_telemetry))
    html_content = html_content.replace('__LONG_DATA__', json.dumps(long_telemetry))
    html_content = html_content.replace('__LONG_TABLE_ROWS__', '\n'.join(long_rows))
    html_content = html_content.replace('__SHORT_TABLE_ROWS__', '\n'.join(short_rows))

    # Mature ACC
    html_content = html_content.replace('__BP_RMEAN__', f"{bp['Rmean']:.2f}")
    html_content = html_content.replace('__BP_ERR_RMEAN__', f"{bp.get('std_err_Rmean', 0.0):.2f}")
    html_content = html_content.replace('__BP_STDDEV__', f"{bp['stdDev']:.2f}")
    html_content = html_content.replace('__BP_ERR_STDDEV__', f"{bp.get('std_err_stdDev', 0.0):.2f}")
    html_content = html_content.replace('__BP_POLY_PCT__', f"{poly_pct:.1f}")
    html_content = html_content.replace('__BP_PHI__', f"{bp['phi'] * 1e4:.2f}")
    html_content = html_content.replace('__BP_ERR_PHI__', f"{bp.get('std_err_phi', 0.0) * 1e4:.2f}")

    # Growing ACC
    html_content = html_content.replace('__BP_R2__', f"{bp.get('R2', 0.0):.2f}")
    html_content = html_content.replace('__BP_ERR_R2__', f"{bp.get('std_err_R2', 0.0):.2f}")
    html_content = html_content.replace('__BP_STDDEV2__', f"{bp.get('stdDev2', 0.0):.2f}")
    html_content = html_content.replace('__BP_ERR_STDDEV2__', f"{bp.get('std_err_stdDev2', 0.0):.2f}")
    html_content = html_content.replace('__BP_POLY2_PCT__', f"{poly2_pct:.1f}")
    html_content = html_content.replace('__BP_PHI2__', f"{bp.get('phi2', 0.0) * 1e4:.2f}")
    html_content = html_content.replace('__BP_ERR_PHI2__', f"{bp.get('std_err_phi2', 0.0) * 1e4:.2f}")

    # Macroscopic Droplets
    html_content = html_content.replace('__BP_RLARGE__', f"{bp.get('mean_radius_large', 0.0):.0f}")
    html_content = html_content.replace('__BP_ERR_RLARGE__', f"{bp.get('std_err_mean_radius_large', 0.0):.0f}")
    html_content = html_content.replace('__BP_STDLARGE__', f"{bp.get('std_radius_large', 0.0):.0f}")
    html_content = html_content.replace('__BP_ERR_STDLARGE__', f"{bp.get('std_err_std_radius_large', 0.0):.0f}")
    html_content = html_content.replace('__BP_POLYLARGE_PCT__', f"{poly_large_pct:.1f}")
    html_content = html_content.replace('__BP_PHILARGE__', f"{bp.get('total_volume_fraction_large', 0.0) * 1e5:.2f}")
    html_content = html_content.replace('__BP_ERR_PHILARGE__', f"{bp.get('std_err_total_volume_fraction_large', 0.0) * 1e5:.2f}")

    # Background & Goodness of Fit
    html_content = html_content.replace('__BP_BKG__', f"{bp['bkg'] * 1e5:.2f}")
    html_content = html_content.replace('__BP_ERR_BKG__', f"{bp.get('std_err_bkg', 0.0) * 1e5:.2f}")
    html_content = html_content.replace('__BP_REDCHI__', f"{bp['redchi']:.3f}")

    out_file = os.path.join(dashboard_dir, 'index.html')
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Standardized Interactive Dashboard successfully built: {out_file}")


def main():
    build_dashboard()


if __name__ == '__main__':
    main()

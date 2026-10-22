# Unified SAXS Multiscale Kinetic Analysis: ND-i-005
**System:** 40 mM $\text{CaCO}_3$, pH 11.3 with Cobalt Additive  
**Instrumentation:** High-Speed Synchrotron SAXS (Short Distance 3.0 m & Long Distance 31.0 m)  
**Reaction Timescale:** $t = 0.060\text{ s}$ to $10.877\text{ s}$ (50 time-resolved frames)  

---

## 1. Executive Summary

This repository houses the unified, time-resolved small-angle X-ray scattering (SAXS) kinetic analysis for dataset **ND-i-005** ($40\text{ mM }\text{CaCO}_3$, **pH 11.3** with **Cobalt additive**).

Using high-resolution camera geometries ($3.0\text{ m}$ short distance and $31.0\text{ m}$ long distance), the scattering curves span over 3 orders of magnitude in momentum transfer ($q \in [0.002, 2.5]\text{ nm}^{-1}$) across real reaction times from $60\text{ ms}$ to $10.88\text{ s}$.

### Key Physical Discoveries:
1. **Structural Length Scale Hierarchy:**
   - **Primary Polydisperse Spherical Building Blocks:** First-principles building blocks ($\langle R_{\rm prim} \rangle \approx 7.0 \to 9.7\text{ nm}$, $\sigma_{\rm prim} \approx 3.0 - 4.5\text{ nm}$) fitted in the first place across the primary Guinier knee ($q \le 0.40\text{ nm}^{-1}$).
   - **Viscoelastic Network Fluctuation Component:** Dynamic concentration fluctuations modelled via the generalized Fisher-Burford / Tanaka formulation ($\mu_{\rm fluc} = 4.0$) with correlation mesh $\xi_{\rm fluc} \approx 8.0 \to 9.1\text{ nm}$ and amplitude $A_{\rm fluc} \approx 0.52 \to 0.89\text{ cm}^{-1}$, capturing the remaining high-$q$ intensity.
   - **Mesoscale Mass Fractal Network:** Percolated cluster network with correlation cutoff $\xi \approx 1000 - 1100\text{ nm}$ and flowing fractal dimension $D_f \approx 2.15 \to 2.38$.
   - **Secondary Macroscopic Droplets:** Dilute secondary dense liquid droplets with mean radius $\langle R_{\rm large} \rangle \approx 3.1 - 3.8\ \mu\text{m}$ and volume fraction $\phi_{\rm large} \sim 2.5 - 4.0 \times 10^{-5}$.
2. **Two-Step Fitting Methodology:**
   - Step 1: Fit polydisperse spherical building blocks + mass fractal network in the first place on $q \le 0.40\text{ nm}^{-1}$.
   - Step 2: Explain remaining intensity in the fit with the viscoelastic fluctuation component ($A_{\rm fluc}, \xi_{\rm fluc}, \mu_{\rm fluc} = 4.0$).
   - Step 3: Joint refinement across the full range with sequential continuity regularization to prevent unphysical jumping.
   - Step 4: Full covariance matrix calculation ($\pm 1\sigma$ standard errors) at every single frame.
3. **Rigorous Parameter Uncertainties:**
   - Every single parameter at every frame is reported with exact $\pm 1\sigma$ standard errors computed from the covariance matrix of the converged fit.

---

## 2. Directory Structure

```
unified_fitting_ND-i-005/
├── data/
│   ├── raw_short_distance/          # 50 raw short-distance .dat files (scan 00252)
│   ├── raw_long_distance/           # 50 raw long-distance .dat files (scan 00720)
│   ├── raw_final_stages/            # Final stage files (00253, 00721, 00722)
│   ├── dataNDi005_fitting_results_short_distance.csv # Wide-format short data (q: 0.022-2.54 nm^-1)
│   ├── dataNDi005_fitting_results_long_distance.csv  # Wide-format long data (q: 0.002-0.25 nm^-1)
│   ├── dataNDi005_fitting_results_composite.csv      # Stitched wide-angle composite anchor
│   ├── composite_final.dat                           # Stitched ASCII curve
│   └── frame_metadata.csv                            # Frame indices, timestamps (0.060s to 10.877s)
├── models/
│   ├── scattering_models.py         # Vectorized form factors, structure factors, decompositions
│   └── fitting_engine.py            # Hierarchical optimizer & objective functions
├── scripts/
│   ├── 00_prepare_data_and_composite.py       # Data staging & composite generation
│   ├── 00_fit_composite_baseline.py          # 4-stage hierarchical baseline fit
│   ├── 01_fit_short_distance_models.py        # Short-distance time-resolved modelling
│   ├── 02_fit_long_distance_models.py         # Long-distance time-resolved modelling
│   ├── 03_generate_all_frame_plots.py         # Parallel diagnostic frame renderer (100 plots)
│   ├── 04_compare_parameter_trajectories.py   # Cross-scale trajectory synthesis
│   ├── 05_build_unified_dashboard.py          # Master interactive dashboard builder
│   └── run_unified_pipeline.py                # End-to-end pipeline orchestrator
├── output/
│   ├── composite/
│   │   ├── baseline_parameters.json           # Baseline parameters with +/- 1 sigma
│   │   ├── baseline_parameters.csv
│   │   └── composite_baseline_fit.png         # Publication baseline plot
│   ├── short_distance/
│   │   ├── parameters.csv                     # 50 frames with +/- 1 sigma uncertainties
│   │   ├── fitted_curves.csv                  # Calculated full-range scattering profiles
│   │   ├── trajectories.png                   # 4-panel kinetic trajectory figure
│   │   └── fits/                              # 50 diagnostic frame fit PNGs
│   ├── long_distance/
│   │   ├── parameters.csv                     # 50 frames with +/- 1 sigma uncertainties
│   │   ├── fitted_curves.csv                  # Calculated full-range scattering profiles
│   │   ├── trajectories.png                   # 4-panel kinetic trajectory figure
│   │   └── fits/                              # 50 diagnostic frame fit PNGs
│   ├── unified_trajectories/
│   │   └── unified_cross_scale_trajectories.png # 4-panel cross-scale synthesis figure
│   └── dashboard/
│       └── index.html                         # Master standalone interactive reviewer UI
└── README.md
```

---

## 3. Converged Baseline Parameters ($\pm 1\sigma$)

| Parameter | Symbol | Value | $\pm 1\sigma$ Error | Physical Meaning |
| :--- | :--- | :--- | :--- | :--- |
| Primary Core Radius | $\langle R_{\rm prim} \rangle$ | $6.59\text{ nm}$ | $\pm 0.53\text{ nm}$ | Mean radius of amorphous building blocks |
| Primary Core Dispersion | $\sigma_{\rm prim}$ | $3.32\text{ nm}$ | $\pm 0.13\text{ nm}$ | Log-normal standard deviation |
| Primary Volume Fraction | $\phi_{\rm prim}$ | $9.7 \times 10^{-5}$ | $\pm 0.3 \times 10^{-5}$ | Volume fraction of primary fractal network |
| Mass Fractal Dimension | $D_f$ | $2.102$ | $\pm 0.013$ | Space-filling scaling exponent |
| Network Correlation Cutoff | $\xi$ | $1194\text{ nm}$ | $\pm 2347\text{ nm}$ | Mesoscale network cluster size |
| Droplet Radius | $\langle R_{\rm large} \rangle$ | $3530\text{ nm}$ | $\pm 48\text{ nm}$ | Mean radius of secondary dense liquid droplets |
| Droplet Dispersion | $\sigma_{\rm large}$ | $504\text{ nm}$ | $\pm 6.3\text{ nm}$ | Log-normal standard deviation ($\sigma/R = 0.143$) |
| Droplet Volume Fraction | $\phi_{\rm large}$ | $2.5 \times 10^{-5}$ | $\pm 0.5 \times 10^{-5}$ | Volume fraction of secondary droplets |
| Cluster Radius | $R_{\rm pnc}$ | $2.40\text{ nm}$ | $\pm 0.03\text{ nm}$ | Radius of sub-nanometer clusters / DOLLOPs |
| Cluster Dispersion | $\sigma_{\rm pnc}$ | $1.32\text{ nm}$ | $\pm 0.01\text{ nm}$ | Log-normal standard deviation |
| Cluster Volume Fraction | $\phi_{\rm pnc}$ | $4.87 \times 10^{-4}$ | $\pm 0.08 \times 10^{-4}$ | Volume fraction of cluster population |
| Incoherent Background | $bkg$ | $8.43 \times 10^{-4}$ | $\pm 0.16 \times 10^{-4}$ | Constant flat background [cm$^{-1}$] |
| Goodness of Fit | $\chi^2_{\rm red}$ | **2.67** | : | Reduced chi-square across 1831 points |

---

## 4. Master Interactive Dashboard

The interactive reviewer dashboard is located at:
`output/dashboard/index.html`

Features:
- **Tab 1: Overview & Baseline Anchor:** Full composite scattering curve and parameter summary table.
- **Tab 2: Cross-Scale Kinetics:** 4-decade hierarchy of length scales and cross-camera $D_f(t)$ congruence.
- **Tab 3 & 4: Time-Resolved Inspectors:** Real-time slider and playback controls for all 50 frames with live telemetry tables reporting $\pm 1\sigma$ standard errors.
- **Tab 5: Parameter Tables:** Scrollable tables with frame-by-frame values and uncertainties.
- **Tab 6: Physical Formulations:** Native MathJax rendering of scattering equations.

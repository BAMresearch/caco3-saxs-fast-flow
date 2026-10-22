# Unified SAXS Multiscale Kinetic Analysis: ND-i-002
**System:** 20 mM $\text{CaCO}_3$, pH 11.3 (Pure System without Additive)  
**Instrumentation:** High-Speed Synchrotron SAXS (Short Distance 3.0 m & Long Distance 31.0 m, ESRF ID02)  
**Reaction Timescale:** $t = 0.060\text{ s}$ to $10.877\text{ s}$ (50 time-resolved frames)  

---

## 1. Executive Summary

This repository houses the unified, time-resolved small-angle X-ray scattering (SAXS) kinetic analysis suite for dataset **ND-i-002** ($20\text{ mM }\text{CaCO}_3$, **pH 11.3**, pure homogeneous precipitation without additive).

Using synchronized high-speed synchrotron camera geometries ($3.0\text{ m}$ short distance and $31.0\text{ m}$ long distance), the scattering profiles span over 3 orders of magnitude in momentum transfer ($q \in [0.002, 2.54]\text{ nm}^{-1}$) across real reaction times from $60\text{ ms}$ to $10.88\text{ s}$.

### Key Physical Discoveries:
1. **Structural Hierarchy Across Multiple Spatial Scales:**
   - **Primary Polydisperse Spherical Building Blocks:** First-principles building blocks ($\langle R_{\rm prim} \rangle \approx 8.5 \to 11.4\text{ nm}$, $\sigma_{\rm prim} \approx 4.0 - 4.4\text{ nm}$) fitted across the primary knee ($q \le 0.40\text{ nm}^{-1}$).
   - **Viscoelastic Network Fluctuation Component:** Dynamic concentration fluctuations modelled via the generalized Fisher-Burford / Tanaka formulation ($\mu_{\rm fluc} = 4.0$) with correlation mesh $\xi_{\rm fluc} \approx 13.5 \to 15.0\text{ nm}$ and amplitude $A_{\rm fluc} \approx 0.03 \to 0.16\text{ cm}^{-1}$, fully explaining the high-$q$ excess intensity.
   - **Mesoscale Mass Fractal Network:** Percolated cluster network with correlation cutoff $\xi \approx 400\text{ nm}$ and fractal dimension $D_f \approx 1.30 \to 1.97$.
   - **Secondary Macroscopic Dense Liquid Droplets:** Macroscopic dense phase droplets with $\langle R_{\rm large} \rangle \approx 5,200 \to 5,880\text{ nm}$ under constant relative polydispersity constraint ($\sigma_{\rm large}/R_{\rm large} = 0.0760$).
2. **Sequential Evolutionary Regularization & Robustness:**
   - Point-to-point parameter continuity regularization strictly prevents unphysical frame-to-frame jumping while capturing genuine physical kinetics, and gracefully resists single-frame beamline transient events (e.g. bubble/dust spike at $t \approx 1.85\text{ s}$).
3. **Rigorous Parameter Uncertainties:**
   - Every single parameter at every frame is reported with exact $\pm 1\sigma$ standard errors computed from the covariance matrix of the non-linear least-squares estimator (`lmfit`).

---

## 2. Directory Structure

```
unified_fitting_ND-i-002/
├── README.md                                          # Comprehensive documentation
├── data/
│   ├── raw_short_distance/                            # 50 raw short-distance .dat files (scan 00128)
│   ├── raw_long_distance/                             # 50 raw long-distance .dat files (scan 00564)
│   ├── raw_final_stages/                              # Final stage files (00129, 00565)
│   ├── dataNDi002_fitting_results_short_distance.csv  # Wide-format short data (q: 0.022-2.54 nm^-1)
│   ├── dataNDi002_fitting_results_long_distance.csv   # Wide-format long data (q: 0.002-0.25 nm^-1)
│   ├── dataNDi002_fitting_results_composite.csv       # Stitched wide-angle composite anchor
│   ├── composite_final.dat                            # Stitched ASCII curve
│   └── frame_metadata.csv                             # Frame indices and timestamps (0.060s to 10.877s)
├── models/
│   ├── __init__.py
│   ├── scattering_models.py                           # Vectorized form factors, structure factors, decompositions
│   └── fitting_engine.py                              # Evolutionary parameter inheritance & objective functions
├── scripts/
│   ├── 00_prepare_data_and_composite.py              # Data staging & composite generation
│   ├── 00_fit_composite_baseline.py                  # 4-stage hierarchical baseline fit
│   ├── 01_fit_short_distance_models.py               # Short-distance time-resolved modelling
│   ├── 02_fit_long_distance_models.py                # Long-distance time-resolved modelling
│   ├── 03_generate_all_frame_plots.py                # Parallel diagnostic frame renderer (100 plots)
│   ├── 04_compare_parameter_trajectories.py          # Cross-scale trajectory synthesis
│   ├── 05_build_unified_dashboard.py                 # Master interactive dashboard builder
│   └── run_unified_pipeline.py                       # End-to-end pipeline orchestrator
└── output/
    ├── composite/
    │   ├── baseline_parameters.json                  # Baseline parameters with +/- 1 sigma
    │   ├── baseline_parameters.csv
    │   └── composite_baseline_fit.png                # Publication baseline plot
    ├── short_distance/
    │   ├── parameters.csv                            # 50 frames with +/- 1 sigma uncertainties
    │   ├── fitted_curves.csv                         # Calculated full-range scattering profiles
    │   ├── trajectories.png                          # 4-panel kinetic trajectory figure
    │   └── fits/                                     # 50 diagnostic frame fit PNGs
    ├── long_distance/
    │   ├── parameters.csv                            # 50 frames with +/- 1 sigma uncertainties
    │   ├── fitted_curves.csv                         # Calculated full-range scattering profiles
    │   ├── trajectories.png                          # 4-panel kinetic trajectory figure
    │   └── fits/                                     # 50 diagnostic frame fit PNGs
    ├── unified_trajectories/
    │   └── unified_cross_scale_trajectories.png        # 4-panel cross-scale synthesis figure
    └── dashboard/
        └── index.html                                # Master standalone interactive reviewer UI
```

---

## 3. Converged Baseline Parameters ($\pm 1\sigma$)

| Parameter | Symbol | Value | $\pm 1\sigma$ Error | Physical Meaning |
| :--- | :--- | :--- | :--- | :--- |
| Primary Building Block Radius | $\langle R_{\rm prim} \rangle$ | $11.24\text{ nm}$ | $\pm 0.13\text{ nm}$ | Mean radius of amorphous building blocks |
| Primary Core Dispersion | $\sigma_{\rm prim}$ | $3.99\text{ nm}$ | $\pm 0.02\text{ nm}$ | Log-normal standard deviation |
| Primary Network Volume Fraction | $\phi_{\rm net}$ | $4.74 \times 10^{-5}$ | $\pm 0.02 \times 10^{-5}$ | Volume fraction of primary fractal network |
| Mass Fractal Dimension | $D_m$ | $1.856$ | $\pm 0.007$ | Space-filling network scaling exponent |
| Network Correlation Cutoff | $\xi$ | $356.0\text{ nm}$ | $\pm 24.8\text{ nm}$ | Mesoscale network cluster correlation length |
| Porod Interface Scale | $A_{\rm porod}$ | $5.28 \times 10^{-8}$ | $\pm 0.02 \times 10^{-8}$ | Macroscopic surface interface scattering |
| Prenucleation Cluster Radius | $\langle R_{\rm pnc} \rangle$ | $1.46\text{ nm}$ | $\pm 0.39\text{ nm}$ | Mean radius of sub-nm PNC clusters |
| PNC Volume Fraction | $\phi_{\rm pnc}$ | $1.50 \times 10^{-5}$ | $\pm 0.12 \times 10^{-5}$ | Volume fraction of PNC cluster population |
| Incoherent Background | $bkg$ | $1.37 \times 10^{-4}$ | $\pm 0.08 \times 10^{-4}$ | Constant flat background [cm$^{-1}$] |
| Goodness of Fit | $\chi^2_{\rm red}$ | **0.964** | : | Reduced chi-square across 1879 points |

---

## 4. Pipeline Execution

To run the complete workflow from the project directory:

```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --all
```

To run individual steps:
```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 0   # Data preparation & composite anchor
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 1   # Baseline composite fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 2   # Short-distance time-resolved fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 3   # Long-distance time-resolved fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 4   # Universal frame diagnostic plots (100 plots)
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 5   # Cross-scale trajectory synthesis
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 6   # Compile master interactive dashboard
```

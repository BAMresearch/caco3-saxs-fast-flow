# Unified SAXS multi-scale kinetic analysis: ND-i-001
**System:** 40 mM $\text{CaCO}_3$, pH 11.3 (Pure System without Additive)  
**Instrumentation:** High-Speed Synchrotron SAXS (Short Distance 3.0 m & Long Distance 31.0 m, ESRF ID02)  
**Reaction Timescale:** $t = 0.060\text{ s}$ to $10.877\text{ s}$ (50 time-resolved frames)  

---

## 1. Executive summary

This repository houses the unified, time-resolved small-angle X-ray scattering (SAXS) kinetic analysis suite for dataset **ND-i-001** ($40\text{ mM }\text{CaCO}_3$, **pH 11.3**, pure homogeneous precipitation without additive).

Using synchrotron camera geometries ($3.0\text{ m}$ short distance and $31.0\text{ m}$ long distance), the scattering profiles span over 3 orders of magnitude in momentum transfer ($q \in [0.002, 2.54]\text{ nm}^{-1}$) across real reaction times from $60\text{ ms}$ to $10.88\text{ s}$.

### Key physical discoveries:
1. **Structural Length Scale Hierarchy Across 4 Decades:**
   - **Primary Polydisperse Spherical Building Blocks:** First-principles building blocks ($\langle R_{\rm prim} \rangle \approx 6.16 \to 6.78\text{ nm}$, $\sigma_{\rm prim} \approx 2.5 - 3.4\text{ nm}$) fitted across the primary Guinier knee ($q \le 0.40\text{ nm}^{-1}$).
   - **Viscoelastic Network Fluctuation Component:** Dynamic concentration fluctuations modelled via the generalized Fisher-Burford / Tanaka formulation ($\mu_{\rm fluc} = 4.0$) with correlation mesh $\xi_{\rm fluc} \approx 3.03 \to 3.05\text{ nm}$ and amplitude $A_{\rm fluc} \approx 0.08 \to 0.12\text{ cm}^{-1}$, fully explaining the high-$q$ excess intensity.
   - **Mesoscale Mass Fractal Network:** Percolated cluster network with correlation cutoff $\xi \approx 735\text{ nm}$ and flowing fractal dimension $D_f \approx 1.86 \to 1.98$.
   - **Secondary Macroscopic Dense Liquid Droplets:** Coarsening macroscopic dense phase droplets evolving from $\langle R_{\rm large} \rangle \approx 3,280\text{ nm}$ to $3,912\text{ nm}$ under constant relative polydispersity constraint ($\sigma_{\rm large}/R_{\rm large} = 0.0870$).
2. **Sequential Evolutionary Regularization:**
   - Point-to-point parameter continuity regularization strictly prevents unphysical frame-to-frame jumping while capturing genuine physical kinetics.
3. **Rigorous Parameter Uncertainties:**
   - Every single parameter at every frame is reported with exact $\pm 1\sigma$ standard errors computed from the covariance matrix of the non-linear least-squares estimator (`lmfit`).

---

## 2. Directory structure

```
unified_fitting_ND-i-001/
├── README.md                                          # Comprehensive documentation
├── data/
│   ├── raw_short_distance/                            # 50 raw short-distance .dat files (scan 00093)
│   ├── raw_long_distance/                             # 50 raw long-distance .dat files (scan 00679)
│   ├── raw_final_stages/                              # Final stage files (00674, 00675, 00680, 00681, 00093_00050)
│   ├── dataNDi001_fitting_results_short_distance.csv  # Wide-format short data (q: 0.022-2.54 nm^-1)
│   ├── dataNDi001_fitting_results_long_distance.csv   # Wide-format long data (q: 0.002-0.25 nm^-1)
│   ├── dataNDi001_fitting_results_composite.csv       # Stitched wide-angle composite anchor
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

## 3. Converged baseline parameters ($\pm 1\sigma$)

| Parameter | Symbol | Value | $\pm 1\sigma$ error | Physical meaning |
| :--- | :--- | :--- | :--- | :--- |
| Primary core radius | $\langle R_{\rm prim} \rangle$ | $3.50\text{ nm}$ | $\pm 0.28\text{ nm}$ | Mean radius of amorphous building blocks |
| Primary core dispersion | $\sigma_{\rm prim}$ | $5.18\text{ nm}$ | $\pm 0.06\text{ nm}$ | Log-normal standard deviation |
| Primary volume fraction | $\phi_{\rm prim}$ | $5.0 \times 10^{-5}$ | $\pm 0.1 \times 10^{-5}$ | Volume fraction of primary fractal network |
| Mass fractal dimension | $D_f$ | $2.328$ | $\pm 0.003$ | Space-filling scaling exponent |
| Network correlation cutoff | $\xi$ | $735.4\text{ nm}$ | $\pm 42.4\text{ nm}$ | Mesoscale network cluster size |
| Droplet radius | $\langle R_{\rm large} \rangle$ | $3015.0\text{ nm}$ | $\pm 8.8\text{ nm}$ | Mean radius of secondary dense liquid droplets |
| Droplet dispersion | $\sigma_{\rm large}$ | $262.3\text{ nm}$ | $\pm 15.2\text{ nm}$ | Log-normal standard deviation ($\sigma/R = 0.087$) |
| Droplet volume fraction | $\phi_{\rm large}$ | $7.0 \times 10^{-6}$ | $\pm 0.3 \times 10^{-6}$ | Volume fraction of secondary droplets |
| Cluster radius | $R_{\rm pnc}$ | $6.37\text{ nm}$ | $\pm 0.23\text{ nm}$ | Radius of sub-nanometer clusters / high-$q$ excess |
| Cluster dispersion | $\sigma_{\rm pnc}$ | $3.13\text{ nm}$ | $\pm 0.01\text{ nm}$ | Log-normal standard deviation |
| Cluster volume fraction | $\phi_{\rm pnc}$ | $3.80 \times 10^{-4}$ | $\pm 0.05 \times 10^{-4}$ | Volume fraction of cluster population |
| Incoherent background | $bkg$ | $1.28 \times 10^{-3}$ | $\pm 0.36 \times 10^{-3}$ | Constant flat background [cm$^{-1}$] |
| Goodness of fit | $\chi^2_{\rm red}$ | **0.398** | : | Reduced chi-square across 1540 points |

---

## 4. Pipeline execution

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

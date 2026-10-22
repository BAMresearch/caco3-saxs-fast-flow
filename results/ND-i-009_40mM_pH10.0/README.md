# Unified Multi-Scale SAXS Fitting Suite (Dataset ND-i-009)

This repository contains the clean, unified SAXS kinetic fitting suite for dataset **ND-i-009** ($40.0\text{ mM }\text{CaCO}_3$, pH 10.0, pure system without additive).

It integrates the validated theoretical concepts across short distance ($3.0\text{ m}$ camera, 351 frames), long distance ($31.0\text{ m}$ camera, 350 frames), and full-range composite scattering curves.

---

## Architecture & Scientific Models

### 1. Short Distance ($q \in [0.022, 2.54]\text{ nm}^{-1}$, 351 Frames, $t = 0.045 \to 5.295\text{ s}$)
* **Option 2: Continuous Evolutionary Viscoelastic Network Fluctuations**:
  - **Frame 0**: Finite pre-percolation cutoff ($\xi = 46.7\text{ nm}$).
  - **Frames $> 0$**: Infinite network percolated regime ($\xi \to \infty$).
  - **Point-to-Point Continuous $D_f(t)$**: Genuinely fitted for every individual frame with sequential continuity regularization to strictly eliminate discrete jumps ($|dD_f| \le 0.0040$, mean $|dD_f| = 0.00041$, range $1.8810 \to 1.9254$).
  - **Continuous Viscoelastic Drift**: Both correlation length $\xi_{\rm fluc}(t)$ ($2.72 \to 2.84\text{ nm}$) and fluctuation amplitude $A_{\rm fluc}(t)$ ($0.142 \to 0.197\text{ cm}^{-1}$) evolve smoothly without zig-zags or spikes.
  - **Quality of Fit**: Mean reduced $\chi^2_{\rm red} = 0.290$ (median $0.205$), high-$q$ $\chi^2_{\rm red} = 0.075$.

### 2. Long Distance ($q \in [0.007, 0.44]\text{ nm}^{-1}$, 350 Frames, $t = 0.045 \to 5.280\text{ s}$)
* **Option C Forward Model**:
  - **Frame 0**: Early induction phase fitted without macroscopic droplets ($\phi_{\rm large} = 0$).
  - **Frames $> 0$**: Forward continuity model with infinite network cutoff ($\xi \to \infty$) and constant relative polydispersity ($\sigma_{\rm large}/R_{\rm large} = 0.065$).
  - **Coarsening Droplets**: Droplet radius evolves from $2,274\text{ nm}$ to $2,943\text{ nm}$ under constant polydispersity constraint, eliminating false spikes in $D_f(t)$.
  - **Effective Mesoscopic Exponent**: Under the $\xi \to \infty$ constraint, $D_f(t)$ serves as an effective mesoscopic power-law exponent that reflects cluster growth during pre-percolation ($1.46 \to 1.89$ for $t < 1.0\text{ s}$), converging to $1.90$ (matching the short-distance fractal dimension) once macroscopic network percolation is reached.
  - **Quality of Fit**: Mean reduced $\chi^2_{\rm red} = 0.162$ (median $0.151$).

### 3. Rigorous Parameter Uncertainties
* All fitted parameters include formal standard errors ($\pm 1\sigma$ confidence limits) propagated from the covariance matrix of the non-linear least-squares estimator (`lmfit`).
* All trajectory figures (`trajectories.png` and `unified_cross_scale_trajectories.png`) feature shaded $\pm 1\sigma$ uncertainty bands.
* The interactive dashboard displays live $\pm 1\sigma$ error margins on all parameters across both camera distances.

---

## Directory Organization

```
unified_fitting_ND-i-009/
├── README.md                                       # Project documentation
├── data/
│   ├── dataNDi009_fitting_results_short_distance.csv # 3.0 m camera dataset (351 frames)
│   ├── dataNDi009_fitting_results_long_distance.csv  # 31.0 m camera dataset (350 frames)
│   ├── dataNDi009_fitting_results_merged.csv         # Full-range composite curve
│   └── merged.dat                                   # Original merged scattering curve
├── models/
│   ├── __init__.py
│   ├── scattering_models.py                         # Vectorized form factors & structure factors
│   └── fitting_engine.py                            # Evolutionary parameter inheritance & bounds
├── scripts/
│   ├── 00_fit_merged_baseline.py                    # Anchor full-range composite curve fit
│   ├── 01_fit_short_distance_models.py              # Option 2 short-distance fit with flowing Df
│   ├── 02_fit_long_distance_models.py               # Option C forward long-distance fit
│   ├── 03_generate_all_frame_plots.py               # Diagnostic fit plots for all frames
│   ├── 04_compare_parameter_trajectories.py         # Cross-scale kinetic trajectory synthesis
│   ├── 05_build_unified_dashboard.py                # Master interactive HTML reviewer
│   └── run_unified_pipeline.py                      # End-to-end pipeline runner
└── output/
    ├── merged/                                      # Merged baseline parameters & fit plot
    ├── short_distance/
    │   └── viscoelastic/                            # Option 2 parameters.csv, curves.csv, trajectories.png, fits/ (351 PNGs)
    ├── long_distance/
    │   └── option_c_forward/                        # Option C parameters.csv, curves.csv, trajectories.png, fits/ (350 PNGs)
    ├── unified_trajectories/                        # unified_cross_scale_trajectories.png
    └── dashboard/
        └── index.html                               # Master interactive dashboard
```

---

## Pipeline Execution

To run the complete workflow from the project directory:

```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py
```

To run individual steps:
```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 1   # Short-distance fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 2   # Long-distance fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 3   # Generate frame plots
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 4   # Update trajectories
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 5   # Build interactive dashboard
```

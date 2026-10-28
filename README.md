# Multi-scale in-situ SAXS analysis of calcium carbonate precipitation

Publication repository and unified multi-scale kinetic modelling scripts for time-resolved synchrotron Small-Angle X-ray Scattering (SAXS) measurements of calcium carbonate ($\mathrm{CaCO}_3$) precipitation collected at the **ESRF ID02 SAXS beamline** (European Synchrotron Radiation Facility, Grenoble, France).

This research is conducted within the **DFG-funded Heisenberg project (WO 1712/7-1)** and the associated **DFG project (WO 1712/8-1)**, under the scientific leadership of **Stephan E. Wolf (Principal Investigator, FAU Erlangen–Nürnberg)**. **Tomasz M. Stawski (BAM)** is responsible for SAXS data analysis and multi-scale modelling, with experimental support from **Thi Quynh Nhi Dao (FAU Erlangen–Nürnberg)**.

The project resolves the structural hierarchy and kinetic pathway across four decades in momentum transfer ($q \in [0.002, 2.54]\ \mathrm{nm}^{-1}$, real-space dimensions $2.5\ \mathrm{nm}$ to $3.1\ \mu\mathrm{m}$) and milliseconds to seconds in reaction time ($t = 0.045\ \mathrm{s}$ to $10.88\ \mathrm{s}$), discriminating between sub-spinodal nucleation of amorphous calcium carbonate (ACC) nanospheres, continuous spinodal gelation of viscoelastic mass-fractal networks, and prenucleation cluster (PNC) fluctuation mesh dynamics.

---

## Table of contents
1. [Experimental and beamline parameters](#experimental-and-beamline-parameters)
2. [Reaction series and experimental matrix](#reaction-series-and-experimental-matrix)
3. [Multi-scale physical scattering formalism](#multi-scale-physical-scattering-formalism)
4. [Standardised component palette and visual register](#standardised-component-palette-and-visual-register)
5. [Master dashboards and publication diagrams](#master-dashboards-and-publication-diagrams)
6. [Repository directory structure](#repository-directory-structure)
7. [Reproducibility and execution guide](#reproducibility-and-execution-guide)
8. [Data processing and citation](#data-processing-and-citation)
9. [License](#license)

---

## Experimental and beamline parameters

| Parameter | Specification | Details / operation notes (verified from HDF5 metadata) |
| :--- | :--- | :--- |
| **Synchrotron Light Source** | **ESRF - The European Synchrotron** (Grenoble, France) | Extremely Brilliant Source (EBS) storage ring |
| **Beamline** | **ID02 (Time-Resolved Ultra-Small-Angle & Small-Angle X-ray Scattering)** | High-brilliance undulator beamline |
| **X-ray Photon Energy** | **$12.23\ \mathrm{keV}$ ($12,230.0\ \mathrm{eV}$)** | $\mathrm{Si}(111)$ double crystal monochromator (`egy=12230.0eV`, `theta=9.308deg`) |
| **X-ray Wavelength ($\lambda$)** | **$0.101377\ \mathrm{nm}$** ($1.0138\ \mathrm{Å}$) | Calibrated pyFAI wavelength parameter (`WaveLength = 1.01377e-10 m`) |
| **Detector System** | **Dectris Eiger2 4M** | Hybrid photon counting (HPC) pixel detector |
| **Detector Dimensions** | $2068 \times 2162\ \mathrm{pixels}$ | Active sensor area $\approx 155.1 \times 162.15\ \mathrm{mm}^2$ |
| **Pixel Pitch** | $75.0 \times 75.0\ \mu\mathrm{m}^2$ | Isotropic square pixel geometry |
| **Short-distance camera length (SD)** | **$3.0\ \mathrm{m}$ nominal ($2.9997\ \mathrm{m}$ calibrated)** | HDF5 `SampleDistance = 2.9997 m`; covers $q \in [0.022, 2.54]\ \mathrm{nm}^{-1}$ |
| **Long-distance camera length (LD)** | **$31.0\ \mathrm{m}$ nominal ($30.7059\ \mathrm{m}$ calibrated)** | HDF5 `SampleDistance = 30.7059 m`; covers $q \in [0.002, 0.25]\ \mathrm{nm}^{-1}$ |
| **Combined composite $q$-range** | **$0.0020\ \mathrm{nm}^{-1}$ to $2.54\ \mathrm{nm}^{-1}$** | $>3$ orders of magnitude in reciprocal space; stitched across cameras |
| **Sample Environment / Mixer** | **BioLogic SFM-400 / 4-Syringe Stopped-Flow Mixer** | High-speed pneumatic injection with $1.3\ \mathrm{mm}$ optical path capillary |
| **Volumetric Injection Rate** | **$3.33\ \mathrm{mL/s}$** per syringe | $100\ \mu\mathrm{L}$ to $200\ \mu\mathrm{L}$ shot volumes into quartz flow capillary |
| **Mixer Dead Time** | $\sim 2.0 - 4.5\ \mathrm{ms}$ | Rapid turbulent mixing preceding exposure window |
| **Time Resolution** | $\Delta t \ge 15.0\ \mathrm{ms}$ | Geometric and logarithmic temporal spacing ($t = 0.045\ \mathrm{s} \to 10.88\ \mathrm{s}$) |
| **Intensity Calibration Scale** | **Absolute Differential Cross-Section ($\mathrm{cm}^{-1}$)** | Normalised against calibrated water capillary baseline standards |
| **Data Reduction Pipeline** | **ESRF ID02 SAXS Standard Pipeline (`SAXSutilities`)** | Azimuthal integration, detector gap mask, flux and transmission normalisation |

---

## Reaction series and experimental matrix

The experimental beamtime systematically investigates the boundary between classical nucleation and non-classical spinodal decomposition across concentration, pH, and additive perturbation:

| Series ID | Calcium reagent | Carbonate reagent | Buffer / initial pH | Physical regime | Converged model description |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **ND-i-003** | $10.0\ \mathrm{mM}\ \mathrm{CaCl}_2$ | $10.0\ \mathrm{mM}\ \mathrm{Na}_2\mathrm{CO}_3$ | Pure, $\mathrm{pH}\ 11.3$ | Sub-spinodal nucleation | Discrete ACC spheres ($R \approx 18.4\ \mathrm{nm}$, 100% volume) + macroscopic emulsion droplets |
| **ND-i-004** | $15.0\ \mathrm{mM}\ \mathrm{CaCl}_2$ | $15.0\ \mathrm{mM}\ \mathrm{Na}_2\mathrm{CO}_3$ | Pure, $\mathrm{pH}\ 11.3$ | Crossover boundary | Coexisting early mass fractal network ($D_f = 1.95$) + emerging ACC spheres |
| **ND-i-002** | $20.0\ \mathrm{mM}\ \mathrm{CaCl}_2$ | $20.0\ \mathrm{mM}\ \mathrm{Na}_2\mathrm{CO}_3$ | Pure, $\mathrm{pH}\ 11.3$ | Spinodal gelation | Percolated mass fractal network ($D_f = 1.94$, 76% volume) + Tanaka PNC mesh (24% volume) |
| **ND-i-009** | $40.0\ \mathrm{mM}\ \mathrm{CaCl}_2$ | $40.0\ \mathrm{mM}\ \mathrm{NaHCO}_3/\mathrm{Na}_2\mathrm{CO}_3$ | Buffer, $\mathrm{pH}\ 10.0$ | Primary pure benchmark | Continuous spinodal network ($D_f = 1.91$, 74% volume) + Tanaka PNC mesh (26% volume) |
| **ND-i-001** | $40.0\ \mathrm{mM}\ \mathrm{CaCl}_2$ | $40.0\ \mathrm{mM}\ \mathrm{Na}_2\mathrm{CO}_3$ | Pure, $\mathrm{pH}\ 11.3$ | Rapid quench spinodal | Instantaneous spinodal gel ($D_f = 1.98$, 88% volume) + coarsening liquid droplets |
| **ND-i-005** | $40.0\ \mathrm{mM}\ \mathrm{CoCl}_2$ | $40.0\ \mathrm{mM}\ \mathrm{Na}_2\mathrm{CO}_3$ | Pure, $\mathrm{pH}\ 11.3$ | Additive perturbation | Extended fractal network with dopant-expanded nodes ($R_{\mathrm{prim}} \approx 8.7\ \mathrm{nm}$) |
| **ND-i-006** | $40.0\ \mathrm{mM}\ \mathrm{CaCl}_2$ | $40.0\ \mathrm{mM}\ \mathrm{NaHCO}_3/\mathrm{Na}_2\mathrm{CO}_3$ | Buffer, $\mathrm{pH}\ 9.75$ | pH-moderated spinodal | Dynamic spinodal gelation with suppressed droplet coarsening rate |

---

## Multi-scale physical scattering formalism

The total coherent scattering cross-section $I_{\mathrm{tot}}(q)\ [\mathrm{cm}^{-1}]$ over $q \in [0.002, 2.54]\ \mathrm{nm}^{-1}$ is modelled as an analytical linear superposition of decoupled structural contributions:

$$I_{\mathrm{tot}}(q) = I_{\mathrm{droplet/ACC}}(q) + I_{\mathrm{network}}(q) + I_{\mathrm{fluc}}(q) + I_{\mathrm{bkg}}$$

### 1. Primary mass fractal network ($I_{\mathrm{network}}$)
Formed by polydisperse spherical building blocks obeying a lognormal distribution $n(R)$ with Teixeira mass-fractal structure factor $S(q, D_f, \xi)$:
$$I_{\mathrm{network}}(q) = (\Delta\rho)^2 \phi_{\mathrm{prim}} \int_0^\infty n(R) V(R) P_{\mathrm{sphere}}(q, R)\, dR \cdot S(q, D_f, \xi)$$
$$S(q, D_f, \xi) = 1 + \frac{D_f \Gamma(D_f - 1)}{\left[1 + \frac{1}{(q\xi)^2}\right]^{(D_f - 1)/2}} \frac{\sin\left[(D_f - 1)\arctan(q\xi)\right]}{(q R_{\mathrm{prim}})^{D_f}}$$

### 2. Discrete ACC nanospheres ($I_{\mathrm{ACC}}$)
Observed exclusively in the sub-spinodal regime ($10\ \mathrm{mM}$, ND-i-003) where dynamic network percolation is absent:
$$I_{\mathrm{ACC}}(q) = (\Delta\rho)^2 \phi_{\mathrm{ACC}} \int_0^\infty n_{\mathrm{ACC}}(R) V(R) P_{\mathrm{sphere}}(q, R)\, dR$$

### 3. Prenucleation cluster / fluctuation mesh ($I_{\mathrm{fluc}}$)
Dynamic concentration fluctuations within the semi-dilute solvent mesh, modelled using the generalized Fisher-Burford / Tanaka formulation ($\mu_{\mathrm{fluc}} = 4.0$):
$$I_{\mathrm{fluc}}(q) = \frac{A_{\mathrm{fluc}}}{\left[1 + (q \xi_{\mathrm{fluc}})^2\right]^{\mu_{\mathrm{fluc}} / 2}}$$
where $\xi_{\mathrm{fluc}} \approx 2.7 - 3.1\ \mathrm{nm}$ corresponds to the characteristic correlation length of the prenucleation coordination mesh.

### 4. Secondary macroscopic droplets / Porod tail ($I_{\mathrm{droplet}}$)
Ultra-low-$q$ power law scattering originating from coarsening phase-separated liquid-condensed droplets ($\langle R_{\mathrm{large}} \rangle \approx 2.2 - 3.9\ \mu\mathrm{m}$):
$$I_{\mathrm{droplet}}(q) = \frac{K_{\mathrm{Porod}}}{q^4} \quad \text{or} \quad (\Delta\rho)^2 \phi_{\mathrm{large}} \langle V_{\mathrm{large}} \rangle P_{\mathrm{sphere}}(q, R_{\mathrm{large}})$$

---

## Standardised component palette and visual register

The following table summarises the model components palette and visual register:

| Model component | Hex code | Colour sample | Stroke style | Technical representation |
| :--- | :--- | :--- | :--- | :--- |
| **Total unified fit** | `#dc2626` | ![#dc2626](https://via.placeholder.com/14/dc2626/000000?text=+) | Solid (`-`, $2.2\ \mathrm{pt}$) | Complete multi-component model sum |
| **Spinodal fractal network** | `#059669` | ![#059669](https://via.placeholder.com/14/059669/000000?text=+) | Dashed (`--`, $1.8\ \mathrm{pt}$) | Teixeira mass fractal building block aggregate |
| **Discrete ACC spheres** | `#0284c7` | ![#0284c7](https://via.placeholder.com/14/0284c7/000000?text=+) | Dash-dot (`-.`, $1.8\ \mathrm{pt}$) | Polydisperse sub-spinodal ACC nanospheres |
| **Tanaka PNC mesh** | `#9333ea` | ![#9333ea](https://via.placeholder.com/14/9333ea/000000?text=+) | Dotted (`:`, $1.8\ \mathrm{pt}$) | Viscoelastic dynamic concentration fluctuations |
| **Macroscopic droplets / Porod** | `#d97706` | ![#d97706](https://via.placeholder.com/14/d97706/000000?text=+) | Long-dash (`--`, $1.5\ \mathrm{pt}$) | Secondary liquid-condensed emulsion droplets |
| **Solvent background** | `#64748b` | ![#64748b](https://via.placeholder.com/14/64748b/000000?text=+) | Dotted (`:`, $1.0\ \mathrm{pt}$) | Invariant flat residual background |
| **Experimental data** | `#334155` | ![#334155](https://via.placeholder.com/14/334155/000000?text=+) | Markers (`o`, $\alpha = 0.55$) | Calibrated synchrotron scattering profiles |

---

## Master dashboards and publication diagrams

Interactive tools and SI assets are precompiled and directly accessible:

1. **Master synthesis dashboard**:
   - Web application: [`results/index.html`](results/index.html) and [`results/dashboard/index.html`](results/dashboard/index.html)
   - Features dataset hub cards, physicochemical continuum timelines, live cross-series kinetic trajectory charts, comparative parameter matrices, composite curve galleries, and mathematical formalisms.
2. **Model transformation SVG vector diagram**:
   - File: [`results/output/model_transformation_diagram.svg`](results/output/model_transformation_diagram.svg)
   - Scalable vector schema depicting the phase diagram, model equations, component proportion tallies (including 26% PNC mesh for 40 mM benchmark), and morphological transitions.
3. **Cross-series trajectory synthesis**:
   - File: [`results/output/cross_series_trajectories_comparison.png`](results/output/cross_series_trajectories_comparison.png)
   - 4-panel kinetic comparison displaying fractal dimension $D_f(t)$, correlation length $\xi(t)$, volume and invariant component tallies, and mean squared residuals.

---

## Repository directory structure

```text
project - CaCO3/
├── README.md                                          # Comprehensive project and beamline documentation
├── .gitignore                                         # Git ignore rules for caches and auxiliary files
├── model_description_and_interpretation.md           # Core model formulation and physical interpretation
├── LICENSE                                            # Standard MIT License
├── project notes/                                     # Beamtime logbook, proposal, and reference notes
│   ├── Logbook Beamtime 11-11-2023.docx               # Primary beamtime operational logbook
│   ├── Mathematical Framework for VPS Scattering.docx # Theoretical framework of viscoelastic phase separation
│   ├── proposal ESRF/CH-6828.pdf                      # Original ESRF beamtime proposal
│   ├── refs/                                          # Historical literature references
│   ├── dataNDi009/                                    # Primary short/long-distance frame folders
│   └── drawing.svg, overlay.pdf, ND-009.opj           # Working diagrams and project files
└── results/                                           # Centralised analysis suite
    ├── saxs_colors.py                                 # Invariant 9-component color palette & styling
    ├── build_master_dashboard.py                      # Master interactive HTML synthesis generator
    ├── generate_model_transformation_diagram.py       # SVG vector transformation diagram generator
    ├── generate_cross_series_trajectories.py          # 4-panel kinetic synthesis trajectory builder
    ├── standardize_composite_plots.py                 # Uniform baseline curve re-rendering script
    ├── verify_vps_framework.py                        # Theoretical self-consistency validator
    ├── index.html                                     # Master interactive synthesis dashboard
    ├── dashboard/index.html                           # Standalone synthesis dashboard copy
    ├── output/                                        # Master publication figures and diagrams
    │   ├── model_transformation_diagram.svg
    │   ├── cross_series_trajectories_comparison.png
    │   └── cross_series_trajectories_comparison.json
    ├── ND-i-001_40mM_pH11.3/                          # Rapid quench 40 mM spinodal dataset
    ├── ND-i-002_20mM_pH11.3/                          # 20 mM spinodal gelation dataset
    ├── ND-i-003_10mM_pH11.3/                          # 10 mM sub-spinodal nucleation dataset
    ├── ND-i-004_15mM_pH11.3/                          # 15 mM crossover boundary dataset
    ├── ND-i-005_40mM_pH11.3_CoCl2/                    # 40 mM additive perturbation dataset (+Co2+)
    ├── ND-i-006_40mM_pH9.75/                          # 40 mM buffered pH 9.75 dataset
    └── ND-i-009_40mM_pH10.0/                          # 40 mM buffered pH 10.0 primary benchmark
```

---

## Reproducibility and execution guide

### Environment prerequisites
The analysis suite requires Python 3.10+ with standard scientific libraries:
```bash
pip install numpy scipy matplotlib pandas lmfit
```

### Running the end-to-end analysis
All scripts can be executed using the system Python interpreter:

1. **Re-render the standardised composite baseline fits**:
   ```bash
   python results/standardize_composite_plots.py
   ```

2. **Regenerate the cross-series trajectory synthesis**:
   ```bash
   python results/generate_cross_series_trajectories.py
   ```

3. **Rebuild the vector model transformation diagram**:
   ```bash
   python results/generate_model_transformation_diagram.py
   ```

4. **Verify the viscoelastic phase separation consistency**:
   ```bash
   python results/verify_vps_framework.py
   ```

5. **Recompile the master synthesis dashboards**:
   ```bash
   python results/build_master_dashboard.py
   ```

---

## Data processing and citation
- **Project leader / responsible**: Stephan E. Wolf
- **Affiliation**: Friedrich-Alexander University Erlangen-Nürnberg (FAU) 
- **Contact**: `stephan.e.wolf@fau.de`
- **Funding**: DFG WO1712/7-1 and 8-1 (Heisenberg Project)
- **Data processing responsible**: Tomasz M. Stawski
- **Affiliation**: Federal Institute for Materials Research and Testing (BAM)
- **Contact**: `tomasz.stawski@bam.de`
- **Beamtime**: ESRF - The European Synchrotron (Proposal CH-6828)
- **Primary Beamline Reference**:
  > Narayanan, T., Sztucki, M., Van Vaerenbergh, P., Léonardon, J., Gorini, J., Claustre, L., Sever, F., Morse, J., & Boesecke, P. (2018). *The ID02 high-brilliance small-angle X-ray scattering beamline at the ESRF*. **Journal of Applied Crystallography**, 51(6), 1511–1524.

---

## License

This project is licensed under the standard [MIT License](LICENSE).

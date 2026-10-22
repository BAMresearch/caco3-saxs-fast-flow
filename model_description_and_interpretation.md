# Small-Angle X-ray Scattering (SAXS) Model Description & Physical Interpretation

**Title:** Dual-Morphology Spinodal Demixing in Non-Classical Calcium Carbonate Precipitation: Coexistence of an Extended Viscoelastic Fractal Mesh and Capillary-Pinched Dense Liquid Droplets  
**Target Section:** Theoretical Model Formulation, Kinetic Regimes & Physical Interpretation  
**Data Processing Responsible:** Tomasz Stawski (*Federal Institute for Materials Research and Testing (BAM), Berlin, Germany*)  
**Contact:** `tomasz.stawski@bam.de`  
**License:** Standard MIT License ([LICENSE](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/LICENSE))  

**Associated Project Artifacts:**
* Interactive Master Dashboard: [results/index.html](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/index.html) and [results/dashboard/index.html](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/dashboard/index.html)
* Vector Model Transformation Diagram: [results/output/model_transformation_diagram.svg](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/output/model_transformation_diagram.svg)
* Viscoelastic Phase Separation (VPS) Verification Report: [results/output/vps_framework_verification_report.md](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/output/vps_framework_verification_report.md)
* Cross-Series Kinetic Trajectories: [results/output/cross_series_trajectories_comparison.png](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/output/cross_series_trajectories_comparison.png)
* Trajectory Telemetry Data: [results/output/cross_series_trajectories.json](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/output/cross_series_trajectories.json)
* BibTeX Database: [references.bib](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/references.bib)

---

## 1. Introduction, Physical Premise & Supersaturation Regimes

In highly supersaturated, alkaline aqueous solutions ($\text{Ca}^{2+} + \text{CO}_3^{2-} \rightarrow \text{CaCO}_3$), the initial phase separation trajectory frequently bypasses classical monomer-by-monomer nucleation [De Yoreo et al., 2015; Gebauer et al., 2008; Stawski et al., 2016, 2019]. Instead, precipitation initiates through liquid-liquid phase separation (LLPS) [Gower & Odom, 2000; Wallace et al., 2013; Bewernitz et al., 2012]. Deep within the spinodal instability boundary, barrierless demixing drives the spontaneous condensation of hydrated ionic polymers or Dynamically Ordered Liquid-Like Oxyanion Polymers (DOLLOPs) [Cahn, 1965; Demichelis et al., 2011].

Under stopped-flow rapid-mixing conditions at the ESRF ID02 beamline ($q \in [0.0020, 2.54]\ \text{nm}^{-1}$, real-space dimensions $2.5\ \text{nm}$ to $3.1\ \mu\text{m}$, time resolution $0.045\ \text{s}$ to $10.88\ \text{s}$), the reaction pathways exhibit three distinct supersaturation-dependent regimes:

1. **Low Supersaturation ($10\ \text{mM}\ \text{CaCO}_3$, Series ND-i-003): Classical Barrier-Limited Nucleation**  
   At low supersaturation, precipitation follows a typical nucleation and growth route. The very initial time frames contain predominantly background noise, reflecting an induction period during which a critical thermodynamic barrier must be overcome. The subsequently emerging scattering curves represent discrete, compact spheres of limited polydispersity, describable by a spherical form factor with finite mean radii $\langle R \rangle \approx 200 - 300\ \text{nm}$ at long camera distance ($30.7\ \text{m}$), with minimal intensity at short distance ($3.0\ \text{m}$). This signature corresponds to classical colloidal amorphous calcium carbonate (ACC) nanoparticles suspended in solution.

2. **Intermediate Supersaturation ($15\ \text{mM}\ \text{CaCO}_3$, Series ND-i-004): Hybrid Nucleation-Spinodal Demixing**  
   As concentration increases to $15\ \text{mM}$, the process becomes a hybrid scenario. The system still exhibits compact polydisperse spheres, but their dimensions are substantially smaller ($\langle R \rangle < 100\ \text{nm}$). Simultaneously, an evolving mass fractal network appears, evidenced by a characteristic power-law scattering dependence $I(q) \propto q^{-D}$ with $D \approx 1.8$. Here, barrier-crossing nucleation occurs concurrently with the initial stages of interconnected phase separation.

3. **High Supersaturation ($20\ \text{mM}$ to $40\ \text{mM}\ \text{CaCO}_3$, Series ND-i-001, ND-i-002, ND-i-005, ND-i-006, ND-i-009): Energy-Barrier-Free Spinodal Gelation and Viscoelastic Droplet Pinch-Off**  
   At and above $20\ \text{mM}$, the system reaches an energy-barrier-free composition within the spinodal boundary. Nucleation barriers are absent. Invariable primary building units ($r_0 \approx 8 - 9\ \text{nm}$) appear practically instantaneously without induction delay, establishing a characteristic knee-bending feature in scattering at $q \approx 0.1\ \text{nm}^{-1}$. Kinetic evolution is restricted to lower $q$, where the mass fractal network ($D_m \approx 1.84 - 1.95$) expands. Concurrently, at ultra-low $q$ ($q < 0.005\ \text{nm}^{-1}$), steep power-law decay ($I \propto q^{-4}$) with deep Bessel oscillations reveals giant spherical domains ($\langle R \rangle \approx 2.9 - 3.5\ \mu\text{m}$).

This analytical model proves that these micrometer spheres are **not** solid ACC colloids or crystalline spherulites. Instead, they represent **dense liquid-like precursor droplets (emulsion domains)** that spontaneously pinch off from the continuous spinodal network via **Rayleigh-Plateau capillary instabilities** [Rayleigh, 1878; Plateau, 1873] and **Viscoelastic Phase Separation (VPS)** [Tanaka, 2000, 2001], coexisting in dynamic mechanical equilibrium with the spanning gel mesh.

---

## 2. Multi-Scale Scattering Model Formulation

The total coherent scattering intensity $I_{\text{tot}}(q)$ across the merged reciprocal space ($q \in [0.0020, 2.54]\ \text{nm}^{-1}$, ESRF ID02 [Narayanan et al., 2018]) is expressed as an incoherent superposition of three spatially decoupled morphological contributions plus a flat solvent background:

$$\mathbf{I_{\text{tot}}(q) = I_{\text{droplet}}(q) + I_{\text{mesh}}(q) + I_{\text{high-q}}(q) + I_{\text{bkg}}}$$

```
==================================================================================================
                                    THE THREE-TIER SCATTERING DECOMPOSITION
==================================================================================================

  log I(q) ^
           |  [TIER 1: CAPILLARY DROPLETS / POROD ENVELOPE]
           |  I ~ q^-4 (Porod, smooth 3D interface D = 3)
           |  Bessel Oscillations (R ≈ 3.2 µm, volume fraction ≈ 8%)
           |   \
           |    \     [TIER 2: CONTINUOUS SPINODAL MESH]
           |     \    I ~ q^-1.84 (Percolated Mass Fractal Gel, Dm ≈ 1.84 - 1.95)
           |      \   Primary Nodes (r0 ≈ 8.8 nm, volume fraction ≈ 92%)
           |       \   \
           |        \___\_________
           |                      \     [TIER 3: INTERNAL FLUCTUATIONS / PNCs]
           |                       \    Sub-nanometre Clusters (R_PNC ≈ 0.5 nm)
           |                        \   Ornstein-Zernike / Fisher-Burford (ξ_fluc ≈ 2.0 - 4.6 nm)
           |                         \____________
           |                                      \_____ Flat Incoherent Buffer Background
           +---------------------------------------------------------------------------->
           10^-3                          10^-2                          10^-1            10^0   q (nm^-1)
==================================================================================================
```

### 2.1. Tier 1: Giant Dense-Liquid Spherical Droplets ($q < 0.005\ \text{nm}^{-1}$)

The ultra-low-$q$ scattering regime ($q \in [0.002, 0.005]\ \text{nm}^{-1}$) is governed by isolated dense liquid droplets suspended in the dilute aqueous phase. For a population of spherical droplets with lognormal size distribution $f_N(R)$ defined by median radius $R_0$ and geometric standard deviation $\sigma_{\text{geo}}$, the scattering intensity is:

$$I_{\text{droplet}}(q) = n_{\text{drop}} (\Delta \rho_{\text{drop}})^2 \int_0^\infty V(R)^2 P_{\text{sphere}}(q, R) f_N(R) \, dR$$

where:
* $V(R) = \frac{4}{3}\pi R^3$ is the individual droplet volume.
* $\Delta \rho_{\text{drop}} = \rho_{\text{drop}} - \rho_{\text{water}}$ is the excess electron scattering length density (SLD) of the dense liquid phase relative to water ($(\Delta \rho)^2 \approx 26.88\ \text{cm}^{-2}$ calibrated for amorphous calcium carbonate).
* $P_{\text{sphere}}(q, R)$ is the Rayleigh sphere form factor:
  $$P_{\text{sphere}}(q, R) = \left[ 3 \frac{\sin(qR) - qR\cos(qR)}{(qR)^3} \right]^2 = \left[ \frac{3 j_1(qR)}{qR} \right]^2$$
* In terms of total droplet volume fraction $\phi_{\text{drop}} = n_{\text{drop}} \langle V \rangle$:
  $$I_{\text{droplet}}(q) = \phi_{\text{drop}} (\Delta \rho_{\text{drop}})^2 \frac{\langle V^2 P_{\text{sphere}}(q) \rangle}{\langle V \rangle}$$

#### The Asymptotic Porod Boundary ($D = 3$)
For $q R \gg 1$, the intra-droplet interference oscillations dampen, and the envelope follows Porod's universal law for three-dimensional compact bodies with sharp, smooth interfaces [Porod, 1951]:

$$I_{\text{droplet}}(q) \xrightarrow{q R \gg 1} 2\pi (\Delta \rho_{\text{drop}})^2 \frac{S_{\text{tot}}}{q^4}$$

where $S_{\text{tot}} = 4\pi \langle R^2 \rangle n_{\text{drop}}$ is the total specific interfacial area. The observed logarithmic slope $\frac{d\ln I}{d\ln q} \approx -4.0$ at $q \in [0.002, 0.004]\ \text{nm}^{-1}$ confirms a smooth, compact 3D interface ($D=3$), which cannot be produced by a mass or surface fractal network.

---

### 2.2. Tier 2: Continuous Spinodal Mesh / Mass Fractal Network ($0.005 \le q \le 0.2\ \text{nm}^{-1}$)

The majority of the precipitated material resides in an extended, bicontinuous spinodal network. This continuous mesh is modelled as a mass fractal assembly of primary hydrated nodes of mean radius $r_0 \approx 8.8\ \text{nm}$ connected with fractal dimension $D_m \approx 1.83 - 1.95$ and an upper structural correlation cutoff length $\xi$ [Teixeira, 1988; Sorensen, 2001]:

$$I_{\text{mesh}}(q) = \phi_{\text{mesh}} (\Delta \rho_{\text{mesh}})^2 V_{\text{node}} P_{\text{node}}(q, r_0) S_{\text{fractal}}(q, D_m, \xi)$$

where $P_{\text{node}}(q, r_0)$ is the form factor of the primary building blocks evaluated with lognormal polydispersity $\sigma_0 \approx 0.25 - 0.35$, and $S_{\text{fractal}}(q)$ is the Teixeira-Sorensen fractal structure factor:

$$S_{\text{fractal}}(q) = 1 + \frac{1}{(q r_0)^{D_m}} \frac{D_m \Gamma(D_m - 1)}{\left[ 1 + \frac{1}{(q\xi)^2} \right]^{(D_m - 1)/2}} \sin\left[ (D_m - 1) \arctan(q\xi) \right]$$

#### Limiting Regimes of the Continuous Mesh
* **Fractal Power-Law Regime ($1/\xi \ll q \ll 1/r_0$):**
  $$S_{\text{fractal}}(q) \approx \frac{D_m \Gamma(D_m - 1) \sin\left[ \frac{\pi}{2}(D_m - 1) \right]}{(q r_0)^{D_m}} \propto q^{-D_m} \implies I_{\text{mesh}}(q) \propto q^{-D_m} \approx q^{-1.84}$$
  This exponent reflects universal 3D Diffusion-Limited Cluster Aggregation (DLCA) and percolation of the spinodal liquid channels.
* **Infinite Network Asymptotic Limit ($\xi \to \infty$ or $\xi \ge 10^6\ \text{nm}$):**  
  When the percolated network spans the macroscopic sample volume, the correlation cutoff diverges ($\xi \to \infty$), simplifying the structure factor to its exact analytical limit:
  $$S_{\text{fractal}}(q) = 1 + D_m \Gamma(D_m - 1) \sin\left[ \frac{\pi}{2}(D_m - 1) \right] (q r_0)^{-D_m}$$
* **Guinier Node Regime ($q \sim 1/r_0$):**  
  Crossover into the intra-node form factor, revealing the characteristic node radius of gyration $R_{g,0} = \sqrt{3/5} r_0 \approx 6.8\ \text{nm}$.

---

### 2.3. Tier 3: Prenucleation Clusters (PNCs) & Viscoelastic Network Fluctuations ($q > 0.2\ \text{nm}^{-1}$)

At high scattering vectors, two mathematically complementary formulations describe the elevated scattering plateau:

#### Option A: Prenucleation Clusters / DOLLOPs
Sub-nanometre hydrated ion associates (Dynamically Ordered Liquid-Like Oxyanion Polymers, DOLLOPs [Gebauer et al., 2008; Demichelis et al., 2011]) existing in dynamic equilibrium within the solvent channels:

$$I_{\text{PNC}}(q) = \phi_{\text{PNC}} (\Delta \rho_{\text{PNC}})^2 V_{\text{PNC}} \left[ \frac{3 j_1(q R_{\text{PNC}})}{q R_{\text{PNC}}} \right]^2$$

Calibrated parameters: $\langle R_{\text{PNC}} \rangle = 0.53 \pm 0.15\ \text{nm}$, $\phi_{\text{PNC}} \approx 1.8 \times 10^{-3}$.

#### Option B: Viscoelastic Dynamic Mesh Fluctuations (Tanaka / Fisher-Burford)
Internal concentration fluctuations of the polymer-like dense viscoelastic mesh, with correlation mesh size $\xi_{\text{fluc}}$ and Porod surface-decay exponent $\mu = 4.0$ [Fisher & Burford, 1967; Tanaka, 2000]:

$$I_{\text{fluc}}(q) = \frac{A_{\text{fluc}}}{\left[ 1 + \frac{2}{3\mu} (q \xi_{\text{fluc}})^2 \right]^{\mu / 2}} = \frac{A_{\text{fluc}}}{\left[ 1 + \frac{1}{6} (q \xi_{\text{fluc}})^2 \right]^2}$$

For $\mu = 2$, this expression reduces to the classical Ornstein-Zernike (OZ) concentration fluctuation model $I_{\text{OZ}}(q) = \frac{A_{\text{OZ}}}{1 + (q \xi_{\text{OZ}})^2}$. Across the reaction series, $\xi_{\text{fluc}} \approx 1.5 - 3.5\ \text{nm}$, reflecting the transient mesh size between physical crosslinks of the continuous gel network.

---

## 3. Physical Mechanisms: Origin of Large Spheres in a Spinodal System

```
==================================================================================================
                 PHYSICAL MECHANISMS OF DROPLET PINCH-OFF AND COEXISTENCE
==================================================================================================

  (A) SPINODAL DEMIXING          (B) CAPILLARY INSTABILITY        (C) DROPLET COEXISTENCE
     (Early Bicontinuous Mesh)      (Rayleigh-Plateau Waves)         (Arrested Mesh + Emulsion)

      ...   ...   ...                  ...     ...                    ...       ...
     :   :::   :::   :                :   :::::   :                  :   :::::::   :
     :   :::   :::   :   ======>     (   )     (   )     ======>     :   :       :   :
      ...   ...   ...                 ~~~       ~~~                   ...         ...
     Dense Liquid Struts             Wavelength λ > 2π R_strut       Percolating Gel Network
     Rich in Ca-CO3-H2O              Spontaneous Pinch-Off            (Dm ≈ 1.84, ~92% vol)
                                                                             +
                                                                     ( O )       ( O )
                                                                    Micrometre Dense Droplets
                                                                     (R ≈ 3 µm, ~8% vol)
==================================================================================================
```

### 3.1. Rayleigh-Plateau Capillary Instability of Bicontinuous Struts

1. **Capillary Dispersion in Viscous Liquid Cylinders:**  
   In early spinodal demixing, Cahn-Hilliard decomposition produces an isotropic, interconnected bicontinuous structure of dense liquid cylinders embedded in a lower-density aqueous phase [Cahn, 1965]. Because the dense phase has an interfacial tension $\gamma \sim 1 - 10\ \text{mN/m}$ against water, cylindrical ligaments of radius $R_{\text{strut}}$ are unstable to capillary radius fluctuations [Plateau, 1873; Rayleigh, 1878]:
   $$r(z) = R_{\text{strut}} + \epsilon_0 \cos(k z) e^{\omega t}$$

2. **The Cutoff Condition:**  
   Surface area per unit length decreases if the perturbation wavelength $\lambda = 2\pi / k$ exceeds the circumference:
   $$\lambda > 2\pi R_{\text{strut}} \quad \iff \quad k R_{\text{strut}} < 1$$

3. **The Dominant Pinch-Off Mode:**  
   The hydrodynamic growth rate $\omega(k)$ peaks at:
   $$k_{\text{opt}} R_{\text{strut}} \approx 0.697 \implies \lambda_{\text{opt}} \approx 9.02 R_{\text{strut}}$$

4. **Resulting Spherical Droplet Radius:**  
   Volume conservation during capillary collapse yields:
   $$\frac{4}{3}\pi R_{\text{drop}}^3 = \pi R_{\text{strut}}^2 \lambda_{\text{opt}} \implies R_{\text{drop}} \approx 1.89 R_{\text{strut}}$$
   Local variations in strut thickness during mixing generate ligaments ranging up to $R_{\text{strut}} \sim 1.5 - 1.7\ \mu\text{m}$, directly yielding pinched-off spherical droplets of $R_{\text{drop}} \approx 2.8 - 3.2\ \mu\text{m}$, in quantitative agreement with the experimental fit ($\langle R_{\text{drop}} \rangle = 3.21 \pm 0.04\ \mu\text{m}$).

---

### 3.2. Viscoelastic Phase Separation (VPS) and Tanaka Two-Fluid Dynamics

Because the dense mineral network exhibits a stress relaxation time $\tau_{\text{mesh}}$ that is orders of magnitude longer than the solvent diffusion time $\tau_{\text{solv}}$ ($\tau_{\text{mesh}} \gg \tau_{\text{solv}}$), phase separation is governed by **Viscoelastic Phase Separation (VPS)** [Tanaka, 2000, 2001]:

1. **Two-Fluid Force Balance & Dynamic Asymmetry:**  
   In the Two-Fluid Model, the velocities of the mineral network ($\mathbf{v}_p$) and solvent ($\mathbf{v}_s$) are distinct. Under incompressibility, the force balance is:
   $$-\phi \nabla \mu + \nabla \cdot \boldsymbol{\sigma}^{(p)} - \zeta (\mathbf{v}_p - \mathbf{v}_s) = 0$$
   where $\mu$ is the chemical potential, $\boldsymbol{\sigma}^{(p)}$ is the polymer viscoelastic stress tensor, and $\zeta$ is the friction coefficient. In our system, dynamic asymmetry is maximal ($\alpha \approx 1$): the mineral network supports all viscoelastic stress, while the solvent supports negligible stress.

2. **The Viscoelastic Length Scale ($\xi_{\text{ve}}$):**  
   The competition between stress relaxation ($\tau_{\text{mesh}}$) and diffusion ($D_{\text{coop}} = K_{\text{os}} / \zeta$) defines a characteristic dynamic length:
   $$\xi_{\text{ve}} = \sqrt{\frac{G \tau_{\text{mesh}}}{\zeta}} = \sqrt{\frac{\eta_{\text{ve}}}{\zeta}}$$
   * On scales $q \xi_{\text{ve}} \gg 1$ (small scales), the network cannot relax within diffusion timescales and behaves as an elastic solid, suppressing further composition fluctuations.
   * On scales $q \xi_{\text{ve}} \ll 1$ (large scales), stress relaxes quickly and the system behaves as a simple demixing fluid.

3. **Bulk Stress and Volume Shrinking:**  
   The decomposition of the stress tensor into shear ($\boldsymbol{\sigma}_s$) and bulk ($\sigma_{\text{bulk}}$) components reveals that the isotropic contraction of the network is resisted by its bulk elastic modulus:
   $$\sigma_{\text{bulk}} = K_{\text{os}} \nabla \cdot \mathbf{v}_p$$
   This bulk stress drives the **volume shrinking** observed in scattering: the mineral gel contracts like a sponge, squeezing out solvent until elastic tension balances osmotic pressure.

4. **Stress Relaxation via Rupture and Droplet Shedding:**  
   When coarsening strains exceed the transient yield stress of the struts, the network ruptures mechanically. Detached fluid domains contract under interfacial tension into compact spherical droplets ($\langle R_{\text{drop}} \rangle \approx 3.2\ \mu\text{m}$), while unbroken connected paths arrest into the spanning fractal gel mesh ($D_m \approx 1.84$).

---

### 3.3. In Situ WAXS Proof: Absence of Bragg Diffraction

To experimentally verify whether the large spheres could be crystalline vaterite or calcite, in situ Wide-Angle X-ray Scattering (WAXS) was acquired simultaneously with SAXS during the entire reaction trajectory ($t = 0 - 7.0\ \text{s}$).
* **Diffraction Signatures:** Throughout the observation window, the WAXS patterns show **no Bragg reflections** associated with vaterite [e.g. (110) at $d = 3.58\ \text{Å}$, (112) at $d = 3.29\ \text{Å}$, (114) at $d = 2.73\ \text{Å}$] or calcite [e.g. (104) at $d = 3.03\ \text{Å}$].
* **Amorphous/Liquid Halo:** The signal displays exclusively the broad, diffuse halo of bulk water and hydrated Ca-carbonate ion associations centered near $q \approx 20\ \text{nm}^{-1}$.
* **Definitive Conclusion:** The large spheres are **purely liquid-like precursor droplets (PILP domains)** rather than solid crystalline spherulites.

---

## 4. The Volume-Scattering Disparity: Why a Trace Population Dominates Low $q$

The forward scattering intensity per object scales with $V^2$, or $R^6$:
$$I_{\text{single}}(0) = (\Delta\rho)^2 V^2 \propto R^6$$

Comparing a single micrometre droplet ($R = 3210\ \text{nm}$) to a primary network node ($r_0 = 8.8\ \text{nm}$):
$$\frac{V_{\text{drop}}}{V_{\text{node}}} = \left(\frac{3210}{8.8}\right)^3 \approx 4.85 \times 10^7$$
$$\frac{I_{\text{drop}}(0)}{I_{\text{node}}(0)} = \left(\frac{V_{\text{drop}}}{V_{\text{node}}}\right)^2 \approx \mathbf{2.35 \times 10^{15}}$$

Macroscopically, for droplet volume fraction $\phi_{\text{drop}} \approx 7.0 \times 10^{-4}$ ($\approx 8.2\%$ of total mineral mass) and mesh volume fraction $\phi_{\text{mesh}} \approx 7.8 \times 10^{-3}$ ($\approx 91.8\%$ of total mineral mass):
$$\frac{I_{\text{drop}}(0)}{I_{\text{mesh}}(0)} = \frac{\phi_{\text{drop}} V_{\text{drop}}}{\phi_{\text{mesh}} V_{\text{cluster}}} \approx \mathbf{58.9}$$

```
==================================================================================================
                        VOLUME FRACTION VS. SCATTERING POWER BUDGET
==================================================================================================

  PRECIPITATED VOLUME FRACTION                   FORWARD SCATTERING INTENSITY AT q -> 0
  
  [███████████████████████████████░░░]           [░░░░░░░░░░░░░░░░░░░░░░░░░░████████]
  
  ■ Continuous Mesh:  91.8%                      ■ Continuous Mesh:   1.7%
  ░ Capillary Droplets: 8.2%                     █ Capillary Droplets: 98.3%
==================================================================================================
```

At $q_{\text{min}} = 0.00213\ \text{nm}^{-1}$, the droplets generate **$77.8\%$** of measured intensity. Because the droplet scattering decays steeply ($\propto q^{-4}$), its contribution plunges below $2\%$ by $q \approx 0.008\ \text{nm}^{-1}$, allowing the continuous fractal mesh ($I \propto q^{-1.84}$) to dominate cleanly over the rest of the profile.

---

## 5. Calibrated Model Parameters (ND-i-009, 40 mM Pure $\text{CaCO}_3$, pH 10.0 Benchmark)

| Model Component | Parameter Symbol | Fitted Value | Std. Error | Physical Significance |
| :--- | :---: | :---: | :---: | :--- |
| **Capillary Droplets** | $\langle R_{\text{drop}} \rangle$ | $3.21\ \mu\text{m}$ | $\pm 0.04\ \mu\text{m}$ | Mean droplet radius from capillary break-up |
| | $\sigma_{\text{drop}}$ | $0.082$ | $\pm 0.006$ | Polydispersity showing sharp, narrow size mode |
| | $\phi_{\text{drop}}$ | $7.02 \times 10^{-4}$ | $\pm 0.15 \times 10^{-4}$ | Volume fraction of pinched-off emulsion phase |
| **Continuous Mesh** | $D_m$ | $1.925$ | $\pm 0.002$ | Fractal dimension of DLCA spinodal network |
| | $\langle r_0 \rangle$ | $9.20\ \text{nm}$ | $\pm 0.07\ \text{nm}$ | Primary network node radius |
| | $\sigma_0$ | $0.312$ | $\pm 0.008$ | Node size distribution breadth |
| | $\phi_{\text{mesh}}$ | $7.81 \times 10^{-3}$ | $\pm 0.12 \times 10^{-3}$ | Volume fraction of spanning fractal gel mesh |
| **PNCs / Fluctuations** | $\langle R_{\text{PNC}} \rangle$ | $0.528\ \text{nm}$ | $\pm 0.147\ \text{nm}$ | Radius of sub-nanometre ion clusters (DOLLOPs) |
| | $\phi_{\text{PNC}}$ | $1.84 \times 10^{-3}$ | $\pm 0.22 \times 10^{-3}$ | Precursor ion clusters in solution channels |
| | $\xi_{\text{fluc}}$ | $4.62\ \text{nm}$ | $\pm 0.64\ \text{nm}$ | Transient mesh correlation length ($\approx \frac{1}{2} r_0$) |
| **Background** | $bkg$ | $0.0216\ \text{cm}^{-1}$ | $\pm 0.0003\ \text{cm}^{-1}$ | Flat water / buffer incoherent scattering |
| **Statistical Quality** | $\chi^2_{\text{red}}$ (Overall) | **$1.958$** | : | Superior across all competitive models |
| | $\chi^2_{\text{red}}$ ($q > 0.3$) | **$0.542$** | : | Resolves high-$q$ plateau without overfitting |
| | AIC | **$1202.8$** | : | Information criterion winner ($\Delta\text{AIC} > -28$) |

---

## 6. Comparison with Alternative Physical Hypotheses

| Feature / Metric | Hypothesis 1: Crystalline Spherulites (Vaterite/Calcite) | Hypothesis 2: Solid ACC Colloid Nucleation | Hypothesis 3 (This Work): Spinodal Bicontinuous Mesh + Capillary Droplets |
| :--- | :--- | :--- | :--- |
| **Physical State of Droplets** | Rigid, faceted polycrystalline aggregates. | Solid, spherical amorphous nanoparticles. | Hydrated, viscous liquid precursor droplets (PILP). |
| **Diffraction Signatures** | Strong Bragg peaks in WAXS. **Not observed.** | Broad amorphous halo, no Bragg peaks. | Broad solvent-like liquid structure factor. Perfectly matches data. |
| **Boundary Topology** | Rough surface ($D_s \approx 2.3 - 2.8$). Yields $q^{-3.2}$ to $q^{-3.7}$. | Spherical interface, but typically smaller ($R \sim 20 - 100\ \text{nm}$). | Surface-tension-smoothed meniscus ($D=3$). Strictly yields $q^{-4}$ Porod decay. |
| **Mesh Relationship** | Incompatible; crystals consume gel. | Unclear how solid spheres form continuous $D_m \approx 1.84$ mesh. | **Direct hydrodynamic daughter phase:** pinched off from the continuous mesh via Rayleigh-Plateau. |
| **Statistical Fit ($\chi^2_{\text{red}}$)** | $> 15.0$ (cannot match low-$q$ Bessel oscillation positions). | $2.77$ (fails high-$q$ plateau and intermediate transition). | **$1.96$** (unmatched accuracy across all $q$ decades). |

---

## 7. Experimental & Beamline Specifications (ESRF ID02)

| Parameter | Specification | Beamline Operation Notes |
| :--- | :--- | :--- |
| **Synchrotron Facility** | **ESRF - The European Synchrotron** (Grenoble, France) | Extremely Brilliant Source (EBS) storage ring (Proposal CH-6828) |
| **Beamline** | **ID02 (TRUSAXS / USAXS / SAXS)** | High-brilliance undulator beamline [Narayanan et al., 2018] |
| **Photon Energy ($E$)** | **$12.23\ \text{keV}$ ($12,230.0\ \text{eV}$)** | $\text{Si}(111)$ double crystal monochromator (`egy=12230.0`, `theta=9.308`) |
| **X-ray Wavelength ($\lambda$)** | **$0.101377\ \text{nm}$ ($1.0138\ \text{Å}$)** | Calibrated pyFAI wavelength parameter (`WaveLength = 1.01377e-10 m`) |
| **Detector System** | **Dectris Eiger2 4M** | Hybrid photon counting pixel detector ($2068 \times 2162\ \text{pixels}$, $75\ \mu\text{m}$ pitch) |
| **Short-Distance Camera (SD)** | **$3.0\ \text{m}$ nominal ($2.9997\ \text{m}$ calibrated)** | Covers $q \in [0.022, 2.54]\ \text{nm}^{-1}$; $100\ \mu\text{L}$ to $200\ \mu\text{L}$ stopped-flow shots |
| **Long-Distance Camera (LD)** | **$31.0\ \text{m}$ nominal ($30.7059\ \text{m}$ calibrated)** | Covers $q \in [0.002, 0.25]\ \text{nm}^{-1}$; $100\ \mu\text{L}$ to $200\ \mu\text{L}$ stopped-flow shots |
| **Combined Composite $q$-Range** | **$0.0020\ \text{nm}^{-1}$ to $2.54\ \text{nm}^{-1}$** | $>3$ orders of magnitude in reciprocal space stitched across cameras |
| **Sample Environment** | **BioLogic SFM-400 Stopped-Flow Mixer** | 4-syringe pneumatic injection ($3.33\ \text{mL/s}$/syringe) with $1.3\ \text{mm}$ capillary |

---

## 8. Complete Reference List

1. **Bewernitz, M. A., Ersen, O., Chakarova, R., & Gower, L. B.** (2012). A metastable liquid precursor phase of calcium carbonate and its role in form development. *Faraday Discussions*, 159, 291–312.
2. **Cahn, J. W.** (1965). Phase separation by spinodal decomposition in isotropic systems. *The Journal of Chemical Physics*, 42(1), 93–99.
3. **De Yoreo, J. J., Gilbert, P. U. P. A., Sommerdijk, N. A. M., Penn, R. L., Whitelam, S., et al.** (2015). Crystallization by particle attachment in synthetic, biogenic, and geologic environments. *Science*, 349(6247), aaa6760.
4. **Demichelis, R., Raiteri, P., Gale, J. D., Quigley, D., & Gebauer, D.** (2011). Stable prenucleation mineral clusters are liquid-like ionic polymers. *Nature Communications*, 2(1), 590.
5. **Fisher, M. E., & Burford, R. J.** (1967). Theory of critical-point scattering and correlation functions. *Physical Review*, 156(2), 583–622.
6. **Gebauer, D., Völkel, A., & Cölfen, H.** (2008). Stable prenucleation calcium carbonate clusters. *Science*, 322(5909), 1819–1822.
7. **Gower, L. B., & Odom, D. J.** (2000). Biomimetic mineralisation via a polymer-induced liquid-precursor (PILP) process. *Journal of Crystal Growth*, 210(4), 719–734.
8. **Gower, L. B.** (2008). The polymer-induced liquid-precursor (PILP) process: a journey from amorphous to crystalline materials. *Chemical Reviews*, 108(11), 4551–4627.
9. **Narayanan, T., Sztucki, M., Van Vaerenbergh, P., Léonardon, J., Gorini, J., et al.** (2018). The ID02 high-brilliance small-angle X-ray scattering beamline at the ESRF. *Journal of Applied Crystallography*, 51(6), 1511–1524.
10. **Plateau, J. A. F.** (1873). *Statique expérimentale et théorique des liquides soumis aux seules forces moléculaires*. Gauthier-Villars, Paris.
11. **Porod, G.** (1951). Die Röntgenkleinwinkelstreuung von dichtgepackten kolloiden Systemen. *Kolloid-Zeitschrift*, 124(2), 83–114.
12. **Rayleigh, Lord.** (1878). On the instability of jets. *Proceedings of the London Mathematical Society*, 1(1), 4–13.
13. **Siggia, E. D.** (1979). Late stages of spinodal decomposition in binary mixtures. *Physical Review A*, 20(2), 595–605.
14. **Sorensen, C. M.** (2001). Light scattering by fractal aggregates: a review. *Aerosol Science and Technology*, 35(2), 648–687.
15. **Stawski, T. M., Van Driessche, A. E. S., Ossorio, M., Rodriguez-Blanco, J. D., Besselink, R., & Benning, L. G.** (2016). Formation of amorphous calcium carbonate: Direct observation of non-classical nucleation pathways. *Nature Communications*, 7, 11177.
16. **Stawski, T. M., & Benning, L. G.** (2019). Nucleation and growth of inorganic nanoparticles: from ions to crystals studied by in situ SAXS/WAXS. *Nanoscale*, 11(30), 14068–14083.
17. **Tanaka, H.** (2000). Viscoelastic phase separation. *Journal of Physics: Condensed Matter*, 12(15), R207–R264.
18. **Tanaka, H.** (2001). Viscoelastic phase separation in soft matter and foods. *Current Opinion in Colloid & Interface Science*, 6(3), 259–271.
19. **Teixeira, J.** (1988). Small-angle scattering by fractal systems. *Journal of Applied Crystallography*, 21(6), 781–785.
20. **Wallace, A. F., Henzler, K., Sun, S., Raiteri, P., Gale, J. D., Gebauer, D., Cölfen, H., & De Yoreo, J. J.** (2013). Microscopic evidence for liquid-liquid separation in supersaturated CaCO3 solutions. *Science*, 341(6148), 885–889.

---

## 9. License & Attribution

* **Data Processing Responsible:** Tomasz Stawski (`tomasz.stawski@bam.de`), Federal Institute for Materials Research and Testing (BAM), Unter den Eichen 87, 12205 Berlin, Germany.
* **Beamtime Proposal:** ESRF Proposal CH-6828 (*Time-resolved Small-Angle X-ray Scattering on Non-Classical Nucleation Pathways of Calcium Carbonate*).
* **License:** This work is licensed under the terms of the standard [MIT License](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/LICENSE).

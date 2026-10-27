# Mathematical Verification of Tanaka Viscoelastic Phase Separation (VPS) Concepts

## 1. Dynamic Asymmetry and the Viscoelastic Length Scale

The Tanaka two-fluid model introduces the dynamic asymmetry parameter:

$$\alpha = \frac{\zeta_p}{\zeta_s} - 1 \gg 0$$

where $\zeta_p$ is the friction coefficient of the polymer-like mineral network (DOLLOPs / hydrated ionic clusters), and $\zeta_s$ is the solvent friction coefficient. Because ionic oligomers relax slowly relative to water, the viscoelastic length scale emerges as:

$$\xi_{\text{ve}} = \sqrt{\frac{K_{\text{os}} + \frac{4}{3}G}{\zeta}} \tau_r \sim \sqrt{D_{\text{coop}} \tau_r}$$

where $K_{\text{os}}$ is the osmotic bulk modulus, $G$ is the shear modulus, $\zeta$ is the inter-fluid friction, and $\tau_r$ is the stress relaxation time.

In our scattering data, the multi-scale hierarchy verifies this crossover:
- **For $q > 1/\xi_{\text{ve}}$ ($r < \xi_{\text{ve}}$):** Fluctuations probe internal viscoelastic gel modes ($\xi_{\text{fluc}} \in [1.5, 4.6]\text{ nm}$). The network resists hydrodynamic deformation, leading to frozen domain size or volume shrinking.
- **For $q < 1/\xi_{\text{ve}}$ ($r > \xi_{\text{ve}}$):** Stress has relaxed, permitting capillary flow and phase inversion into macroscopic liquid droplets ($\langle R_{\text{large}} \rangle \in [2.2, 3.5]\ \mu\text{m}$).

---

## 2. Quantitative Proof of Volume Shrinking (ND-i-009 Benchmark)

In classical Model B/H spinodal decomposition, phase separation proceeds with constant phase volumes. In contrast, Viscoelastic Phase Separation requires **volume shrinking** of the transient network, accompanied by solvent expulsion. In ND-i-009 (351 frames, 15 ms cadence):
- **Primary building block radius is invariant:** $R_{\text{prim}} \approx 8.35 \to 8.69\text{ nm}$ (mean $8.34\text{ nm}$, composite $9.40\text{ nm}$) throughout the entire reaction run, while the entangled PNC fluctuation mesh correlation length drifts smoothly as $\xi_{\text{fluc}} \approx 2.72 \to 2.84\text{ nm}$.
- **Network fractal dimension increases monotonically:** $D_m$ drifts from $1.881$ to $1.925$.
- **Densification mechanism:** The increase in $D_m$ at constant node radius directly represents the compaction and densification of the network strands, confirming Tanaka volume shrinking.


---

## 3. Phase Inversion into Macroscopic Droplets

As the viscoelastic network contracts, elastic stress $\sigma_{ij}$ accumulates along interconnected struts. When elastic tension exceeds the local yield threshold, struts rupture and retract into spherical droplets.

This explains why all spinodal series (ND-i-004, ND-i-002, ND-i-009, ND-i-001, ND-i-005) exhibit a stationary macroscopic droplet population ($\langle R_{\text{large}} \rangle \in [2.2, 3.5]\ \mu\text{m}$) coexisting with the space-filling gel mesh.

Under Rayleigh-Plateau capillary instability:

$$\lambda_{\text{opt}} \approx 9.02 R_{\text{strut}} \implies R_{\text{drop}} \approx 1.89 R_{\text{strut}} \approx 2.8 - 3.5\ \mu\text{m}$$

---

## 4. Distinction from Arrested Spinodal Decomposition (ASD)

Unlike Arrested Spinodal Decomposition, which permanently freezes domain growth via repulsive glass jamming, the calcium carbonate system shows:
1. Continuous structural evolution and progressive volume shrinking without sudden structural freezing.
2. Dynamic stress relaxation through mechanical rupture rather than kinetic trapping.
3. Coexistence of an open percolating gel mesh with liquid macroscopic droplets.

---

## 5. Bridging Nanoscale Scattering Invariants to Macroscopic Rheology

The dynamic asymmetry parameter $\alpha = (\zeta_p / \zeta_s) - 1 \gg 0$ in Tanaka's two-fluid model can be quantitatively evaluated by coupling structural parameters extracted from the SAXS data to viscoelastic scaling relations.

### 5.1 Hierarchical Length Scales and Moduli
The multi-scale scattering model reveals two distinct mesh correlation scales:
1. **Inter-Strand Network Pore Mesh ($\xi_{\text{pore}} \approx 40\text{ nm}$):**
   The open space-filling mass fractal network ($D_m \approx 1.88 - 1.93$) possesses an internal pore correlation length of $\xi_{\text{pore}} \approx 40\text{ nm}$, corresponding to an inorganic phase volume fraction of $\phi \approx 12 - 14\%$ (a swollen gel containing $\sim 87\%$ trapped solvent). The cooperative diffusion coefficient of the solvent through this mesh is:
   $$D_{\text{coop}} = \frac{k_B T}{6\pi \eta_s \xi_{\text{pore}}} \approx \frac{4.11 \times 10^{-21}\text{ J}}{6\pi (0.89 \times 10^{-3}\text{ Pa}\cdot\text{s}) (40 \times 10^{-9}\text{ m})} \approx 6.1 \times 10^{-12}\text{ m}^2/\text{s}$$
   The macroscopic elastic plateau shear modulus of the swollen gel network scales as:
   $$G_{\text{gel}} \approx \frac{k_B T}{\xi_{\text{pore}}^3} \approx \frac{4.11 \times 10^{-21}\text{ J}}{(40 \times 10^{-9}\text{ m})^3} \approx 64\text{ Pa}$$

2. **Intra-Strand PNC Fluctuation Mesh ($\xi_{\text{fluc}} \approx 2.8\text{ nm}$):**
   Within the dense network strands, high-$q$ concentration fluctuations probe the entangled prenucleation cluster (DOLLOP) mesh ($\xi_{\text{fluc}} = 2.7 - 2.8\text{ nm}$). The local mechanical stiffness of individual struts scales as:
   $$G_{\text{strand}} \approx \frac{k_B T}{\xi_{\text{fluc}}^3} \approx \frac{4.11 \times 10^{-21}\text{ J}}{(2.8 \times 10^{-9}\text{ m})^3} \approx 1.87 \times 10^5\text{ Pa} \approx 180\text{ kPa}$$
   This dual-modulus hierarchy explains how individual struts are sufficiently rigid to resist immediate capillary coalescence, while the macroscopic percolating gel remains soft and compliant.

---

### 5.2 Macroscopic Relaxation Time and Effective Viscosity
Because no spinodal correlation peak $q_m$ is detected down to the low-$q$ instrumental resolution limit ($q_{\text{min}} = 0.002\text{ nm}^{-1}$), the characteristic domain spacing and viscoelastic length scale must strictly exceed the real-space resolution boundary:
$$d_{\text{max}} = \frac{2\pi}{q_{\text{min}}} \approx 3.14\ \mu\text{m} \implies \xi_{\text{ve}} > 3.14\ \mu\text{m}$$

Combining $\xi_{\text{ve}} > 3.14\ \mu\text{m}$ with the cooperative diffusion coefficient $D_{\text{coop}} \approx 6.1 \times 10^{-12}\text{ m}^2/\text{s}$ yields a lower bound for the structural stress relaxation time:
$$\tau_r = \frac{\xi_{\text{ve}}^2}{D_{\text{coop}}} > \frac{(3.14 \times 10^{-6}\text{ m})^2}{6.1 \times 10^{-12}\text{ m}^2/\text{s}} \approx 1.6 - 1.8\text{ s}$$

This relaxation time is orders of magnitude longer than the molecular rearrangement time of bulk water ($\sim 10^{-11}\text{ s}$), formally proving extreme dynamic asymmetry. 

The effective zero-shear viscosity of the transient dense phase is:
$$\eta_{\text{eff}} \approx G_{\text{gel}} \tau_r > (64\text{ Pa}) \times (1.8\text{ s}) \approx 115\text{ Pa}\cdot\text{s}$$
This value is more than five orders of magnitude ($> 10^5$) higher than the viscosity of bulk water ($\eta_{\text{water}} \approx 10^{-3}\text{ Pa}\cdot\text{s}$), quantifying the physical basis of kinetic arrest during liquid-liquid spinodal demixing.

---

### 5.3 Absence of a Spinodal Correlation Peak ($q_m < q_{\text{min}}$)
In classical fluid spinodal decomposition (Cahn-Hilliard / Model H), a scattering correlation peak $q_m(t)$ emerges at intermediate $q$ and shifts to lower $q$ as bicontinuous domains coarsen. In our time-resolved SAXS patterns, no such peak is observed.

This absence directly corroborates the VPS mechanism:
1. The dynamic asymmetry slows down domain growth at early times, locking domain periodicity at scales $\Lambda > 3.14\ \mu\text{m}$.
2. The entire accessible SAXS window ($q \in [0.002, 0.4]\text{ nm}^{-1}$) probes the *internal structure* of the sponge-like bicontinuous network, rather than the inter-domain periodicity.
3. Rather than macroscopic hydrodynamic coarsening, stress relaxation proceeds via **volume shrinking** ($D_m$ densification at invariant $R_{\text{prim}}$) and localized strut rupture into coexisting spherical droplets ($\langle R_{\text{large}} \rangle \approx 2.2 - 3.5\ \mu\text{m}$).


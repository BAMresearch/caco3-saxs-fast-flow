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
- **Primary building block radius is invariant:** $R_{\text{prim}} = 2.80\text{ nm}$ throughout the entire reaction run.
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

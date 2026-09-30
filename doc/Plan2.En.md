# U–Zr Fuel Constituent Redistribution under a Temperature Gradient: A Stepwise MOOSE Reproduction and Extension Plan (Revision 2)

> Subject: constituent redistribution of U–Zr metallic fuel under a temperature gradient  
> Numerical framework: MOOSE, phase-field module, split Cahn–Hilliard equation  
> Scope: strictly limited to the physical models explicitly used by Mohanty et al. (2011), Wen et al. (2022), and Jung et al. (2025)  
> Recommended implementation order: one-dimensional radial (axisymmetric) item-by-item verification $\rightarrow$ one-dimensional fully coupled multiphysics $\rightarrow$ two-dimensional $r$–$z$ extension  
> Revision status: Revision 2, reviewed and corrected by the supervisor. The original student version is preserved as `Plan1.En.md`. All changes are marked in the text with `> Revision R#` blocks and summarized in Section 0.

## 0. Supervisor's Revision Log

The student's Revision 1 is well organized and its staged philosophy is sound. The corrections below concern physics that was stated incorrectly, priorities that were inverted, and traps that the plan did not warn about. Each item gives the location, what was wrong, and why it matters.

| ID | Location | What was wrong or missing | Why it matters |
| --- | --- | --- | --- |
| R1 | §1.3 | The unit audit did not warn that Mohanty's $\beta_i$ (Fig. 2) and $\widetilde Q_i^*$ are per-atom quantities, while Wen and Jung multiply $\beta_i$ by a molar volume $V_m$ in $\mathrm{m^3\,mol^{-1}}$. | Mixing per-atom and per-mole bases silently changes every mobility by a factor $N_A\approx6\times10^{23}$; this is the single most likely source of a wrong time scale. |
| R2 | §2.1, §3.1, §9 | Revision 1 treated the printed Mohanty Eq. (3), which has no composition weights inside the bracket, as suspicious, and treated the student's added Darken weights as a plausible "correction". | Eq. (3) as printed is exactly the Jaffe–Shewmon/Darken thermotransport result for a binary substitutional alloy. The added weights have no thermodynamic basis. The README's justification, divergence of $M_Q/M_c$ as $c\to0$, is false: the ratio is finite. |
| R3 | §2.1 | The plan did not discuss how the thermodynamic driving force $\nabla(\partial f/\partial c)$ is split between an "isothermal" part and a $\nabla T$ part, nor that the pure-element end-member terms $G_U^0(T)-G_{\mathrm{Zr}}^0(T)$ contribute a $\nabla T$-proportional flux. | In a temperature gradient $\nabla\mu$ contains $\partial^2f/\partial c\,\partial T\,\nabla T$. Dropping end-member terms in one branch and keeping them in another changes the effective heat of transport by 5–20 %. This must be a documented choice, not an accident. |
| R4 | §2.2, §2.3, §9 | The sign convention of the Zr heat of transport differs between the papers: Mohanty uses $\widetilde Q_{\mathrm{Zr}}^*<0$ with a minus sign in Eq. (3); Wen Eq. (31) and Jung Eq. (18) use a plus sign with positive tabulated magnitudes. Revision 1 did not mention this. | Copying Mohanty's negative $\widetilde Q_{\mathrm{Zr}}^*$ into the Wen/Jung formula reverses the Zr migration direction in the $\gamma$ phase. |
| R5 | §2.2, Stage 4, §9 | Revision 1 treated "should $\xi$ also multiply $M_T$?" purely as a literal-reproduction question. | Physically, a vacancy-mediated mechanism enhances $\beta_i$, hence both $M_c$ and $M_T$. Enhancing only $M_c$ divides the steady-state Soret ratio $M_T/M_c$ by $\xi$, i.e. it *weakens* the final redistribution while accelerating the transient. Wen's $1.3\times10^5\ \mathrm s$ (1.5 day) run is far from steady state, so the paper's Fig. 2/3 comparison does not test this. A `xi_both` branch is required as a physics-motivated sensitivity, and the paper's claim that RED removes the need for artificial diffusivity enlargement must be read against the 1.5-day time scale. |
| R6 | §2.2, §2.3 | The `literal` $\alpha/\beta$ switching branch was given equal standing with the physically consistent one; Jung's Eqs. (3)–(4) were not identified as carrying the same inversion. | Both papers' result figures unambiguously show $\alpha+\gamma$ at the cold surface and $\beta+\gamma$ inside. `phase_consistent` is the primary branch; `literal` is a documented check only. The switching half-width of 2 K also needs a resolution and sensitivity statement. |
| R7 | §2.3, Stage 5, Stage 8, §8.1 | Revision 1 made the slab (Cartesian) 1-D operator the primary geometry for Wen/Jung and demoted the axisymmetric form to a sensitivity study. | The physical object is a cylindrical fuel slug; Wen calls the domain axisymmetric and Jung converts LHGR to $\dot q$ with $\pi r_0^2$. For the same $\dot q$ and $k$ the centerline temperature rise in a slab is exactly twice that of a cylinder. MOOSE `coord_type = RZ` is a one-line change. RZ is the primary geometry. |
| R8 | Stage 5, §7.2 | The plan hesitated between running a transient heat equation without $\rho c_p$ data and a steady one. | The thermal diffusion time of the slug, $r_0^2\rho c_p/k\sim1\ \mathrm s$, is 5–8 orders of magnitude shorter than redistribution times. Quasi-static heat conduction at every time step is the physically correct, parameter-free choice; the transient form adds no physics. |
| R9 | §2.3, Stage 7 | The typographic ambiguity of the sodium correction factor $P_c$ (whether $(1-P)$ is in the numerator or the denominator of the prefactor) was not listed, and no physical sanity limits were given. | $P_c(P_{\mathrm{Na}}=0)=1$ and $P_c\ge1$ whenever $k_{\mathrm{Na}}>k_0$ are necessary checks; at 900 K $k_{\mathrm{Na}}/k_0\approx2$, so an implementation error here changes the sign of the sodium effect. |
| R10 | §2.3, Stage 3, §9 | Jung's $\gamma$-phase free energy evaluates the ideal-mixing and Redlich–Kister terms at a fixed $T_{\mathrm{ref}}=1010\ \mathrm K$, and Jung's interface width $l=0.5\ \mathrm{mm}$ is 23 % of the fuel radius. Neither was stated. | With $T_{\mathrm{ref}}$ fixed, the only explicit $T$ dependence of $f_\gamma$ is in the SGTE end-members, which changes the chemical part of the Soret driving force (see R3). An interface 12 times wider than Wen's $\approx40\ \mathrm{\mu m}$ smears the zone boundaries and changes the required mesh. |
| R11 | Stage 3 | The KKS stage verified the constraints but never checked the resulting phase diagram. | Jung admits over-prediction of surface $\alpha$. The common-tangent construction from $f_{\alpha\beta}$ and $f_\gamma$ at 900, 950, 1000 K compared with the assessed U–Zr diagram is the decisive thermodynamic test and costs nothing. |
| R12 | §3.1, §9 | The factor `scale = 0.035` was recorded but no physical hypotheses were listed for a 30-fold discrepancy. | Candidate causes are ranked so that the student tests them in order instead of treating `scale` as a free parameter. |
| R13 | §9 | Mohanty's Eq. (8) imposes both Dirichlet temperatures and zero heat flux on a 1-D interval, which is over-determined; not listed. | Only the two Dirichlet values can be used in 1-D. |

## 1. Planning Principles

### 1.1 Reproduction Objectives

The main goal of this project is not to build a “maximal model” that immediately includes every possible physical process. It is to build a U–Zr constituent-redistribution framework that is traceable, verifiable, and expandable layer by layer. Each layer adds only one physical process, or one tightly related group of processes, and the following must be completed before moving to the next layer:

1. Equation and symbol checks;
2. Unit and dimension checks;
3. Conservation checks;
4. Mesh and time-step convergence checks;
5. Quantitative comparison with paper curves or experimental data;
6. Classified records of parameter sources, fitted parameters, and parameters that cannot be determined.

The final framework should be able to answer three levels of questions:

- Mohanty 2011: in single-phase $\gamma$-U–Zr, how do the chemical-potential gradient and the temperature gradient jointly determine Zr migration?
- Wen 2022: how do multiphase thermodynamics, KKS phase equilibrium, and irradiation-enhanced diffusion change the macroscopic radial constituent distribution?
- Jung 2025: how does the feedback among composition, phase fractions, porosity, sodium infiltration, and thermal conductivity change the temperature field and constituent redistribution?

### 1.2 Coupling Boundaries That Do Not Go Beyond the Papers

This plan does not introduce the following physics outside the paper models:

- Do not explicitly solve rate-theory equations for vacancies and interstitials;
- Do not use the phase field to explicitly resolve pore nucleation, growth, coalescence, or migration;
- Do not let porosity directly change the diffusivity, free energy, or mechanical properties;
- Do not add stress, cracking, fuel–cladding interaction, or fission-gas transport;
- Do not interpret the two-dimensional extension as a new physical model; the two-dimensional stage only checks geometric and numerical effects.

In Wen 2022, irradiation acts only by modifying the chemical mobility through the irradiation-enhancement factor $\xi$. In Jung 2025, porosity is computed only from phase fractions and affects the temperature field only through the effective thermal conductivity.

### 1.3 Unit Strategy

Different papers and reproduction directories may use different unit systems. A single unit system is not required for the whole project:

- `s1_mohanty` may continue to use $\mathrm{\mu m}$–s–eV/atom–K;
- Wen 2022 and Jung 2025 should follow the papers and use m–s–J/mol or J/m³–K;
- Do not mix eV/atom, J/atom, J/mol, and J/m³ in the same input file without an explicit conversion;
- Each reproduction directory must contain a standalone `UNITS.md` or a units section in the README;
- Each Material property name should carry its physical meaning, and units should be noted in comments when needed.

The unit audit for each directory must at least answer:

1. Is the free energy per atom, per mole, or per volume?
2. What are the units of the chemical-potential variable?
3. Does the product of mobility and chemical-potential gradient give the correct flux?
4. Does the product of thermal mobility and $\nabla T/T$ have the same dimensions as the chemical flux?
5. Do coordinates and gradients use m or $\mathrm{\mu m}$?
6. Do thermal conductivity, volumetric heat source, density, and specific heat form a consistent heat equation?
7. Are the atomic mobilities and heats of transport on a per-atom or a per-mole basis, and is that basis the same as the free energy's?

> Revision R1. The seventh question is added because the three papers do not share a basis. Mohanty's Fig. 2 mobilities $\beta_i\sim10^7$–$10^8\ \mathrm{m^2\,J^{-1}\,s^{-1}}$ are per atom, $\beta_i=D_i/(k_BT)$, and the heats of transport $\widetilde Q_i^*$ are in J/atom. Wen Eq. (27) and Jung Eq. (17) write $M_c^\gamma=V_m c_\gamma(1-c_\gamma)[\ldots\beta\ldots]$ with $V_m$ in $\mathrm{m^3\,mol^{-1}}$ and free energies in J/mol; there $\beta_i$ must be on a molar basis, $\beta_i^{\mathrm{mol}}=D_i/(RT)=\beta_i^{\mathrm{atom}}/N_A$, and $\widetilde Q_i^{*,\mathrm{mol}}=N_A\widetilde Q_i^{*,\mathrm{atom}}$. Check of magnitude: Mohanty's $\widetilde Q_{\mathrm{Zr}}^*=-17.5\times10^{-20}\ \mathrm J$ per atom is $-105\ \mathrm{kJ\,mol^{-1}}$, which is the magnitude Wen tabulates as $1.05\times10^5\ \mathrm{J\,mol^{-1}}$; Mohanty's $\widetilde Q_U^*=+2.5\times10^{-20}\ \mathrm J$ is $+15\ \mathrm{kJ\,mol^{-1}}$, which does not match Wen's $1.1\times10^5\ \mathrm{J\,mol^{-1}}$ (Hofman 1996 value). A basis error here changes every mobility by $N_A$ and will masquerade as a "time-scale problem". Each `UNITS.md` must contain a table with, for every kinetic parameter, the paper value, the paper basis, the code value, and the conversion factor.

## 2. Relationship Among the Models in the Three Papers

### 2.1 Mohanty et al. (2011): Single-Phase Thermotransport Benchmark

The paper considers only single-phase bcc-$\gamma$ U–Zr. The conserved variable is the Zr atom fraction $c$, and the governing equation is

$$
\frac{1}{V_m}\frac{\partial c}{\partial t}
=
\nabla\cdot\left[
V_m M_c\nabla\left(
\frac{\partial f}{\partial c}-2\kappa_c\nabla^2c
\right)
-M_Q\frac{\nabla T}{T}
\right].
$$

The chemical mobility and thermal mobility are

$$
M_c=\frac{1}{V_m}c(1-c)
\left[c\beta_U+(1-c)\beta_{\mathrm{Zr}}\right],
$$

$$
M_Q=\frac{1}{V_m}c(1-c)
\left[\beta_U\widetilde Q_U^*
-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right].
$$

Atomic mobilities follow

$$
\beta_i=\beta_{0,i}\exp\left(-\frac{H_i}{RT}\right).
$$

> Revision R2. Equations (2) and (3) as printed are correct and must be implemented verbatim. They follow from the lattice-frame fluxes $\mathbf J_i=-\beta_ic_i\left[\nabla\mu_i+\widetilde Q_i^*\nabla T/T\right]$ and the Darken interdiffusion flux $\widetilde{\mathbf J}_{\mathrm{Zr}}=(1-c)\mathbf J_{\mathrm{Zr}}-c\,\mathbf J_U$. For the thermal part this gives directly
>
> $$
> \widetilde{\mathbf J}_{\mathrm{Zr}}^{T}
> =-\frac{c(1-c)}{T}
> \left[\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*-\beta_U\widetilde Q_U^*\right]\nabla T
> =+\frac{c(1-c)}{T}
> \left[\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right]\nabla T,
> $$
>
> which is Eq. (3) with no further composition weights. The weights $c\beta_U+(1-c)\beta_{\mathrm{Zr}}$ in Eq. (2) arise only in the chemical part, from eliminating $\nabla\mu_U$ with the Gibbs–Duhem relation and converting $\nabla\mu_{\mathrm{Zr}}$ to $\nabla(\partial f/\partial c)$. The Darken-type weighting that the student inserted into $M_Q$ in `s1_mohanty.i` therefore has no thermodynamic basis and is not a "correction" of the paper. The README's stated reason, that $\rho=M_Q/M_c$ "diverges" as $c\to0$, is false:
>
> $$
> \lim_{c\to0}\frac{M_Q}{M_c}
> =\frac{\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*}{\beta_{\mathrm{Zr}}}
> =\frac{\beta_U}{\beta_{\mathrm{Zr}}}\widetilde Q_U^*+\lvert\widetilde Q_{\mathrm{Zr}}^*\rvert,
> $$
>
> which is finite ($\approx2.3\ \mathrm{eV}$ at 1400 K with the student's numbers). The deep cold-end depletion the student observed is a genuine prediction of Eq. (3) with those parameters, not a singularity, and must be addressed by revisiting the parameters or the thermodynamics, not by re-weighting the Onsager coefficient. Only `paper_equation.i` is admissible as a Mohanty reproduction; `calibrated.i` is a phenomenological fit and must be labeled as such everywhere.

> Revision R3. Eq. (1) drives the chemical flux with the *full* gradient $\nabla(\partial f/\partial c)$, evaluated in a temperature field. That gradient contains
>
> $$
> \nabla\frac{\partial f}{\partial c}
> =\frac{\partial^2f}{\partial c^2}\nabla c
> +\frac{\partial^2f}{\partial c\,\partial T}\nabla T,
> $$
>
> and the second term is itself a $\nabla T$-driven flux. For an ideal solution it equals $k_B\ln\frac{c}{1-c}\nabla T$, about 5 % of the $\widetilde Q^*\nabla T/T$ term at $c=0.39$; the pure-element end-member difference $G^{0}_{\mathrm{Zr}}(T)-G^{0}_{U}(T)$ contributes a comparable amount through its entropy. In the irreversible-thermodynamics definition of the heat of transport the driving force is the *isothermal* gradient $\nabla_T\mu$, and literature $\widetilde Q^*$ values (Campbell–Huntington, D'Amico–Huntington, Hofman's fit) are defined against that convention. Consequences for this plan:
>
> - The strict Mohanty branch uses the full $\nabla(\partial f/\partial c)$, because that is what Eq. (1) literally says and what FiPy would have computed;
> - The end-member terms $G_U^{0}(T)$ and $G_{\mathrm{Zr}}^{0}(T)$ must be kept in that branch. The current `s1_mohanty.i` drops them and the README says they are "absorbed into the effective $\widetilde Q^*$"; that is a model change and must be moved to the calibrated branch;
> - A documented sensitivity branch, `isothermal_gradient`, may replace $\nabla(\partial f/\partial c)$ by $\frac{\partial^2f}{\partial c^2}\nabla c$ to quantify the ambiguity of the $\widetilde Q^*$ convention;
> - The same question recurs in Wen and Jung, where the KKS chemical potential $\partial f_\gamma/\partial c_\gamma$ carries temperature through the end-members (and, in Jung, only through the end-members; see R10).

The paper explicitly states that the compositional gradient-energy coefficient vanishes in a single-phase alloy. The strict paper branch should therefore use

$$
\kappa_c=0,
\qquad
\mu=\frac{\partial f}{\partial c}.
$$

The current `s1_mohanty.i` uses a nonzero `kappa_c = 1e-5` as numerical regularization. That setting may be kept only in the calibration/regularization branch, and a $\kappa_c\to0$ convergence check should confirm that it does not control the boundary layer or the interior minimum in Fig. 6.

Key benchmarks:

- One-dimensional length of 300 $\mathrm{\mu m}$;
- Initial composition of 39 at.% Zr;
- End temperatures of about 1050 K and 1405 K;
- Concentration profiles at 5, 10, and 30 days;
- After 30 days, about 11 at.% Zr enrichment at the hot end and about 7 at.% Zr depletion at the cold end;
- Diffusion-couple initial compositions of 39 at.% and 22.5 at.% Zr;
- In the paper’s diffusion-couple case, the thermotransport flux is about 4–5 orders of magnitude larger than the chemical-potential-gradient flux.

The paper temperature field is obtained from $\nabla^2T=0$, and thermal conductivity is assumed independent of temperature and composition. Composition-dependent thermal conductivity should therefore not be introduced at this stage.

### 2.2 Wen et al. (2022): Multiphase KKS and Irradiation-Enhanced Diffusion

Wen 2022 adds the following on top of a Mohanty-type conservation equation:

- Macroscopic representations of the $\alpha$, $\beta$, and $\gamma$ phases;
- One nonconserved $\gamma$-phase order parameter;
- KKS intra-phase concentration constraints;
- Temperature-dependent phase-interface mobility;
- Irradiation-enhanced chemical mobility;
- Analysis of the macroscopic interface width and interface parameters.

The constituent evolution equation is

$$
\frac{\partial c}{\partial t}
=
\nabla\cdot\left(
M_c^I\nabla\frac{\delta F}{\delta c}
-M_T\frac{\nabla T}{T}
\right),
$$

where

$$
M_c^I=\xi M_c,
\qquad
\xi=\frac{C_v^e+C_v^r}{C_v^e},
$$

$$
C_v^e=A_v\exp\left(-\frac{E_v^f}{k_BT}\right).
$$

> Revision R4. The $\gamma$-phase mobilities in Wen Eqs. (27) and (31), reused verbatim in Jung Eqs. (17)–(18), are
>
> $$
> M_c^\gamma=V_m\,c_\gamma(1-c_\gamma)\left[c_\gamma\beta_U+(1-c_\gamma)\beta_{\mathrm{Zr}}\right],
> \qquad
> M_T^\gamma=V_m\,c_\gamma(1-c_\gamma)\left[\beta_U\widetilde Q_U^*+\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right],
> $$
>
> with a *plus* sign between the two terms and *positive* tabulated magnitudes ($\widetilde Q_U^*=1.1\times10^5$, $\widetilde Q_{\mathrm{Zr}}^*=1.05\times10^5\ \mathrm{J\,mol^{-1}}$). Mohanty Eq. (3) has a *minus* sign and a *negative* $\widetilde Q_{\mathrm{Zr}}^*=-17.5\times10^{-20}\ \mathrm J$. The two conventions coincide for the Zr term only if one interprets Wen's positive $\widetilde Q_{\mathrm{Zr}}^*$ as the magnitude of a heat of transport that is negative in Mohanty's convention. They do *not* coincide for the U term, whose sign is the same in both papers but whose magnitude differs by a factor of seven (R1). The plan must never carry a signed Mohanty $\widetilde Q_i^*$ into a Wen/Jung input file or vice versa: doing so reverses the Zr thermal flux in the $\gamma$ phase. Each stage's `UNITS.md` states the sign convention of its $M_T$ explicitly, and the Stage 2 verification test (the constant-$\nabla T$ Soret steady state) is run once with each convention so the sign of the resulting Zr enrichment is checked against the paper's own figures: Zr enriches at the *hot* center in all three papers.

The paper uses $E_v^f=1.20$ eV and cites $C_v^r=7.04\times10^{-7}$ at about 1000 K and 30 kW/m. The paper does not explicitly solve defect-concentration evolution in the phase-field calculation, so the reproduction should treat $C_v^r$ as a prescribed value or a parameter function, not as a new defect PDE.

> Revision R5. With $E_v^f=1.20\ \mathrm{eV}$, $C_v^e=\exp(-E_v^f/k_BT)$ is $1.9\times10^{-7}$ at 900 K and $7.6\times10^{-7}$ at 988 K, so $\xi\approx4.7$ at the surface and $\approx1.9$ at the center. Two physical consequences that Revision 1 did not draw:
>
> 1. In Wen's Eqs. (24) and (29)–(31) the factor $\xi$ multiplies $M_c$ only. In a vacancy-mediated mechanism the same excess vacancies enhance every $\beta_i$, and $M_T\propto\beta_i\widetilde Q_i^*$ would scale identically. Enhancing $M_c$ alone changes the steady-state balance $M_c\nabla\mu=M_T\nabla T/T$ to $\nabla\mu=(M_T/\xi M_c)\nabla T/T$, i.e. it reduces the final Zr redistribution by a factor $\xi$ while shortening the time to reach it. Wen's Fig. 2 versus Fig. 3 shows *more* redistribution with $\xi$ only because the simulated time, $1.3\times10^5\ \mathrm s\approx1.5\ \mathrm{day}$, is far from steady state. The literal branch must follow Wen; a `xi_both` branch with $M_T^I=\xi M_T$ is mandatory as a physics-motivated sensitivity, and every plot must state which branch it uses.
> 2. The DP-81 pin was irradiated for on the order of $10^2$ effective full-power days. Matching its PIE profile after 1.5 simulated days means the model still contains an implicit acceleration of order $10^2$, whatever its origin. Wen's claim that RED "replaces the artificial increase of the diffusion coefficient" is therefore not established by the paper, and this plan must not repeat it as a conclusion. The $C_v^r$ value was also computed for 30 kW/m, whereas Jung quotes 24 kW/m for DP-81.

The phase order parameter satisfies the Allen–Cahn equation:

$$
\frac{\partial\gamma}{\partial t}
=-L(T,c_\gamma)\frac{\delta F}{\delta\gamma},
$$

$$
L=L_0\exp\left(-\frac{Q_\gamma}{RT}\right),
$$

$$
Q_\gamma=128000-107000c_\gamma+174000c_\gamma^2.
$$

The total free energy contains the phase bulk free energies, a double-well potential, the compositional gradient energy, and the phase-field gradient energy:

$$
F=\int_V
\left\{
\frac{(1-h)f_{\alpha\beta}+hf_\gamma}{V_m}
+\omega g(\gamma)
+\frac{\kappa_c}{2}\lvert\nabla c\rvert^2
+\frac{\kappa_\gamma}{2}\lvert\nabla\gamma\rvert^2
\right\}\,\mathrm dV.
$$

The KKS constraints are

$$
c=(1-h)c_{\alpha\beta}+hc_\gamma,
$$

$$
\frac{\partial f_{\alpha\beta}}{\partial c_{\alpha\beta}}
=
\frac{\partial f_\gamma}{\partial c_\gamma}.
$$

The $\alpha/\beta$ weights in Wen equations (7)–(8) contradict the immediately adjacent text and the U–Zr phase diagram: according to the printed formulas, the high-temperature region would be labeled $\alpha$ and the low-temperature region would be labeled $\beta$. Both of the following must therefore be retained:

- `literal`: use the printed formulas of the paper verbatim;
- `phase_consistent`: swap the $\alpha$ and $\beta$ weights so that $T<935\,\mathrm K$ is $\alpha$ and $T>935\,\mathrm K$ is $\beta$.

> Revision R6. The decision is not open. Wen's Fig. 2–4 and Jung's Fig. 5–7 all show the $\alpha+\gamma$ two-phase zone at the cold surface and $\beta+\gamma$ inside, in agreement with the U–Zr phase diagram (the $\alpha\to\beta$ transition of the U-rich phase lies at 935 K). `phase_consistent` is the primary branch of every KKS stage. `literal` is retained only as a one-off check whose purpose is to document that the printed formulas invert the zone order, and it must not be used for any curve compared with experiment. Two further points:
>
> - The switching function $h_{\alpha\beta}=\tfrac12\left[1+\tanh\frac{T-935}{2}\right]$ has an implicit half-width of 2 K. The mesh must resolve this: with a slug temperature drop of about 90 K over 1.085 mm the isotherm band is roughly 50 $\mathrm{\mu m}$ wide, so the radial mesh must be at most 10–20 $\mathrm{\mu m}$ near the 935 K isotherm or the switch must be smoothed. A sensitivity on the width (2, 5, 10 K) belongs in Stage 3.
> - The switch acts on the local temperature, i.e. on position via $T(r)$. Where $T(r)$ evolves (Jung), the $\alpha/\beta$ boundary moves; the plan should check that the KKS phase concentration $c_{\alpha\beta}$ does not jump when the switch passes through a node.

The phase-region order in the paper’s main figures decides which of the two is used. They must not be swapped in the code without a record.

Wen 2022 uses a prescribed radial temperature profile and does not solve fully coupled heat conduction with composition-dependent thermal conductivity. The irradiation-enhancement stage should therefore keep the temperature field prescribed and avoid introducing the thermal feedback of Jung 2025 at the same time.

> Revision R7. Wen describes the domain as a one-dimensional *axisymmetric* fuel slug of radius 2.17 mm with a prescribed radial $T(r)$. The divergence operator in the Cahn–Hilliard equation is therefore the cylindrical one, $\frac1r\frac{\partial}{\partial r}\left(r\,J_r\right)$, and mass conservation must be checked with the weight $2\pi r\,\mathrm dr$, not $\mathrm dr$. In MOOSE this is `Problem/coord_type = RZ` with the 1-D mesh along $r$ (or a thin $r$–$z$ strip), and all `ElementIntegralVariablePostprocessor` values automatically acquire the $2\pi r$ weight. A Cartesian slab mesh gives a different $c(r)$ for the same fluxes because the outward-facing area grows with $r$. The slab form is *not* an acceptable primary geometry for Wen or Jung; see Stage 5 for the quantitative consequence on the temperature field.

### 2.3 Jung et al. (2025): Heat Conduction, Porosity, and Sodium-Infiltration Feedback

Jung 2025 continues the KKS multiphase phase-field model and solves heat conduction:

$$
\rho c_p\frac{\partial T}{\partial t}
=
\nabla\cdot\left[
k(P,P_{\mathrm{Na}},T,c)\nabla T
\right]
+\dot q(B).
$$

> Revision R8. Jung writes the transient form, but the thermal time scale makes it irrelevant. With $k\approx30\ \mathrm{W\,m^{-1}\,K^{-1}}$, $\rho\approx16\times10^3\ \mathrm{kg\,m^{-3}}$, $c_p\approx170\ \mathrm{J\,kg^{-1}\,K^{-1}}$ the thermal diffusivity is $\approx1.1\times10^{-5}\ \mathrm{m^2\,s^{-1}}$ and the slug relaxation time $r_0^2/\alpha_T\approx0.4\ \mathrm s$, against redistribution times of $10^5$–$10^7\ \mathrm s$. The temperature field is therefore in quasi-steady equilibrium with the instantaneous $k(c,\phi,T)$ and $\dot q(B)$ at every step. The primary implementation is the *steady* heat equation
>
> $$
> 0=\nabla\cdot\left[k\,\nabla T\right]+\dot q(B),
> $$
>
> solved in the same nonlinear system as $c,w,\phi$ (it costs one extra variable, needs no $\rho c_p$, and eliminates the "unsourced parameter" problem the student flagged in §7.2). If a transient form is wanted for comparison, $\rho=15.8$–$16.0\ \mathrm{g\,cm^{-3}}$ and $c_p=150$–$200\ \mathrm{J\,kg^{-1}\,K^{-1}}$ for U–10Zr at 900–1000 K may be used, tagged `cited`, and the result must be indistinguishable from the steady branch after the first step.

> Revision R7 (continued). Jung defines the heat source through $\pi r_0^2$, i.e. as a *cylindrical* average, and the temperature drop reported (990 K center, 900 K surface for DP-81) is only reproducible in cylindrical coordinates. For uniform $\dot q$ and constant $k$, the centerline rise is $\dot q r_0^2/(4k)$ in a cylinder and $\dot q r_0^2/(2k)$ in a slab of half-width $r_0$ with the same surface temperature. With $\dot q\approx1.6\times10^9\ \mathrm{W\,m^{-3}}$, $r_0=2.17\ \mathrm{mm}$, $k\approx30\ \mathrm{W\,m^{-1}\,K^{-1}}$ this is 63 K versus 126 K. Only the cylindrical value is consistent with the 90 K the paper reports once the composition dependence of $k$ is included. Stage 5 and all later stages therefore run in `coord_type = RZ`; the slab is not a permitted primary geometry.

The effective thermal conductivity is

$$
k(P,P_{\mathrm{Na}},T,c)
=P_c(P,P_{\mathrm{Na}},T,c)(1-P)^{1.5}k_0(T,c).
$$

The matrix thermal conductivity is written as

$$
k_0(T,c)=a_1(c)+a_2(c)T+a_3T^2.
$$

For U–Zr, $W_{\mathrm{Pu}}=0$, and paper equations (28)–(30) should be arranged as

$$
a_1(c)=17.5
\frac{1-2.23W_{\mathrm{Zr}}(c)}
{1+1.61W_{\mathrm{Zr}}(c)},
$$

$$
a_2(c)=1.54\times10^{-2}
\frac{1+0.061W_{\mathrm{Zr}}(c)}
{1+1.61W_{\mathrm{Zr}}(c)},
$$

$$
a_3=9.38\times10^{-6}.
$$

The paper text first calls the coefficients $a_0,a_1,a_2$, while the subsequent formulas are numbered $a_1,a_2,a_3$. The code should use one consistent naming and keep a mapping to the paper numbering in the README.

Porosity in the paper is not an independently evolving variable. It is an algebraic function of the phase fractions:

$$
P=0.4\eta_\gamma+0.13(\eta_\alpha+\eta_\beta).
$$

The sodium thermal conductivity is

$$
k_{\mathrm{Na}}(T)
=93-0.0581(T-273.15)
+1.173\times10^{-5}(T-273.15)^2.
$$

The porosity correction factor is implemented according to paper equation (31). 

> Revision R9. Equation (31) is typographically ambiguous in the PDF: it is not clear whether the factor $(1-P)$ that multiplies the bracket sits in the numerator or the denominator of the prefactor. The two readings differ by $(1-P)^2\approx0.4$–$0.75$ and one of them is unphysical. Before coding, both readings are written out and subjected to the following checks, which any correct sodium correction must satisfy:
>
> 1. $P_c\to1$ when the sodium-filled fraction $P_{\mathrm{Na}}\to0$ (no sodium, no correction);
> 2. $P_c\ge1$ whenever $k_{\mathrm{Na}}>k_0$, because filling an insulating pore with a better conductor cannot lower the conductivity. At 900 K, $k_{\mathrm{Na}}\approx61\ \mathrm{W\,m^{-1}\,K^{-1}}$ and $k_0(W_{\mathrm{Zr}}=0.10)\approx31\ \mathrm{W\,m^{-1}\,K^{-1}}$, so $P_c>1$ is required over the whole slug;
> 3. $P_c(P_{\mathrm{Na}}=P)$, all pores filled, should reproduce a standard two-phase (e.g. Maxwell–Eucken) bound to within a few percent.
>
> The reading that passes all three is the primary one; the other is recorded in the README as rejected with the numbers. This test also fixes the definition of $P_{\mathrm{Na}}$: whether it is the sodium fraction of the pore volume (bounded by 1) or of the total volume (bounded by $P$). Jung's burnup ramp (1 % $\to$ 1.5 % B) and the two-branch treatment in Stage 7 (`na_literal`, `na_pore`) depend on this choice.

The treatment of porosity with burnup is:

- 0–1% burnup: porosity increases linearly from 0 to the target value;
- Above 1% burnup: porosity is held at the target value;
- After porosity has stabilized, sodium infiltration begins to increase;
- Sodium infiltration reaches its target value at 1.5% burnup and then remains constant.

The volumetric heat source is approximated by

$$
\dot q(B)=\frac{\mathrm{LHGR}(1-B)}{\pi r_0^2}
$$

Burnup increases linearly over the total simulation time to 3.6 at.% for DP-81 or 7.7 at.% for DP-11.

The feedback loop in this model is:

$$
c,\phi,T
\longrightarrow
\eta_\alpha,\eta_\beta,\eta_\gamma
\longrightarrow
P
\longrightarrow
k
\longrightarrow
T
\longrightarrow
f,M_c,M_T,L
\longrightarrow
c,\phi.
$$

Note: the paper explicitly states that porosity affects only thermal conductivity and does not directly affect diffusion kinetics or phase equilibrium. This reproduction should keep that restriction.

> Revision R6 (Jung). Jung's Eqs. (3)–(4) reproduce Wen's $\alpha/\beta$ switching formulas verbatim, including the inversion; Jung's Fig. 5–7 again show $\alpha$ at the cold surface. The same `phase_consistent` primary branch applies. Jung's double-well is printed as $g(\phi)=\phi^2(1-\phi^2)$, which is not a double well on $[0,1]$ (it has a single interior maximum and is negative for $\phi>1$); the standard $g=\phi^2(1-\phi)^2$ that Wen uses is almost certainly intended and is the primary form. Record the printed form as a suspected typo.

> Revision R10. Two model choices peculiar to Jung must be made explicit because they change the physics relative to Wen:
>
> 1. Jung evaluates the ideal-mixing and Redlich–Kister parts of $f_\gamma$ at a *fixed* $T_{\mathrm{ref}}=1010\ \mathrm K$, so that $\partial^2 f_\gamma/\partial c\,\partial T$ receives contributions only from the SGTE end-members $G_U^0(T)$ and $G_{\mathrm{Zr}}^0(T)$. This removes the configurational-entropy part of the chemical $\nabla T$ driving force discussed in R3 while keeping the end-member part. Since $T_{\mathrm{ref}}$ lies 20 K above the hottest point in the slug, this is a small numerical convenience for the $\gamma$ phase, but it means the Wen and Jung $\gamma$ free energies are *different functions* and cannot be interchanged between Stages 3–4 and 5–8 without a note. The `f_gamma` material must carry a `T_mode = {local, fixed_1010}` switch.
> 2. Jung uses a constant $L=2.5\times10^{-10}\ \mathrm{m^3\,J^{-1}\,s^{-1}}$, an interface width $l=0.5\ \mathrm{mm}$ (23 % of $r_0$) and $\sigma=0.5\ \mathrm{J\,m^{-2}}$, whereas Wen uses a thermally activated $L$ with $\kappa_\gamma=7.8\times10^{-4}\ \mathrm{J\,m^{-1}}$ and $\omega=5\times10^5\ \mathrm{J\,m^{-3}}$, i.e. an interface of about 40 $\mathrm{\mu m}$. Jung's interface is an order of magnitude wider; it smears the two-phase zone boundaries and relaxes the mesh requirement from about 8 $\mathrm{\mu m}$ (Wen) to about 100 $\mathrm{\mu m}$. Both sets are kept, under `interface = {wen, jung}`, and the zone-boundary positions are compared at equal time. Jung's constant $L$ also gives an independent order-of-magnitude bound for Wen's unreported $L_0$: $L_0\approx L\exp(Q_\gamma/RT)$ at 950 K and $c\approx0.3$.

## 3. Positioning of the Current First-Step Implementation and Items That Need Correction

The current directory `palm/problem/uzr/s1_mohanty` already implements:

- A one-dimensional 300 $\mathrm{\mu m}$ mesh;
- Split CH for the Zr atom fraction $c$ and the chemical potential $w$;
- A linearly prescribed temperature field;
- A CALPHAD-type $\gamma$-phase free energy;
- Temperature-dependent atomic mobility;
- MOOSE’s built-in `SoretDiffusion`;
- Output at 5, 10, and 30 days;
- Postprocessing of the mass average, endpoint concentrations, and spatial profiles.

This work can serve as a “phenomenon-reproduction branch,” but a “paper-equation branch” should be established before later stages, for the reasons below.

### 3.1 The Current Model Contains Three Artificial Corrections

The current input file uses:

- An overall atomic-mobility scale factor `scale = 0.035`;
- An added composition weight on Mohanty equation (3);
- Omission of the pure-element end-member terms $G_U^0(T)$ and $G_{\mathrm{Zr}}^0(T)$ from the free energy, justified in the README as "absorbed into the effective $\widetilde Q^*$".

The revision notes in the current README still retain `scale = 0.05`, while the parameter table and the actual input file use `0.035`. Later cleanup should take the input file as authoritative and correct the old value in the documentation. The current `kappa_c` comment is written as “thermal conductivity”; it should be the compositional gradient-energy coefficient. The comment is also inconsistent with the paper's statement that $\kappa_c=0$ in the single-phase alloy; the nonzero value $10^{-5}$ is a numerical regularization and must be labeled as such.

> Revision R2/R3. The second and third items are not "stabilizations"; they change the model. The composition weight has no thermodynamic basis (see §2.1, R2), and the end-member omission removes part of the $\nabla T$-driven chemical flux (R3). Both are removed from `paper_equation.i`. The student's argument that the un-weighted ratio $M_Q/M_c$ diverges at $c\to0$ is arithmetically wrong; the limit is $\frac{\beta_U}{\beta_{\mathrm{Zr}}}\widetilde Q_U^*+\lvert\widetilde Q_{\mathrm{Zr}}^*\rvert$. If `paper_equation.i` drives $c$ to values that trip `SplitCHParsed` through $\ln c$, the correct remedy is a bounded logarithm (`log_tol`) or a smaller `dtmax`, not a change of the Onsager coefficient.

> Revision R12. A factor of 0.035, i.e. a 30-fold slow-down of the paper's own mobilities, requires an explanation before it is accepted even in the calibrated branch. Ranked hypotheses to test in Stage 1, in this order:
>
> 1. *Basis error* (R1): if the Fig. 2 $\beta_i$ are per atom and the free energy is per atom, no $N_A$ appears; but if `JeV = 96485` converts a per-mole Redlich–Kister energy to eV per atom while $\beta$ was read as $\mathrm{m^2\,J^{-1}\,s^{-1}}$ per atom, the product $M_c\,\partial^2f/\partial c^2$ is consistent, while any $V_m$ factor would not be. Write out the units of $M_c\nabla\mu$ in $\mathrm{m\,s^{-1}}$ explicitly.
> 2. *Digitization of Fig. 2*: the Arrhenius fit reads $\beta_{0}$ off a log plot; a 0.3-decade reading error gives a factor 2, not 30, so this cannot explain the full factor but must be quantified.
> 3. *Unstated nondimensionalization*: Mohanty's FiPy model is written in dimensionless variables with a length scale $l$, energy scale $\Delta f$ and mobility scale that are not all reported. If the paper's simulated 5/10/30 days are dimensionless times reconverted with an inconsistent scale, no MOOSE input will match without a factor. Estimate the paper's own time scale $L^2/D_U$: with $D_U(1400\ \mathrm K)\approx3.7\times10^{-12}\ \mathrm{m^2\,s^{-1}}$ and $L=300\ \mathrm{\mu m}$ this is 0.28 day, so the paper's profiles should be near steady state after 5 days and essentially frozen between 10 and 30 days; the paper's Figs. 6–7 show continued evolution to 30 days. That is the root of the disagreement, and it argues that the paper's kinetics are internally slower than its stated $\beta_i$ by a factor of order 10–30, consistent with the student's fitted `scale`.
> 4. *Thermodynamic factor*: the $\gamma$-phase $\partial^2f/\partial c^2$ at 39 at.% Zr and 1050–1400 K is reduced by the nearby $\gamma_1/\gamma_2$ miscibility gap; check whether the student's Redlich–Kister parameters (Quaini et al.) give a factor smaller than Mohanty's source (Sheldon & Peterson / Leibowitz). This affects the chemical back-flux but not $M_T$, so it changes the steady-state profile rather than the time scale.
>
> Until hypothesis 3 is confirmed or excluded, `scale` is reported as an empirical time-scale factor of the *paper's* simulation, not as a correction to the alloy's mobility.

These treatments improve the fit to Figs. 6/7, but they are not the original Mohanty 2011 equations. Two parallel cases should be kept:

- `paper_equation.i`: strictly use paper equations (2)–(4), with no composition weight, with the end-member terms retained, with $\kappa_c=0$, and no prior fitting of the time scale;
- `calibrated.i`: keep the current stabilized and fitted version, with all three modifications listed in its header.

Every figure must clearly mark whether the “original equation” or the “calibrated equation” was used, so that a fitted result is not mislabeled as a parameter-free reproduction.

### 3.2 The Current `SoretDiffusion` Is Unsuitable for Later Strong Coupling

The residual of MOOSE’s built-in `SoretDiffusion` is

$$
R_T=
\int_\Omega
\frac{DQc}{k_BT^2}\nabla T\cdot\nabla\psi\,\mathrm dV.
$$

It is currently mapped to the target thermal-migration flux through

$$
D_{\mathrm{Soret}}=\frac{M_ck_BT}{c},
\qquad
Q_{\mathrm{eff}}=-\frac{M_T}{M_c}
$$

The mapping is algebraically feasible, but the built-in Kernel Jacobian explicitly treats only $c$ and $T$ and does not include the material derivatives of $D(c,T)$ and $Q(c,T)$. Using PJFNK can avoid some of the problems caused by an incomplete analytic Jacobian, but in the strongly coupled $c$–$\phi$–$T$ problem of Jung 2025 it leads to:

- Slower Newton convergence;
- A preconditioner matrix that lacks key coupling blocks;
- Difficulty distinguishing physical stiffness from an incomplete Jacobian;
- Harder maintenance once material properties depend on intra-phase concentrations.

A new AD Kernel that uses $M_T$ directly should be added in PALM, for example `ADThermalMigration`, whose weak-form residual is

$$
R_{\mathrm{TM}}
=-\int_\Omega
\frac{M_T}{T}\nabla T\cdot\nabla\psi\,\mathrm dV.
$$

This object should:

- Inherit `ADKernel`;
- Couple the temperature variable `T`;
- Read the AD material property `thermal_mobility`;
- Let automatic differentiation generate the full Jacobian with respect to $c$, $\phi$, $T$, $c_{\alpha\beta}$, and $c_\gamma$;
- Not hard-code the Boltzmann constant or a unit system inside the Kernel;
- Not artificially split $M_T$ into `D_soret` and `Qeff`.

`SoretDiffusion` can continue to be used for the first-step comparison, but the upstream MOOSE source should not be modified. The new object should be placed in PALM:

- `include/kernels/ADThermalMigration.h`
- `src/kernels/ADThermalMigration.C`
- `test/tests/kernels/ad_thermal_migration/`

PALM already registers this application’s objects through `Registry::registerObjectsTo(f, {"palmApp"})`, so using

```cpp
registerMooseObject("palmApp", ADThermalMigration);
```

in the new class source file is enough to include it in the application.

### 3.3 Completion Criteria for the First Step

Before moving on to Wen 2022, the following should be completed:

- A side-by-side comparison of the original-equation and calibrated-equation cases;
- With thermal migration turned off under isothermal conditions, confirm that a uniform composition remains unchanged;
- Set $Q_U^*=Q_{\mathrm{Zr}}^*=0$ and confirm that the thermal-migration term is strictly zero;
- Reverse the temperature gradient and confirm that the migration direction reverses;
- Check the relative drift of $\int_\Omega c\,\mathrm dV$;
- Use at least 150, 300, and 600 elements for mesh convergence;
- Use at least three sets of `dtmax` for time-step convergence;
- Quantify the hot-end, cold-end, and profile errors separately;
- Add output of the chemical flux, the thermal-migration flux, and the total flux.

## 4. Recommended Directory Structure

Independent, runnable stage directories should be created under `palm/problem/uzr`:

```text
uzr/
├── s1_mohanty/
│   ├── paper_equation.i
│   ├── calibrated.i
│   ├── README.md
│   ├── UNITS.md
│   └── scripts/
├── s2_mohanty_diffusion_couple/
├── s3_kks_multiphase/
├── s4_wen_irradiation/
├── s5_jung_heat_baseline/
├── s6_jung_porosity/
├── s7_jung_sodium/
├── s8_dp11_validation/
└── s9_2d_geometry/
```

Each directory should contain at least:

- The main input file;
- A minimal verification input file;
- `README.md`;
- `UNITS.md`;
- A parameter-source file;
- Scripts for automatic extraction and plotting;
- Mesh/time-step convergence scripts;
- A MOOSE tests entry point, or test notes that point to PALM `test/tests`.

## 5. Staged Reproduction Roadmap

## Stage 0: Build a Repeatable Verification Infrastructure

### Objective

Before adding more physics, automate execution, sampling, digitization of paper curves, and error calculation.

### Work

1. Digitize the target curves of each paper into CSV and record the figure number, axis units, and manual-reading error.
2. Report the following metrics uniformly:
   - Total Zr amount;
   - Minimum and maximum $c$;
   - Center and surface $c$;
   - Center and surface temperature;
   - Volume fraction of each phase;
   - Phase-region boundary locations;
   - $L_1$, $L_2$, and maximum relative error.
3. Save the MOOSE version, PALM commit, input-file hash, and command line for each case.
4. Change the plotting scripts so that they read a specified file instead of defaulting to the newest Exodus file in the directory.
5. Add checks for nonphysical values:
   - $0\le c\le1$;
   - $0\le\phi\le1$;
   - $0\le P<1$;
   - $k>0$;
   - $T>0$;
   - The phase fractions sum to 1.

### Pass Criteria

- Repeating the same input gives consistent results;
- Postprocessing scripts do not depend on manual edits;
- Error metrics can be generated by a single command;
- The mass-conservation error should first be kept within a relative magnitude of $10^{-8}$; if solver tolerances prevent that, record how the error changes with tolerance.

## Stage 1: Complete the Mohanty 2011 Single-Phase Homogeneous-Alloy Benchmark

### Subtask 1A: Equation-Level Verification

Build a small manufactured solution or a simplified case and verify separately:

- The chemical-diffusion part of `SplitCHParsed`/`SplitCHWRes`;
- The thermal-migration part of `ADThermalMigration`;
- The total flux after the two terms are superimposed;
- Temperature reversal, reversal of the thermal-migration sign, and the zero-temperature-gradient limit.

### Subtask 1B: Paper Figs. 5–7

Run, in order:

1. Constant $\beta_i$, constant $Q_i^*$;
2. Constant $\beta_i$, temperature-dependent $Q_i^*$;
3. Temperature-dependent $\beta_i$, constant $Q_i^*$;
4. Transients at 5, 10, and 30 days;
5. Hot-end concentration versus time.

If the complete data for temperature-dependent $Q_i^*$ cannot be obtained from the paper, mark that sub-item as “not uniquely reproducible” and do not construct a function independently.

> Revision R13. Two further items about Mohanty's set-up that the student must record rather than discover mid-run. (i) For the "constant $\beta_i$" cases of Fig. 5 the paper does not state at which temperature $\beta_i$ is evaluated; run them with $\beta_i(1225\ \mathrm K)$ (mid-domain) as primary and $\beta_i(1050)$, $\beta_i(1400)$ as brackets, and report the spread. (ii) Mohanty's Eq. (8) prescribes both Dirichlet temperatures at the two ends *and* zero heat flux; on a 1-D interval this is over-determined. Use only the two Dirichlet values, which give the linear $T(x)$ that the paper plots; note the over-determination in the README.

### Verification Quantities

- Migration direction: Zr migrates toward the hot end;
- Composition changes at the hot and cold ends at 30 days;
- Rapid migration in the first 5 days and slower migration later;
- Endpoint-value error and full-profile error;
- Mass conservation;
- Mesh convergence on 150/300/600 elements;
- Convergence on three sets of maximum time steps.

## Stage 2: Mohanty 2011 Diffusion Couple and Flux Decomposition

### Objective

Reproduce paper Fig. 8 and confirm the competition between thermal migration and chemical diffusion.

### Cases

1. Isothermal diffusion couple at 1225 K;
2. Temperature gradient in the same direction as the initial Zr concentration gradient;
3. Temperature gradient opposite to the initial Zr concentration gradient.

The initial compositions are 39 at.% and 22.5 at.% Zr. A step initial condition should first use a narrow, smooth hyperbolic-tangent transition, to avoid mesh dependence caused by an initially discontinuous finite-element field; the transition width must be smaller than the paper length scale and must be subjected to a sensitivity analysis.

### Required Output

Define and output

$$
\mathbf J_c=-M_c\nabla w,
\qquad
\mathbf J_T=M_T\frac{\nabla T}{T},
\qquad
\mathbf J_{\mathrm{total}}=\mathbf J_c+\mathbf J_T.
$$

Verify the order of magnitude of $\lvert J_T/J_c\rvert$ stated in the paper, and check its spatial and temporal variation. Do not compare only the final concentration curves.

## Stage 3: Irradiation-Free KKS Multiphase Skeleton

### Objective

Before adding irradiation enhancement, independently verify the $\alpha/\beta/\gamma$ thermodynamics, the KKS constraints, and Allen–Cahn phase evolution.

### Recommended Variables

- `c`: total Zr atom fraction, the conserved variable;
- `w`: split CH chemical potential;
- `phi`: $\gamma$-phase order parameter;
- `c_ab`: Zr concentration inside the $\alpha\beta$ mixed phase;
- `c_g`: Zr concentration inside the $\gamma$ phase;
- `T`: a prescribed AuxVariable at this stage;
- `eta_a`, `eta_b`, `eta_g`: AuxVariables or material properties for output.

### MOOSE Object Mapping

Prefer to reuse the KKS objects in the phase-field module:

- `KKSPhaseConcentration`
- `KKSPhaseChemicalPotential`
- `KKSSplitCHCRes`
- `KKSACBulkF`
- `KKSACBulkC`
- `KKSGlobalFreeEnergy` (an AuxKernel for output of the total free energy density, not a Kernel; listed here for completeness)
- `TimeDerivative`
- `ACInterface`

References:

- `moose/modules/phase_field/test/tests/KKS_system/two_phase.i`
- `moose/modules/phase_field/test/tests/KKS_system/kks_example_split.i`
- `moose/modules/phase_field/examples/kim-kim-suzuki/`

The free energy may first be implemented with `DerivativeParsedMaterial`, but the domain issues of $\ln c$ and $\ln(1-c)$ should be explicit tests. Do not hide a failed intra-phase concentration solve behind an excessively large numerical cutoff.

`KKSSplitCHCRes` enforces the KKS chemical-potential relation and does not automatically add Wen’s $\kappa_c\nabla^2c$. If the Wen branch uses a nonzero $\kappa_c$, the gradient term must be added explicitly with `MatDiffusion` or an equivalent Kernel; the strict Jung branch does not add that term. If the definition of $2\kappa_c$ in Mohanty equation (1) is retained, the coefficient passed to the split Kernel should be twice the paper $\kappa_c$; in the strict $\kappa_c=0$ branch that factor does not affect the result.

### Stepwise Verification

1. Compare the single-point free energy and its first and second derivatives with independent Python/SymPy results;
2. Given $c,T,\phi$, check the residuals of the two KKS constraints;
3. Isothermal two-phase equilibrium: check equality of chemical potentials;
4. Interface equilibrium without thermal migration;
5. Locations of the three phase regions under a prescribed temperature field;
6. Check $\eta_\alpha+\eta_\beta+\eta_\gamma=1$;
7. Check how the interface width varies with $\kappa_\phi/\omega$;
8. Compare the computed phase diagram with the assessed U–Zr diagram (see R11 below).

> Revision R11. Step 8 is new and is the decisive thermodynamic test of the stage. Using the same `f_ab` and `f_gamma` parsed materials (through a Python re-implementation checked against them in step 1), perform the common-tangent construction at 900, 935, 950, 1000 and 1050 K and plot the resulting $(\alpha\text{ or }\beta)+\gamma$ two-phase boundaries over the assessed U–Zr phase diagram (Sheldon & Peterson 1989; Quaini et al. 2018 CALPHAD). Acceptance: the $\gamma$-side boundary of the $(\alpha,\beta)+\gamma$ field in the 900–1000 K window must be Zr-rich (of order 60–75 at.% Zr in the assessments) and agree with the assessed value to within 5 at.%; the solubility of Zr in $\alpha$-U and $\beta$-U must be at most a few at.%, as in the assessments. Note that the model has no $\delta\text{-}\mathrm{UZr_2}$ phase, so agreement is expected only above the $\delta$ formation temperature (about 890 K); DP-81's 900 K surface sits just above it. Jung reports over-prediction of the surface $\alpha$ zone; a phase diagram with too wide a two-phase field is the most likely cause, and it is far cheaper to detect it here than in Stage 8. Also record the position of the $\gamma_1/\gamma_2$ miscibility gap implied by `f_gamma`, since it controls the "Zr well" of Wen Fig. 2–3.

> Revision R10 (Stage 3 consequence). Two input-file switches are introduced here and carried into every later KKS stage: `interface = {wen, jung}` selecting $(\kappa_\gamma,\omega,L)$ from Wen Table 1 or from Jung's $(l,\sigma,L)$, and `T_mode = {local, fixed_1010}` for the $\gamma$ free energy. The Stage 3 mesh study must be run with the narrower Wen interface ($\approx40\ \mathrm{\mu m}$, needing $\le8\ \mathrm{\mu m}$ elements or adaptivity) because it sets the more demanding requirement.

### Important Ambiguity

The interpolation functions and barrier functions in Wen 2022 and Jung 2025 may differ in typesetting. In particular, Jung equation (8) is printed as $g(\phi)=\phi^2(1-\phi^2)$, while the common double-well form is $\phi^2(1-\phi)^2$. Before implementation:

- Build one version verbatim from the original text;
- Build a reference version from the model it cites;
- Check whether both versions have the expected minima at $\phi=0,1$;
- If the original form cannot produce a reasonable double well, record it as a typesetting or definition ambiguity in the paper.

## Stage 4: Wen 2022 Irradiation-Enhanced Diffusion

### Objective

Add irradiation only through the effective-mobility method of the Wen paper, without introducing extra defect equations.

### Implementation Order

1. Compute the thermal-equilibrium vacancy concentration $C_v^e(T)$;
2. Use the paper’s prescribed $C_v^r$;
3. Compute $\xi(T)$;
4. Set `chemical_mobility_irradiated = xi * chemical_mobility`;
5. Compute the thermal mobility from paper equations (29)–(31). Provide a switch `xi_target = {chemical, both}`: `chemical` is Wen's literal model ($M_c^I=\xi M_c$, $M_T$ unchanged); `both` also sets $M_T^I=\xi M_T$;
6. Compare irradiated and unirradiated results, and `chemical` against `both`.

> Revision R5 (Stage 4). Item 5 is changed from "do not multiply unless re-checked" to a mandatory two-branch comparison. The paper's equations are unambiguous ($\xi$ on $M_c$ only), so the literal branch is settled; the issue is physical. Vacancy-mediated diffusion enhancement multiplies the atomic mobilities $\beta_i$, and $M_T$ is built from the same $\beta_i$, so the physically consistent model is `both`. The two branches predict opposite trends for the *steady-state* Zr redistribution (unchanged for `both`, reduced by $1/\xi$ for `chemical`) and identical trends for the early transient. Because Wen stops at 1.5 days, Fig. 3 cannot discriminate between them; the Stage 4 report must show both curves at $1.3\times10^5\ \mathrm s$ and at a time long enough for the center concentration to plateau, and it must state that the paper's central claim (RED replaces an artificial diffusivity enlargement) is not tested by the paper's own run length (see §2.2, R5).

### Suggested Material Properties

- `thermal_vacancy_concentration`
- `irradiation_vacancy_concentration`
- `irradiation_enhancement_factor`
- `chemical_mobility_unirradiated`
- `chemical_mobility_irradiated`

These relations may first be implemented with `ADParsedMaterial` or an AD parsed material. A C++ Material does not need to be written immediately. Add `ADU_ZrIrradiationMobilityMaterial` only when piecewise functions, numerical safeguards, or performance become a problem.

### Paper Reproduction Settings

- One-dimensional *radial* domain, $0\le r\le2170\ \mathrm{\mu m}$, solved with `coord_type = RZ` (R7); the $r=0$ end carries the natural zero-flux condition of the axisymmetric operator and needs no explicit BC;
- 2170 finite-difference elements is the original paper setting; MOOSE should perform an independent mesh-convergence study and should not mechanically require the same number of elements;
- About 988 K at the center and 900 K at the surface;
- Initial Zr atom fraction of about 0.225;
- Zero constituent flux and zero phase-field flux;
- Paper explicit time step of 1 s and total time of $1.3\times10^5$ s;
- MOOSE may use implicit adaptive time stepping, but it must be shown that the result does not depend on a larger time step.

### Target Figures

- Wen Fig. 2: without irradiation enhancement;
- Wen Fig. 3: with irradiation enhancement;
- Wen Fig. 4: comparison with other models and experimental data;
- Wen Fig. 5: sensitivity to $\kappa_\gamma$ and $\omega$;
- Wen Table 2: interfacial energy and macroscopic interface width.

### Pass Criteria

- The irradiated model clearly enhances kinetics in the low-temperature region;
- The center Zr-rich region, the intermediate “Zr well,” and the outer depleted region described in the paper can form;
- The order of the three phase regions is consistent with the paper;
- The KKS constraints, mass conservation, and phase-fraction normalization pass;
- One-parameter sensitivity analyses are performed separately for $\xi$, $\kappa_\gamma$, and $\omega$, plus the `xi_target` branch comparison (R5);
- The sign of the Zr thermal flux in the $\gamma$ phase is verified against Wen's convention (R4): Zr enriches at the hot center;
- Parameters given in the literature, cited parameters, and trial-fit parameters are clearly distinguished.

## Stage 5: Jung 2025 Heat-Conduction Baseline

### Objective

First build a fully coupled $c$–$\phi$–$T$ baseline without porosity or sodium infiltration, then add thermal-conductivity corrections.

### Heat-Equation Objects

Suggested objects (primary, quasi-static form):

- `ADHeatConduction`
- `ADBodyForce` or an equivalent AD heat-source Kernel
- A surface `ADDirichletBC` at $r=r_0$ (the paper's surface temperature); no BC at $r=0$ in RZ

Optional, comparison branch only (`thermal = transient`):

- `ADHeatConductionTimeDerivative` with `density` and `specific_heat` tagged `cited`

Material properties should include at least:

- `thermal_conductivity`
- `volumetric_heat_source`
- (`transient` branch only) `density`, `specific_heat`

> Revision R8. Revision 1 listed the transient kernel first and then worried about missing $\rho c_p$. The order is reversed: the quasi-static equation is the physically correct primary model because the slug's thermal relaxation time ($\approx0.4\ \mathrm s$) is negligible against every other time scale in the problem, including the smallest phase-field time step the student will use. Solving $0=\nabla\cdot(k\nabla T)+\dot q$ monolithically with $c,w,\phi$ at each step introduces no new parameters and no new stiffness. The transient branch is kept only to demonstrate this equivalence once; if it is run, use $\rho=15.8$–$16.0\ \mathrm{g\,cm^{-3}}$ and $c_p=150$–$200\ \mathrm{J\,kg^{-1}\,K^{-1}}$ (U–10Zr, 900–1000 K) as `cited` values and show that $T(r)$ differs from the steady branch by less than 0.1 K after the first step. Note also that with a quasi-static $T$ the Jacobian block $\partial R_T/\partial c$ through $k(c)$ is what couples the systems; the AD kernels supply it automatically.

### Conversion from Atom Fraction to Weight Fraction

The Jung thermal-conductivity relation uses the Zr weight fraction, while the phase-field variable is the Zr atom fraction. Use

$$
W_{\mathrm{Zr}}
=
\frac{cM_{\mathrm{Zr}}}
{cM_{\mathrm{Zr}}+(1-c)M_U},
$$

$$
W_U=1-W_{\mathrm{Zr}}.
$$

The extracted text of paper equations (34)–(35) contains a circular definition or unclear typesetting. The code must use an independently verifiable atom-fraction-to-weight-fraction conversion and explain it in the README.

### Compositional Gradient Term in the Strict Jung Model

Jung equations (5)–(11) give only the phase-order-parameter gradient energy

$$
\frac{\kappa_\phi}{2}\lvert\nabla\phi\rvert^2,
$$

and do not report $\kappa_c\lvert\nabla c\rvert^2/2$. Therefore:

- $\kappa_c=0$ should be set in `jung_literal`;
- If Wen’s nonzero $\kappa_c$ is carried over, it must be named `wen_inherited` or a sensitivity branch;
- Both branches may use the split CH software structure, but when $\kappa_c=0$ the equation physically reduces to second-order nonlinear diffusion.

### Substeps

1. Verify the steady *cylindrical* analytical solution $T(r)=T_s+\dot q(r_0^2-r^2)/(4k)$ for constant $k$ in `coord_type = RZ`, and check that the integrated heat source equals LHGR per unit length (this verifies the $2\pi r$ weighting);
2. Set only $k=k(T)$;
3. Set $k=k(T,c)$;
4. Turn on one-way $c$–$\phi$–$T$ coupling: temperature affects the phase field, but the phase field does not feed back into thermal conductivity;
5. Turn on two-way coupling: $c$ affects $k$, and $k$ affects $T$;
6. Check the difference between a monolithic solve and a split solve.

The paper does not report the total simulation time, yet it lets burnup increase linearly over the total simulation time and uses Hirschhorn’s optimized mobility to compress a real irradiation process of about 284 EFPD. Jung results can therefore be used to reproduce end-state profiles, but they cannot be used to verify real-time kinetics until the authors’ input files or supplementary data are obtained.

### Geometric Notes

The paper treats the left end of the one-dimensional interval as the fuel center and the right end as the surface. Two clearly distinguished inputs are created:

- `paper_1d.i`: axisymmetric radial formulation, `coord_type = RZ`, primary;
- `slab_check.i`: Cartesian slab operator on the same interval, sensitivity only.

> Revision R7. Revision 1 had these two roles reversed and called the cylindrical result "a geometric correction that cannot replace the reproduction". That is backwards. The fuel slug is a cylinder; Wen says "axisymmetric"; Jung normalizes LHGR by $\pi r_0^2$ and reports a 90 K center-to-surface drop that a slab cannot produce with the same $k$ and $\dot q$ (a slab gives twice the drop, §2.3). If either paper's code actually used a slab operator, that is a defect of the paper, not a target for reproduction, and the discrepancy would show up as a factor-2 error in $\Delta T$ that `slab_check.i` will quantify. The $c(r)$ profile is likewise affected because the outward area grows with $r$: in RZ a given radial flux moves less concentration per unit volume near the surface than near the center. Mass-conservation postprocessors must use the RZ-weighted integrals.

## Stage 6: Jung 2025 Porosity Coupling

### Objective

Map phase fractions to porosity with the paper’s algebraic model, and correct only the thermal conductivity.

### Implementation

Porosity should be an AD material property, not a nonlinear variable:

$$
P_{\mathrm{phase}}
=0.4\eta_\gamma+0.13(\eta_\alpha+\eta_\beta).
$$

Then multiply by the burnup growth function:

$$
P(B)=s_P(B)P_{\mathrm{phase}},
$$

where

$$
s_P(B)=
\begin{cases}
B/0.01, & 0\le B<0.01,\\
1, & B\ge0.01.
\end{cases}
$$

Here $B$ should consistently be a dimensionless fraction, so that 1% is written as 0.01 rather than 1.

### Comparison Cases

- Constant $P=0$;
- Constant $P=0.1$;
- Constant $P=0.2$;
- Constant $P=0.3$;
- Baseline porosity that depends on phase fractions.

### Target Figures and Metrics

Reproduce the following from Jung Fig. 5:

- Zr concentration;
- Temperature;
- Thermal conductivity;
- $\gamma$ phase fraction.

Additional output:

- Porosity profile;
- Zone I/II/III boundary locations;
- Center temperature versus porosity;
- Center Zr concentration versus porosity.

Order-of-magnitude checks given in the paper abstract include: when porosity is neglected, the center Zr concentration decreases by about 0.17 at.%, and the center temperature decreases by about 36 K. Specific comparisons must be based on the definition of the reproduced baseline.

## Stage 7: Jung 2025 Sodium Infiltration

### Objective

Add sodium’s correction to the effective thermal conductivity of the pores according to paper equations (31)–(32).

### Implementation Order

1. Independently verify $k_{\mathrm{Na}}(T)$;
2. Resolve the typographic ambiguity of Eq. (31) with the three physical checks of R9 (§2.3): $P_c(P_{\mathrm{Na}}=0)=1$, $P_c\ge1$ when $k_{\mathrm{Na}}>k_0$, and agreement with a Maxwell–Eucken-type bound when all pores are filled; record the rejected reading;
3. Check the value of $k_{\mathrm{Na}}/k_0$ over the full temperature and composition ranges (expected $\approx2$ at 900 K);
4. Add 5%, 25%, and 50% sodium infiltration;
5. Add $P_{\mathrm{Na}}(B)$ that grows with burnup;
6. Turn on all feedback.

The sodium-infiltration percentage must be explicit about whether it is a fraction from 0 to 1 or a percentage from 0 to 100. Internally, the code should consistently use 0–1 and convert at the input layer.

The paper also does not make clear whether $P_{\mathrm{Na}}$ means “the fraction of pores filled by sodium” or “the absolute volume fraction of sodium in the total fuel volume.” Two explicit branches should be created:

$$
P_{\mathrm{Na}}^{\mathrm{literal}}=f_{\mathrm{Na}},
$$

$$
P_{\mathrm{Na}}^{\mathrm{pore}}=f_{\mathrm{Na}}P.
$$

The main reproduction should first use the literal definition in paper equation (31). The other branch is for uncertainty analysis. The definition must not be switched automatically according to the local porosity inside the same input file.

> Revision R9 (Stage 7). The two definitions are not equally plausible. Jung's burnup ramp raises the sodium quantity from 1 % to 1.5 % burnup *after* the porosity has saturated, and reports cases labelled 5 %, 25 %, 50 % "sodium infiltration"; a 50 % *absolute* volume fraction of sodium in a fuel with $P\le0.4$ is impossible, so those labels must be fractions of the pore volume, i.e. the `pore` reading, whatever the printed formula appears to say. Run `literal` only up to the value where $P_{\mathrm{Na}}^{\mathrm{literal}}\le P$ everywhere, and report explicitly which reading reproduces the paper's stated 0.7 at.% and 20 K center shifts for the 50 % case. That comparison, not the printed equation, settles the definition.

### Target Figures

Reproduce the following from Jung Fig. 6:

- Zr concentration;
- Temperature;
- Thermal conductivity;
- $\gamma$ phase fraction.

The order-of-magnitude check given in the paper is: relative to Baseline, 50% sodium infiltration lowers the center Zr concentration by about 0.7 at.% and lowers the center temperature by about 20 K. The whole profile should be compared at the same time; the center point alone should not be fitted.

## Stage 8: DP-81 and DP-11 Validation

### DP-81

- Radius/one-dimensional length: 2.17 mm;
- Surface temperature: 900 K;
- LHGR: 24 kW/m;
- Target burnup: 3.6 at.%;
- Initial Zr: about 23 at.%;
- Initial $\gamma$ phase fraction: 0.5 at the center, 0.4 at the surface;
- Initial temperature: 990 K at the center, 900 K at the surface;
- Baseline sodium infiltration: 5%.

### DP-11

- Surface temperature: 910 K;
- LHGR: 25 kW/m;
- Target burnup: 7.7 at.%;
- The remaining missing initial conditions should be traced from paper figures, cited reports, or supplementary material, and should not be copied directly from DP-81.

### Validation Order

1. Temperature field;
2. Thermal-conductivity field;
3. Three-phase fractions and phase-region boundaries;
4. Zr concentration;
5. Porosity;
6. Evolution of center quantities with time/burnup.

PIE data should be digitized, and the error should be computed after the experimental scatter and the simulation are interpolated onto the same normalized radius $r/r_0$.

> Revision R7 (Stage 8). All Stage 8 runs are in `coord_type = RZ`. The first validation item, the temperature field, is a direct check of the geometry: with DP-81's 24 kW/m and the composition-dependent $k$, the cylindrical solution gives a center-to-surface drop close to the 90 K Jung reports; the slab gives about twice that. A slab run that "matches" the concentration profile does so only by compensating a wrong $\Delta T$ with a wrong effective time scale, and must not be presented as a validation.

> Revision R5 (Stage 8). Validation against DP-81 and DP-11 tests the *end state* after several hundred EFPD, while the models are run for compressed simulated times with optimized mobilities (Hirschhorn) or with a 1.5-day horizon (Wen). State explicitly in the Stage 8 report which quantities are legitimately compared (steady or near-steady profiles: $T$, $k$, zone boundaries, near-steady $c$) and which are not (any time-resolved quantity). The "evolution of center quantities with time/burnup" item is a model-consistency output, not a validation against experiment, because no time-resolved PIE data exist for these pins.

## Stage 9: Two-Dimensional Extension

### Objective

Without adding physics beyond the papers, check the effects of the one-dimensional assumption, geometric discretization, and nonuniform boundaries on the results.

### Order

1. Reproduce the one-dimensional RZ solution on a two-dimensional $r$–$z$ strip with uniform axial conditions, for code regression (the strip and the 1-D RZ mesh must agree to solver tolerance, since they discretize the same operator);
2. Extend the axial length to a representative slug segment and keep uniform axial conditions;
3. Allow axial power or surface temperature to vary (e.g. a chopped-cosine LHGR and a coolant temperature rise along $z$);
4. Compare the centerline, results at different axial positions, and volume-averaged results;
5. Check two-dimensional mesh anisotropy and parallel scalability.

The first acceptance criterion for the two-dimensional results is: when the boundary and initial values are uniform along the second dimension, the two-dimensional solution should reduce to the one-dimensional solution.

## 6. Suggested New PALM Code

## 6.1 Required Custom Objects

### `ADThermalMigration`

Purpose: directly discretize

$$
-\nabla\cdot\left(M_T\frac{\nabla T}{T}\right)
$$

and retain the automatic-differentiation Jacobian for all material dependencies.

Minimum parameters:

- `variable`: the variable of the split mass-balance equation;
- `T`: the temperature variable;
- `thermal_mobility`: the AD material-property name.

Tests:

- Manufactured solution with constant $M_T$;
- Jacobian test for $M_T(c,T)$;
- Temperature reversal;
- Comparison with `SoretDiffusion` under a constant-coefficient mapping.

## 6.2 Prefer Parsed Materials, and Customize Only When Necessary

The following relations should preferably be implemented with `ADParsedMaterial`, `DerivativeParsedMaterial`, `ParsedFunction`, or `PiecewiseLinear`:

- $\alpha\beta/\gamma$ free energies;
- Chemical mobility of each phase;
- Thermal mobility of each phase;
- $\xi(T)$;
- Phase fractions;
- Porosity;
- $k_0(T,c)$;
- $k_{\mathrm{Na}}(T)$;
- $P_c$;
- Burnup and LHGR histories.

Write a C++ Material only when one of the following occurs:

- A complex piecewise expression is hard to review;
- The same formula is repeated across multiple cases;
- Unit conversion and range checks need to be centralized;
- The parsed expression causes a clear performance problem;
- Clear error messages and parameter-validity checks are needed.

At that point, consider:

- `ADUZrKKSFreeEnergyMaterial`
- `ADUZrMobilityMaterial`
- `ADUZrThermalConductivityMaterial`

These classes should not all be implemented at once before the equations have passed single-point verification.

## 6.3 Upstream MOOSE Code That Should Not Be Modified

Do not modify directly:

- `moose/modules/phase_field/src/kernels/SoretDiffusion.C`
- MOOSE KKS Kernels;
- MOOSE HeatConduction Kernels.

Upstream changes increase the cost of upgrades, regression, and result traceability. PALM custom objects can express the paper equations more accurately while keeping the MOOSE checkout updatable.

PALM currently has no custom C++ physics objects. `HEAT_TRANSFER` and `PHASE_FIELD` are already enabled, and `POROUS_FLOW` is not enabled. Jung’s porosity is only an algebraic material property, so PorousFlow does not need to be enabled for it; module dependencies should be reassessed only if pore flow is actually solved in the future.

## 7. Solution Strategy

### 7.1 First Stage

The current `PJFNK + LU` can continue to be used for the Mohanty comparison. After the new AD Kernel is added, a full Newton solve should be tried:

- Small problems: `NEWTON + LU`;
- Larger problems: PJFNK or JFNK with field-split preconditioning;
- Variable scaling for $c,w,\phi,c_{\alpha\beta},c_\gamma,T$.

### 7.2 Strongly Coupled Stage

The Jung model should first use a monolithic solve, because temperature, phase fractions, and composition form clear feedback. With the quasi-static heat equation (R8) the thermal variable adds one elliptic block whose off-diagonal coupling to $c$ enters only through $k(c,\phi,T)$; this is mild and a monolithic Newton solve with `NEWTON + LU` in 1-D is expected to converge without difficulty. If the monolithic system is hard to converge, then try:

1. At each time step, first solve the thermal steady state, then solve the phase field;
2. A thermal–phase-field Picard outer iteration;
3. Finally, compare the difference between the split scheme and the monolithic scheme.

A split scheme should not be assumed equivalent merely because it converges more easily. Check whether the two schemes approach the same solution as the time step is reduced.

### 7.3 Numerical Safeguards

Allowed numerical safeguards include:

- A recorded smoothing of the initial step;
- Adaptive truncation of the time step;
- Reasonable scaling of nonlinear variables;
- Variable bounds or a line search to keep $c$ inside $(0,1)$.

Not recommended:

- Truncating the free energy without a record;
- Using a very large gradient-energy coefficient to hide an insufficient mesh;
- Using a mobility scale to fit both the curve shape and the time scale at once;
- Freely adjusting several highly correlated parameters in the same fit.

## 8. Verification and Uncertainty Quantification

### 8.1 Conservation

For no-flux boundaries, compute

$$
\epsilon_m(t)
=
\frac{\left|\int_\Omega c(t)\,\mathrm dV
-\int_\Omega c(0)\,\mathrm dV\right|}
{\int_\Omega c(0)\,\mathrm dV}.
$$

The integration measure for the one-dimensional slab, the cylinder, and the two-dimensional model must be consistent with the geometry. In MOOSE, `coord_type = RZ` makes every `ElementIntegral*` postprocessor use $2\pi r\,\mathrm dr\,\mathrm dz$ automatically; a slab input that reports a "conserved" $\int c\,\mathrm dx$ says nothing about conservation in the cylinder. For the split CH formulation, also record the boundary flux $\int_{\partial\Omega}\mathbf J\cdot\mathbf n\,\mathrm dS$ with a `SideIntegral*` postprocessor: with `ADThermalMigration` the natural BC of the split equation is zero *total* flux, $\mathbf J_c+\mathbf J_T=0$, which is the physically required condition and differs from Mohanty's $c_x=0$.

### 8.2 Mesh Convergence

Each key stage should use at least three meshes. Compare:

- Center and surface $c$;
- Maximum/minimum $c$;
- Zone boundaries;
- Center temperature;
- Full-profile $L_2$ difference.

Interface models must also ensure enough elements inside each diffuse interface. Mesh convergence cannot be judged from global averages alone.

### 8.3 Time-Step Convergence

Using the same spatial mesh and solver tolerances, halve the maximum time step at least twice in succession. For burnup-driven piecewise functions, the time step must resolve the turning points near 1% and 1.5% burnup.

### 8.4 Error Against Paper Curves

For a digitized paper curve $y_i^{\mathrm{ref}}$ and a simulated curve $y_i$, report

$$
E_{L_2}
=
\sqrt{
\frac{\sum_i(y_i-y_i^{\mathrm{ref}})^2}
{\sum_i(y_i^{\mathrm{ref}})^2}
},
$$

$$
E_\infty=\max_i\lvert y_i-y_i^{\mathrm{ref}}\rvert.
$$

Also report:

- Number of digitized points;
- Interpolation method;
- Coordinate-normalization method;
- Uncertainty of the image reading;
- Whether parameters were fitted.

### 8.5 Parameter Classification

Label each parameter as one of the following:

- `direct`: given directly by the paper;
- `cited`: taken from literature cited by the paper;
- `digitized`: read from a figure;
- `derived`: converted from known quantities;
- `calibrated`: adjusted to fit results;
- `assumed`: an assumption used when the value is missing.

Only `direct`, `cited`, and verifiable `derived` parameters can be called unfitted inputs.

## 9. Key Gaps and Risks in the Papers

### Mohanty 2011

- Equation (3) as printed, with no composition weight inside the bracket, is the correct Darken/Jaffe–Shewmon thermotransport form; $M_Q/M_c$ is finite at both composition limits (R2). It is *not* a gap in the paper; the gap was in the student's reading;
- The heat-of-transport convention (signed $\widetilde Q_i^*$, minus sign in Eq. (3), per-atom units) differs from Wen's and Jung's (positive magnitudes, plus sign, per-mole units); $\widetilde Q_U^*$ differs between the papers by a factor of about seven (R1, R4);
- The paper uses the full $\nabla(\partial f/\partial c)$ as chemical driving force in a temperature field; the split into isothermal and thermal parts, and hence the exact meaning of its $\widetilde Q_i^*$, is not discussed (R3);
- The paper states that $\kappa_c=0$ in the single phase; the nonzero gradient energy in the current case is extra numerical regularization;
- The atomic-mobility prefactors are not fully tabulated in the main text and must be read from figures or traced through citations;
- The temperature at which the "constant $\beta_i$" of Fig. 5 is evaluated is not stated (R13);
- The free energy comes from a commercial database and cannot be uniquely recovered from the paper text alone;
- The text says 1405 K, while the figure axis may show 1400 K;
- Equation (8) prescribes both Dirichlet temperatures and zero heat flux on a 1-D interval, which is over-determined (R13);
- $c_x=0,c_{xxx}=0$ in paper equation (6) is not strictly equivalent to a total no-flux boundary that includes the thermal-migration term;
- The paper's own diffusion time $L^2/D_U\approx0.3\ \mathrm{day}$ at 1400 K is inconsistent with the continued profile evolution it shows between 10 and 30 days; the paper's kinetics are internally slower than its stated $\beta_i$ by a factor of order 10–30, which is the likely origin of the student's `scale` and cannot be resolved from the final curves alone (R12).

### Wen 2022

- $C_v^r$ comes from an external rate-theory code and is not a derivation that can be repeated inside the paper;
- The irradiation-enhancement factor lacks a complete dependence on irradiation conditions;
- The $\alpha/\beta$ phase-fraction and mobility interpolation weights are opposite to the paper’s phase-region text;
- $L_0$ is not given, and $\kappa_c$ is given only as an order of magnitude;
- Some free-energy expressions use a fixed temperature or have unclear typesetting;
- The abstract says the effect of irradiation on thermal diffusion is discussed, but the governing equation only sets $M_c^I=\xi M_c$; applying $\xi$ to $M_c$ alone is inconsistent with a vacancy mechanism and *reduces* the steady-state Soret redistribution by $1/\xi$ (R5);
- The simulated time of $1.3\times10^5\ \mathrm s$ (1.5 days) is two orders of magnitude shorter than the DP-81 irradiation; the comparison with PIE is therefore a transient matched to an end state, and the claim that RED removes the need for an artificial diffusivity enlargement is not supported by the run length (R5);
- $C_v^r$ was computed for 30 kW/m while DP-81 ran at 24 kW/m;
- The domain is stated to be axisymmetric; whether the code's 1-D operator was cylindrical is not verifiable from the text (R7);
- The time discretizations of explicit finite differences and MOOSE implicit finite elements differ;
- The phase-interface parameters involve a clear trial-and-fit process;
- Zr enrichment near the fuel surface cannot be explained by this model.

### Jung 2025

- Hirschhorn’s optimized mobility is used to compress a real time of 284 EFPD, so the strict time scale cannot be interpreted directly;
- The total simulation time is not reported, so the slope of burnup versus time cannot be uniquely recovered;
- Transient thermal parameters such as $\rho$ and $c_p$ are absent, but they are also irrelevant: the thermal relaxation time is of order 1 s, so the heat equation is quasi-static (R8);
- $\beta_0$, $\widetilde Q_U^*$, $\widetilde Q_{\mathrm{Zr}}^*$, and the optimized mobility scale are incomplete;
- The $\gamma$ free energy is evaluated at a fixed $T_{\mathrm{ref}}=1010\ \mathrm K$ except for the SGTE end-members, so it is a different function from Wen's (R10);
- The interface width $l=0.5\ \mathrm{mm}$ is 23 % of the fuel radius and an order of magnitude wider than Wen's; zone boundaries are correspondingly smeared (R10);
- The thermal-conductivity coefficient numbering and typesetting may be misaligned;
- The prefactor of the sodium correction factor $P_c$ in Eq. (31) is typographically ambiguous; only one reading satisfies $P_c(0)=1$ and $P_c\ge1$ for $k_{\mathrm{Na}}>k_0$ (R9);
- The formula text for converting atom fraction to weight fraction is ambiguous;
- The $\alpha/\beta$ switching direction (same inversion as Wen; figures show $\alpha$ at the surface), the double-well function, and the definition of the sodium-infiltration amount are all ambiguous (R6, R9);
- The paper does not report a compositional gradient energy, so the strict branch should take $\kappa_c=0$;
- The porosity model is an algebraic average over phase regions, not an evolution of pore structure;
- Both LHGR and burnup use a simplified linear history;
- The paper acknowledges deviations in the surface $\alpha$ phase and the surface concentration; a too-wide two-phase field in the model's phase diagram is the first suspect (R11);
- The heat source is normalized by $\pi r_0^2$, which presupposes cylindrical geometry; the reported 90 K drop is consistent with the cylindrical solution and not with a slab (R7).

These gaps should be part of the reproduction report, rather than being hidden by unexplained parameter adjustments.

## 10. Milestones and Completion Criteria

### M1: Single-Phase Thermal Migration Is Credible

- Trends in Mohanty Figs. 5–7 and the diffusion couple are reproduced;
- The original-equation and calibrated-equation cases are clearly separated;
- Mass, mesh, and time-step checks pass;
- The `ADThermalMigration` tests pass.

### M2: Multiphase KKS Is Credible

- Single-point free energy and derivatives pass;
- KKS constraint residuals converge;
- Phase fractions are normalized;
- The computed common-tangent phase boundaries agree with the assessed U–Zr diagram to within 5 at.% in 900–1050 K (R11);
- Interface width is mesh-independent for both the `wen` and `jung` interface parameter sets (R10);
- The unirradiated baseline is stable, and `phase_consistent` places $\alpha$ at the cold surface (R6).

### M3: Wen Irradiation Enhancement Is Reproduced

- The difference between unirradiated and irradiated cases is reproduced in `coord_type = RZ`;
- Three phase regions and the Zr well appear;
- Sensitivity to $\xi$ and the interface parameters is completed, including the `xi_target = both` branch and its long-time steady state (R5);
- The Zr thermal-flux sign follows Wen's convention and Zr enriches at the center (R4);
- No defect equations outside the paper are introduced.

### M4: The Jung Thermal-Feedback Baseline Is Credible

- The cylindrical analytical heat-conduction benchmark passes, including LHGR recovery from the RZ-weighted source integral (R7);
- The quasi-static and transient thermal branches agree to $<0.1\ \mathrm K$ (R8);
- One-way and two-way coupling of $k(T,c)$ pass layer by layer;
- Temperature, phase fractions, and composition of the DP-81 Baseline reach quantifiable agreement.

### M5: Porosity and Sodium Infiltration Are Reproduced

- All four profile types in Jung Figs. 5 and 6 are completed;
- The magnitudes of the center-temperature and center-Zr changes are consistent;
- The burnup piecewise function has passed a time-step convergence check;
- Porosity is coupled only through thermal conductivity.

### M6: The Two-Dimensional Framework Is Complete

- It reduces to the one-dimensional solution when the second dimension is uniform;
- The effect of axisymmetric geometry is quantified independently;
- The two-dimensional calculation preserves mass conservation and mesh convergence;
- The one-dimensional paper-reproduction results remain the primary validation benchmark.

## 11. Recommended Immediate Next Steps

1. Write `UNITS.md` for `s1_mohanty` with the per-atom/per-mole audit table of R1, and correct the README: remove the "divergence" argument, correct `scale = 0.05` to `0.035`, and fix the `kappa_c` comment.
2. Freeze the current `s1_mohanty.i` as the calibrated version, with its three modifications (scale, composition weight, end-member omission) listed in the header.
3. Build the original-equation version of the Mohanty paper: no mobility scale, no composition weight, end-member terms retained, $\kappa_c=0$ (R2, R3).
4. Implement and test `ADThermalMigration` in PALM.
5. Test the ranked `scale` hypotheses of R12, starting with the unit-basis check and the paper's own $L^2/D_U$ time-scale estimate.
6. Build digitized reference CSVs and an automatic error script for Figs. 6/7.
7. Complete the Mohanty mesh, time-step, and flux-decomposition verification.
8. After the Mohanty diffusion couple is reproduced, build the KKS multiphase skeleton in `coord_type = RZ` with `phase_consistent` switching, and pass the phase-diagram check (R6, R7, R11).
9. Add Wen’s $\xi$ only after the irradiation-free KKS model passes; run both `xi_target` branches (R5).
10. Upgrade temperature from a prescribed field to Jung’s quasi-static heat-conduction variable only after the Wen stage passes (R8).
11. Turn on the Jung coupling in the order “constant thermal conductivity $\rightarrow k(T,c)$ $\rightarrow$ porosity $\rightarrow$ sodium infiltration,” resolving the Eq. (31) reading first (R9).
12. Enter two dimensions only after all one-dimensional verification is complete.

## 12. Main Reference Files

- Mohanty et al. (2011): `palm/doc/Mohanty 等 - 2011 - Thermotransport in γ(bcc) U–Zr alloys A phase-field model study.pdf`
- Wen et al. (2022): `palm/doc/Wen 等 - 2022 - A phase-field model with irradiation-enhanced diffusion for constituent redistribution in U-10wt%Zr.pdf`
- Jung et al. (2025): `palm/doc/Jung 等 - 2025 - Investigating Constituent Redistribution in U-Zr Metallic Fuels A Phase-Field Approach Incorporatin.pdf`
- Current first-step notes: `palm/problem/uzr/s1_mohanty/ReadMe.md`
- Current first-step input: `palm/problem/uzr/s1_mohanty/s1_mohanty.i`
- MOOSE Soret Kernel: `moose/modules/phase_field/src/kernels/SoretDiffusion.C`
- MOOSE KKS examples: `moose/modules/phase_field/test/tests/KKS_system/`
- MOOSE heat-conduction examples: `moose/modules/heat_transfer/test/tests/`

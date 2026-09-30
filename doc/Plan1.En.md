# U–Zr Fuel Constituent Redistribution under a Temperature Gradient: A Stepwise MOOSE Reproduction and Extension Plan

> Subject: constituent redistribution of U–Zr metallic fuel under a temperature gradient  
> Numerical framework: MOOSE, phase-field module, split Cahn–Hilliard equation  
> Scope: strictly limited to the physical models explicitly used by Mohanty et al. (2011), Wen et al. (2022), and Jung et al. (2025)  
> Recommended implementation order: one-dimensional item-by-item verification $\rightarrow$ one-dimensional fully coupled multiphysics $\rightarrow$ two-dimensional geometric extension

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

The paper uses $E_v^f=1.20$ eV and cites $C_v^r=7.04\times10^{-7}$ at about 1000 K and 30 kW/m. The paper does not explicitly solve defect-concentration evolution in the phase-field calculation, so the reproduction should treat $C_v^r$ as a prescribed value or a parameter function, not as a new defect PDE.

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

The phase-region order in the paper’s main figures should decide which of the two is used. They must not be swapped in the code without a record.

Wen 2022 uses a prescribed radial temperature profile and does not solve fully coupled heat conduction with composition-dependent thermal conductivity. The irradiation-enhancement stage should therefore keep the temperature field prescribed and avoid introducing the thermal feedback of Jung 2025 at the same time.

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

The porosity correction factor is implemented according to paper equation (31). The treatment of porosity with burnup is:

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

### 3.1 The Current Model Contains Two Artificial Corrections

The current input file uses:

- An overall atomic-mobility scale factor `scale = 0.035`;
- An added composition weight on Mohanty equation (3).

The revision notes in the current README still retain `scale = 0.05`, while the parameter table and the actual input file use `0.035`. Later cleanup should take the input file as authoritative and correct the old value in the documentation. The current `kappa_c` comment is written as “thermal conductivity”; it should be the compositional gradient-energy coefficient.

These two treatments improve the fit to Figs. 6/7, but they are not the original Mohanty 2011 equations. Two parallel cases should be kept:

- `paper_equation.i`: strictly use paper equations (2)–(4), with no composition weight and no prior fitting of the time scale;
- `calibrated.i`: keep the current stabilized and fitted version.

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
- `KKSGlobalFreeEnergy`
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
7. Check how the interface width varies with $\kappa_\phi/\omega$.

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
5. Still compute the thermal mobility from paper equations (29)–(31), and do not multiply it by $\xi$ unless the paper equations and text are checked again;
6. Compare irradiated and unirradiated results.

### Suggested Material Properties

- `thermal_vacancy_concentration`
- `irradiation_vacancy_concentration`
- `irradiation_enhancement_factor`
- `chemical_mobility_unirradiated`
- `chemical_mobility_irradiated`

These relations may first be implemented with `ADParsedMaterial` or an AD parsed material. A C++ Material does not need to be written immediately. Add `ADU_ZrIrradiationMobilityMaterial` only when piecewise functions, numerical safeguards, or performance become a problem.

### Paper Reproduction Settings

- One-dimensional length of 2170 $\mathrm{\mu m}$;
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
- One-parameter sensitivity analyses are performed separately for $\xi$, $\kappa_\gamma$, and $\omega$;
- Parameters given in the literature, cited parameters, and trial-fit parameters are clearly distinguished.

## Stage 5: Jung 2025 Heat-Conduction Baseline

### Objective

First build a fully coupled $c$–$\phi$–$T$ baseline without porosity or sodium infiltration, then add thermal-conductivity corrections.

### Heat-Equation Objects

Suggested objects:

- `ADHeatConductionTimeDerivative`
- `ADHeatConduction`
- `ADBodyForce` or an equivalent AD heat-source Kernel
- A surface `ADDirichletBC`, or a Dirichlet BC consistent with the paper

Material properties should include at least:

- `density`
- `specific_heat`
- `thermal_conductivity`
- `volumetric_heat_source`

If the paper does not give $\rho$ and $c_p$, and the final result is essentially a quasi-steady thermal field, first verify the steady heat equation and mark the transient thermal parameters as missing. Do not assign $\rho c_p$ without a source merely so that the transient heat equation can be run.

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

1. Verify the steady cylinder/slab analytical solution for constant $k$;
2. Set only $k=k(T)$;
3. Set $k=k(T,c)$;
4. Turn on one-way $c$–$\phi$–$T$ coupling: temperature affects the phase field, but the phase field does not feed back into thermal conductivity;
5. Turn on two-way coupling: $c$ affects $k$, and $k$ affects $T$;
6. Check the difference between a monolithic solve and a split solve.

The paper does not report the total simulation time, yet it lets burnup increase linearly over the total simulation time and uses Hirschhorn’s optimized mobility to compress a real irradiation process of about 284 EFPD. Jung results can therefore be used to reproduce end-state profiles, but they cannot be used to verify real-time kinetics until the authors’ input files or supplementary data are obtained.

### Geometric Notes

The paper treats the left end of the one-dimensional interval as the fuel center and the right end as the surface, but a one-dimensional finite-element implementation may use a slab operator. Two clearly distinguished inputs should be created:

- `paper_1d.i`: reproduce the paper’s actual one-dimensional implementation as closely as possible;
- `cylindrical_check.i`: use axisymmetric radial heat conduction or an equivalent radial weight and check the effect of the geometric term.

Results with the geometric correction cannot replace the paper-reproduction results. They should be a separate sensitivity study.

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
2. Independently verify the limit of $P_c$ when $P_{\mathrm{Na}}=0$;
3. Check the value of $k_{\mathrm{Na}}/k_0$ over the full temperature and composition ranges;
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

## Stage 9: Two-Dimensional Extension

### Objective

Without adding physics beyond the papers, check the effects of the one-dimensional assumption, geometric discretization, and nonuniform boundaries on the results.

### Order

1. Reproduce the one-dimensional solution on a two-dimensional rectangular domain, for code regression;
2. A two-dimensional axisymmetric $r$–$z$ domain with uniform axial conditions;
3. Allow axial power or surface temperature to vary;
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

The Jung model should first use a monolithic solve, because temperature, phase fractions, and composition form clear feedback. If the monolithic system is hard to converge, then try:

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

The integration measure for the one-dimensional slab, the cylinder, and the two-dimensional model must be consistent with the geometry.

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

- As printed, equation (3) contains no composition weight, which may produce limiting behavior that differs from intuition or from the numerical results;
- The paper states that $\kappa_c=0$ in the single phase; the nonzero gradient energy in the current case is extra numerical regularization;
- The atomic-mobility prefactors are not fully tabulated in the main text and must be read from figures or traced through citations;
- The free energy comes from a commercial database and cannot be uniquely recovered from the paper text alone;
- The text says 1405 K, while the figure axis may show 1400 K;
- $c_x=0,c_{xxx}=0$ in paper equation (6) is not strictly equivalent to a total no-flux boundary that includes the thermal-migration term;
- The time-scale correction cannot be uniquely determined from the final curves alone.

### Wen 2022

- $C_v^r$ comes from an external rate-theory code and is not a derivation that can be repeated inside the paper;
- The irradiation-enhancement factor lacks a complete dependence on irradiation conditions;
- The $\alpha/\beta$ phase-fraction and mobility interpolation weights are opposite to the paper’s phase-region text;
- $L_0$ is not given, and $\kappa_c$ is given only as an order of magnitude;
- Some free-energy expressions use a fixed temperature or have unclear typesetting;
- The abstract says the effect of irradiation on thermal diffusion is discussed, but the governing equation only sets $M_c^I=\xi M_c$;
- The time discretizations of explicit finite differences and MOOSE implicit finite elements differ;
- The phase-interface parameters involve a clear trial-and-fit process;
- Zr enrichment near the fuel surface cannot be explained by this model.

### Jung 2025

- Hirschhorn’s optimized mobility is used to compress a real time of 284 EFPD, so the strict time scale cannot be interpreted directly;
- The total simulation time is not reported, so the slope of burnup versus time cannot be uniquely recovered;
- Transient thermal parameters such as $\rho$ and $c_p$ are incomplete in the parameter table of the main text;
- $\beta_0$, $\widetilde Q_U^*$, $\widetilde Q_{\mathrm{Zr}}^*$, and the optimized mobility scale are incomplete;
- The thermal-conductivity coefficient numbering and typesetting may be misaligned;
- The formula text for converting atom fraction to weight fraction is ambiguous;
- The $\alpha/\beta$ switching direction, the double-well function, and the definition of the sodium-infiltration amount are all ambiguous;
- The paper does not report a compositional gradient energy, so the strict branch should take $\kappa_c=0$;
- The porosity model is an algebraic average over phase regions, not an evolution of pore structure;
- Both LHGR and burnup use a simplified linear history;
- The paper acknowledges deviations in the surface $\alpha$ phase and the surface concentration;
- Whether the paper’s “one-dimensional” treatment includes the strict cylindrical geometric term needs further confirmation from the input file or the authors’ code.

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
- Interface width is mesh-independent;
- The unirradiated baseline is stable.

### M3: Wen Irradiation Enhancement Is Reproduced

- The difference between unirradiated and irradiated cases is reproduced;
- Three phase regions and the Zr well appear;
- Sensitivity to $\xi$ and the interface parameters is completed;
- No defect equations outside the paper are introduced.

### M4: The Jung Thermal-Feedback Baseline Is Credible

- Analytical heat-conduction benchmarks pass;
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

1. Freeze the current `s1_mohanty.i` as the calibrated version.
2. Build the original-equation version of the Mohanty paper, removing the mobility scale and the extra composition weight.
3. Implement and test `ADThermalMigration` in PALM.
4. Build digitized reference CSVs and an automatic error script for Figs. 6/7.
5. Complete the Mohanty mesh, time-step, and flux-decomposition verification.
6. After the Mohanty diffusion couple is reproduced, build the KKS multiphase skeleton.
7. Add Wen’s $\xi$ only after the irradiation-free KKS model passes.
8. Upgrade temperature from a prescribed field to Jung’s heat-conduction variable only after the Wen stage passes.
9. Turn on the Jung coupling in the order “constant thermal conductivity $\rightarrow k(T,c)$ $\rightarrow$ porosity $\rightarrow$ sodium infiltration.”
10. Enter two dimensions only after all one-dimensional verification is complete.

## 12. Main Reference Files

- Mohanty et al. (2011): `palm/doc/Mohanty 等 - 2011 - Thermotransport in γ(bcc) U–Zr alloys A phase-field model study.pdf`
- Wen et al. (2022): `palm/doc/Wen 等 - 2022 - A phase-field model with irradiation-enhanced diffusion for constituent redistribution in U-10wt%Zr.pdf`
- Jung et al. (2025): `palm/doc/Jung 等 - 2025 - Investigating Constituent Redistribution in U-Zr Metallic Fuels A Phase-Field Approach Incorporatin.pdf`
- Current first-step notes: `palm/problem/uzr/s1_mohanty/ReadMe.md`
- Current first-step input: `palm/problem/uzr/s1_mohanty/s1_mohanty.i`
- MOOSE Soret Kernel: `moose/modules/phase_field/src/kernels/SoretDiffusion.C`
- MOOSE KKS examples: `moose/modules/phase_field/test/tests/KKS_system/`
- MOOSE heat-conduction examples: `moose/modules/heat_transfer/test/tests/`

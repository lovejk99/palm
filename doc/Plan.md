# U–Zr Constituent Redistribution under a Temperature Gradient: Stepwise MOOSE Reproduction Plan

> Code: PALM (MOOSE application; phase-field and heat-transfer modules), split Cahn–Hilliard form  
> Part 1: Mohanty et al. (2011), single-phase γ thermotransport  
> Part 2: Wen et al. (2022), multiphase KKS model with irradiation-enhanced diffusion  
> Part 3: Jung et al. (2025), heat conduction, thermal conductivity, porosity and sodium infiltration  
> Supporting: Hirschhorn et al. (2019), the model Jung builds on (source of Jung's parameters); Li et al. (2010), the γ atomic mobilities used by Hirschhorn and Jung; Hirschhorn et al. (2020), a γ-mobility cross-check  
> Starting point: `problem/uzr/s1_mohanty` (built-in `SoretDiffusion`, documented in its `ReadMe.md`) and the PALM kernel `ADThermalMigration` (5 tests passing)

---

## 1. Strategy

### 1.1 Roadmap

Each row is one small, runnable step. Every step changes exactly one thing relative to the step above it.

| Step | Physics added (or change made) | Key data used for validation | Fitted |
| --- | --- | --- | --- |
| 1.1 | Replace `SoretDiffusion` by `ADThermalMigration` (framework only) | current `s1_mohanty` result | — |
| 1.2 | Single-phase chemical diffusion (isothermal diffusion couple) | Mohanty Fig. 8a, dashed curve | $s_M$ |
| 1.3 | + thermotransport in a linear temperature field | Mohanty Figs. 6 and 7 | — |
| 1.4 | + pre-existing composition gradient | Mohanty Fig. 8a/8b, solid curves | — |
| 1.5 | Closure: convergence, freeze the framework | — | — |
| 2.1 | Two-phase thermodynamics $f_{\alpha\beta}$, $f_\gamma$ (no PDE) | Wen Fig. 2, γ curve | — |
| 2.2 | + KKS multiphase diffusion and Allen–Cahn, prescribed $T(r)$ | Wen Fig. 2, γ curve | — |
| 2.3 | + multiphase thermotransport | Wen Fig. 2, Zr curve | $s_W$ |
| 2.4 | + irradiation-enhanced diffusion $\xi(T)$ | Wen Fig. 3, DP-81 PIE | — |
| 2.5 | Interface width versus $\omega$; closure | Wen Table 2, Fig. 5B | — |
| 3.1 | Heat conduction with a volumetric source, prescribed $k(r)$ | Jung Fig. 5b | — |
| 3.2 | Switch to Jung's thermodynamics and kinetics under a prescribed $T(r)$ | Jung Fig. 3b | $t_{\mathrm{end}}$ |
| 3.3 | + two-way coupling through $k_0(T,c)$ and $\dot q(B)$ | Jung Fig. 5, case $P=0$ | — |
| 3.4 | + porosity from phase fractions | Jung Fig. 5 | — |
| 3.5 | + sodium infiltration | Jung Fig. 6, Fig. 3b, DP-81 PIE | — |
| 3.6 | Second pin DP-11 (no new physics); closure | Jung Fig. 4, DP-11 PIE | — |

### 1.2 Rules for every step

1. One change per step. A step's input file differs from the previous step's input only by the block that implements the new process. Steps 1.1 and 3.2 are the two exceptions: they change the framework or the parameter set, not the physics.
2. Validation by key data. Each step is judged against one to four numbers or features read from one paper figure, each with a stated tolerance, plus at most two sanity checks (a limiting case or a conservation check). Mesh and time-step convergence is done once per part, in its closing step.
3. Printed equations come first. A printed equation is replaced only when it contradicts the paper's own figures or is mathematically impossible. Each such replacement is one row in the part's "Paper issues and decisions" table, with the evidence.
4. One frozen fit per part. Each part has at most one fitted quantity. It is fitted once, on the simplest case of that part, and is never refitted. All later steps of that part are predictions. Part 1 fits the mobility factor $s_M$, Part 2 the γ-phase mobility factor $s_W$ (Wen does not report the γ prefactor), and Part 3 the simulated time $t_{\mathrm{end}}$ (Jung does not report it).
5. No physics beyond the papers. There is no defect rate-theory PDE, no pore nucleation or growth, and no mechanics, FCCI or fission gas. Porosity acts only through the thermal conductivity, and irradiation only through $\xi$.
6. Units follow each paper (§2.4).
7. Values marked "≈" in this plan were read by eye from the PDF figures. They are good enough to plan with, but must be replaced by digitized values (§6) before any error is computed.

### 1.3 Standard checks and outputs

- Mass conservation: $\epsilon_m=\lvert\int_\Omega c\,\mathrm dV-\int_\Omega c_0\,\mathrm dV\rvert/\int_\Omega c_0\,\mathrm dV<10^{-6}$. With `coord_type = RZ`, MOOSE integral postprocessors include the $2\pi r$ weight automatically.
- Bounds: $0<c<1$, $0\le\phi\le1$ (small overshoot allowed and recorded), $T>0$, $k>0$.
- Every step directory contains one input file and a `README.md` of at most one page: what changed from the previous step, the key-data table (paper value, this run, tolerance, pass/fail), the sanity checks, and the run command.

---

## 2. Common Framework (built in Part 1, reused in Parts 2 and 3)

### 2.1 Governing equations

The conserved variable is the Zr atom fraction $c$. The split Cahn–Hilliard system with thermotransport is

$$
\frac{\partial c}{\partial t}=\nabla\cdot\left(M_c\nabla w\right)-\nabla\cdot\left(M_T\frac{\nabla T}{T}\right),
\qquad
w=\frac{\delta F}{\delta c}.
$$

The Zr flux is $\mathbf J=-M_c\nabla w+M_T\nabla T/T$, so $M_T>0$ moves Zr up the temperature gradient, toward the hot end or the fuel center, as in all three papers. Part 2 adds the γ order parameter $\phi$ (Allen–Cahn) and the KKS phase concentrations $c_{\alpha\beta}$, $c_\gamma$. Part 3 turns $T$ into a nonlinear variable.

### 2.2 MOOSE objects

Variable `w` carries the mass balance and variable `c` carries the definition of `w`, as in the standard MOOSE split form.

| Term | Part 1 | Parts 2 and 3 |
| --- | --- | --- |
| $\partial c/\partial t$ | `ADCoupledTimeDerivative` (`variable = w`, `v = c`) | `ADCoupledTimeDerivative` |
| $\nabla\cdot(M_c\nabla w)$ | `ADSplitCHWRes` | `ADSplitCHWRes` |
| $-\nabla\cdot(M_T\nabla T/T)$ | `ADThermalMigration` | `ADThermalMigration` |
| $w=\delta F/\delta c$ | `ADSplitCHParsed` with $\kappa_c=0$ | `KKSSplitCHCRes`, `KKSPhaseConcentration`, `KKSPhaseChemicalPotential` |
| $\phi$ equation | — | `TimeDerivative`, `KKSACBulkF`, `KKSACBulkC`, `ACInterface` |
| Temperature | AuxVariable with `FunctionAux` | Part 2: AuxVariable; Part 3: `ADHeatConduction`, `ADBodyForce`, `DirichletBC` at $r_0$ |

Material property names stay the same in all parts: `M_c` and `M_T` (AD properties), the free energies (`f_chem` in Part 1; `f_ab` and `f_g` in Parts 2 and 3), and the temperature variable `T`. A later step therefore replaces material blocks and leaves the kernel block untouched. In Parts 2 and 3 the KKS kernels are not AD, so the free energies, the switching function $h(\phi)$ and the Allen–Cahn mobility $L$ are `DerivativeParsedMaterial` objects. `M_c` and `M_T` are `ADParsedMaterial` objects of $(\phi,c_{\alpha\beta},c_\gamma,T)$ that write $h(\phi)$ and $h_{\alpha\beta}(T)$ inline, which avoids a name clash between AD and non-AD properties.

Solver: `NEWTON` with LU in all parts (1-D problems are small). Use `IterationAdaptiveDT`.

### 2.3 Why `ADThermalMigration` replaces `SoretDiffusion`

MOOSE's `SoretDiffusion` discretizes $\int_\Omega D\,Q\,c/(k_BT^2)\,\nabla T\cdot\nabla\psi\,\mathrm dV$ with $k_B$ hard-coded in eV/K. Expressing the papers' flux $M_T\nabla T/T$ with it requires the mapping $D=M_ck_BT/c$ and $Q=-M_T/M_c$, which the current `s1_mohanty.i` uses. That mapping has three problems, and each gets worse in Parts 2 and 3:

1. Incomplete Jacobian. The kernel differentiates only its explicit dependence on $c$ and $T$, not the derivatives of $D(c,T)$ and $Q(c,T)$. That is why the current case needs `PJFNK`. In Parts 2 and 3, $M_T$ depends on $\phi$, $c_{\alpha\beta}$, $c_\gamma$ and $T$, and those Jacobian blocks would be missing from the preconditioner entirely.
2. Fixed units. A hard-coded $k_B$ in eV/K forces eV-based units, while Parts 2 and 3 use J/mol and J/m³.
3. Artificial singularity. Splitting $M_T$ into $D\,Q\,c$ puts $1/c$ into $D$, which is singular as $c\to0$. That limit is actually reached: the α-U phase of Parts 2 and 3 has $c\approx0.01$.

`ADThermalMigration` (in `include/kernels/` and `src/kernels/`) discretizes

$$
R_{\mathrm{TM}}=-\int_\Omega\frac{M_T}{T}\,\nabla T\cdot\nabla\psi\,\mathrm dV.
$$

It reads $M_T$ directly as an AD material property, contains no physical constants, and obtains the full Jacobian with respect to every variable $M_T$ depends on through automatic differentiation. Its tests in `test/tests/kernels/ad_thermal_migration/` cover an exact steady state, a reversed gradient, a `PetscJacobianTester` run with $M_T(c,T)$, and equivalence with `SoretDiffusion`. MOOSE's own `ADCHSoretMobility` has the same weak form with the opposite sign (its `mobility` equals $-M_T$). The PALM kernel is kept because its sign matches the papers' $M_T$ directly.

### 2.4 Units and parameter conversion

| Quantity | Part 1 (Mohanty) | Parts 2 and 3 (Wen, Jung) | Conversion |
| --- | --- | --- | --- |
| Length | µm | m | — |
| Free energy | eV/atom | J/m³, i.e. $f_{\mathrm{mol}}/V_m$ | 1 eV/atom = 96485.33 J/mol; $V_m=1.249\times10^{-5}$ m³/mol |
| $\beta_i$ | µm²/(eV s), per atom | m² mol/(J s), per mole | $\beta^{\mathrm{mol}}=\beta^{\mathrm{atom}}/N_A$; the prefactors below are used in Parts 1 and 2, while Part 3 uses Li et al. (2010) (§5.1) |
| $\beta_{0,U}$ | $1.7204\times10^{6}$ | $1.783\times10^{-11}$ | via $1.074\times10^{13}$ m²/(J s) per atom |
| $\beta_{0,\mathrm{Zr}}$ | $6.9045\times10^{7}$ | $7.156\times10^{-10}$ | via $4.31\times10^{14}$ m²/(J s) per atom |
| $\widetilde Q_U^*$ | +0.15604 eV (= +15.1 kJ/mol) | Wen (Hofman): +110 kJ/mol; Jung (Hirschhorn 2019, Mohanty's averages): +15.1 kJ/mol | 1 eV/atom = 96.485 kJ/mol |
| $\widetilde Q_{\mathrm{Zr}}^*$ | −1.09226 eV (= −105.4 kJ/mol) | Wen (Hofman): +105 kJ/mol; Jung: +105.4 kJ/mol (magnitude) | — |
| $M_c$ | µm²/(eV s) | m⁵/(J s) | — |
| $M_T$ | µm²/s | m²/s | — |

Sign convention: Mohanty writes $M_T\propto\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*$ with a negative $\widetilde Q_{\mathrm{Zr}}^*$. Wen and Jung write $\widetilde Q_U^*\beta_U+\widetilde Q_{\mathrm{Zr}}^*\beta_{\mathrm{Zr}}$ with positive magnitudes. The Zr term has the same sign in both. A signed Mohanty value must never be copied into the Wen/Jung form, because that reverses the Zr flux.

### 2.5 Directory layout

```text
palm/problem/uzr/
├── s1_mohanty/                 existing SoretDiffusion case, kept unchanged as a reference
├── data/                       digitized paper curves, one CSV per curve (§6)
├── scripts/                    plot_exodus_cx_profile.py, compare.py (.e + CSV -> key-data table)
├── p1_mohanty/
│   ├── step1_1_kernel_swap/
│   ├── step1_2_isothermal_couple/
│   ├── step1_3_thermotransport/
│   ├── step1_4_couple_in_gradient/
│   └── README.md               filled key-data tables of Part 1 and the value of s_M
├── p2_wen/
│   ├── step2_1_thermo/         Python only: common tangent and lever rule
│   ├── step2_2_kks_diffusion/
│   ├── step2_3_thermotransport/
│   ├── step2_4_irradiation/
│   ├── step2_5_interface/
│   └── README.md
└── p3_jung/
    ├── step3_1_heat/
    ├── step3_2_prescribed_T/
    ├── step3_3_k_coupling/
    ├── step3_4_porosity/
    ├── step3_5_sodium/
    ├── step3_6_dp11/
    └── README.md
```

No C++ beyond `ADThermalMigration` is planned; everything else is parsed materials and functions. A C++ material is written only if a parsed expression becomes unreadable (the most likely candidate is Jung's $k$). MOOSE sources are not modified.

---

## 3. Part 1: Mohanty et al. (2011), Single-Phase γ U–Zr

Goal: build and validate the framework on the simplest system, a single bcc phase with a prescribed linear temperature field, adding chemical diffusion, then thermotransport, then their competition.

### 3.1 Model and setup

$$
M_c=c(1-c)\left[c\,\beta_U+(1-c)\,\beta_{\mathrm{Zr}}\right],
\qquad
M_T=c(1-c)\left[\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right],
\qquad
\beta_i=s_M\,\beta_{0,i}\exp\left(-\frac{H_i}{RT}\right)
$$

$$
f_{\mathrm{chem}}=(1-c)\,G_U^{\mathrm{bcc}}(T)+c\,G_{\mathrm{Zr}}^{\mathrm{bcc}}(T)+RT\left[c\ln c+(1-c)\ln(1-c)\right]+c(1-c)\sum_{k=0}^{4}L_k(T)\,(1-2c)^k
$$

In eV/atom, $w=\partial f_{\mathrm{chem}}/\partial c$ and $\kappa_c=0$. The paper's $V_m$ factors cancel when energies and mobilities are both per atom.

| Item | Value | Source |
| --- | --- | --- |
| Domain | 1-D slab, 300 µm, 300 linear elements | Mohanty §3 |
| Temperature | linear, $T(0)=1400$ K, $T(300\ \mu\mathrm m)=1050$ K | Fig. 6 top axis (text: 1405 K) |
| Initial condition | $c_0=0.39$; diffusion couple 0.39 / 0.225 with a tanh step of half-width 2 µm | Figs. 6, 8 |
| Boundary condition | natural BC of the split form (zero total flux) | — |
| Time | 30 days = $2.592\times10^{6}$ s; output at 5, 10, 30 days | Figs. 6–8 |
| $\beta_{0,U}$, $\beta_{0,\mathrm{Zr}}$ | $1.7204\times10^{6}$, $6.9045\times10^{7}$ µm²/(eV s) | derived from Figs. 2 and 4, see `s1_mohanty/ReadMe.md` §4.3 |
| $H_U$, $H_{\mathrm{Zr}}$ | 128, 195 kJ/mol | Mohanty §2 (Hofman 1996) |
| $\widetilde Q_U^*$, $\widetilde Q_{\mathrm{Zr}}^*$ | $+2.5\times10^{-20}$ J, $-17.5\times10^{-20}$ J | Mohanty §2 |
| $L_0\ldots L_4$ | as in `s1_mohanty.i` | Chevalier et al. 2004 (Mohanty ref. [21]) |
| $G_U^{\mathrm{bcc}}$, $G_{\mathrm{Zr}}^{\mathrm{bcc}}$ | SGTE unary data | Dinsdale 1991; the bcc-U expression changes at 1049 K, so the whole 1050–1400 K domain uses the high-temperature range |

### 3.2 Paper issues and decisions

| Item | In the paper | Decision |
| --- | --- | --- |
| Hot-end temperature | text 1405 K, figure axes 1400 K | 1400 K |
| $\beta_{0,i}$ | not tabulated, only plotted (Fig. 2) | values derived in `s1_mohanty/ReadMe.md` |
| Time scale | with the given $\beta_i$ the interdiffusion time $L^2/\widetilde D$ is ≈ 4 days at 1225 K, so the profiles would freeze within days, yet Figs. 6–8 still evolve at 30 days | one factor $s_M$ multiplying both $\beta_i$, fitted in Step 1.2 |
| Eq. (3) | no composition weights inside the bracket; this is the correct Darken/Jaffe–Shewmon form | used as printed (the weights in the current `s1_mohanty.i` are dropped from Step 1.2 on) |
| Free energy | full bcc energy from the ThermoCalc SGTE database | Redlich–Kister terms plus SGTE end-members. In a temperature gradient the end-members add a $\partial^2f/\partial c\,\partial T$ term to $\nabla w$; its size is measured in Step 1.3 |
| $\kappa_c$ | "vanishes in a single-phase alloy" | 0 |
| Eq. (8) | Dirichlet temperatures at both ends and zero heat flux | only the two end temperatures (linear $T$, AuxVariable) |
| Hot-end value at 5 days | Fig. 6 ≈ 0.455, Fig. 7 ≈ 0.44 | the tolerance in Step 1.3 covers both |

### 3.3 Steps

#### Step 1.1: Kernel swap (framework, no new physics)

- Change: in the current `s1_mohanty.i` model, with its physics unchanged (scale 0.035, weighted $M_T$, no end-members, $\kappa_c=10^{-5}$), replace `SoretDiffusion` and the `Qeff`/`D_soret` materials by `ADThermalMigration` reading $M_T$ directly. Convert $f$ and $M_c$ to AD materials (`ADDerivativeParsedMaterial`) and the kernels to `ADSplitCHParsed` and `ADSplitCHWRes`, and switch `solve_type` to `NEWTON`.
- Key data: the existing result `s1_mohanty_out.e`.

| $t$ (days) | $c(0)$ | $c(300\ \mu\mathrm m)$ |
| --- | --- | --- |
| 5 | 0.457 | 0.355 |
| 10 | 0.470 | 0.336 |
| 30 | 0.492 | 0.266 |

- Pass: profiles identical to $10^{-5}$ (absolute) at all output times; Newton converges in fewer nonlinear iterations than the current `PJFNK` run.
- Sanity: `./run_tests --re ad_thermal_migration` passes.
- After this step `s1_mohanty` is a reference only.

#### Step 1.2: Single-phase chemical diffusion (isothermal diffusion couple, 1225 K)

- Adds: the printed-equation model with chemical diffusion only: $f_{\mathrm{chem}}$ with end-members, $\kappa_c=0$, $M_c$ from Eq. (2), no thermotransport kernel. Uniform $T=1225$ K, diffusion couple 0.39 (left) / 0.225 (right) on $[-150,150]$ µm.
- Fit: $s_M$, once. The calibrated `s1_mohanty` used 0.035, which is a reasonable first guess. $s_M$ describes Mohanty's simulation, not the alloy, so it is not carried into Parts 2 and 3: Hirschhorn et al. (2020) fitted a γ mobility at 1023 K that is 7–14 times the $s_M=1$ value (Step 2.3).
- Key data (Fig. 8a dashed curve, 30 days): $c(-150)\approx0.335$, $c(0)\approx0.31$, $c(+150)\approx0.28$.
- Pass: all three within ±0.01 with a single $s_M$. Record $s_M$.
- Sanity: a uniform $c_0=0.39$ stays uniform at 1225 K; $\epsilon_m<10^{-6}$. Halving the initial tanh width changes the 30-day profile by less than 0.002.

#### Step 1.3: Thermotransport in the homogeneous alloy

- Adds: the linear $T(x)$ and `ADThermalMigration` with $M_T$ from Eq. (3) as printed. $s_M$ stays frozen. Output $J_c=-M_c\,\partial_xw$ and $J_T=M_T\,\partial_xT/T$ as AuxVariables (`VariableGradientComponent` plus `ParsedAux`).
- Key data (Figs. 6, 7; Fig. 5c is the same 30-day case):

| $t$ (days) | $c(0)$, hot end | $c(300\ \mu\mathrm m)$, cold end | Shape |
| --- | --- | --- | --- |
| 5 | ≈ 0.455 (Fig. 7: ≈ 0.44) | ≈ 0.35 | minimum ≈ 0.365 near $x\approx190$ µm, local maximum ≈ 0.375 near $x\approx270$ µm |
| 10 | ≈ 0.47 (Fig. 7: ≈ 0.46) | ≈ 0.34 | same features, weaker |
| 30 | ≈ 0.50 | ≈ 0.28 | shoulder ≈ 0.33 near $x\approx250$ µm |

- Pass: $c(0)$ within ±0.015 at 5, 10 and 30 days; $c(300)$ at 30 days within ±0.03; the interior maximum near the cold end is present at 5 days. Also plot $c(0,t)$ against Fig. 7.
- Sanity: (a) reversing the gradient gives the mirrored profile $c(x)\to c(300-x)$. (b) With $\widetilde Q_U^*=\widetilde Q_{\mathrm{Zr}}^*=0$, only the $\partial^2f/\partial c\,\partial T$ part of $\nabla w$ remains; record $c(0,30\ \mathrm d)-0.39$ as the size of that contribution.

#### Step 1.4: Combined composition and temperature gradients

- Adds: the diffusion-couple initial condition of Step 1.2 placed in the temperature field of Step 1.3, in two orientations. No new parameters.
- Key data (Fig. 8, 30 days):

| Case | Orientation | Key points |
| --- | --- | --- |
| 8a | hot end (1400 K) at $x=-150$ µm, the Zr-rich side | $c(-150)\approx0.44$, $c(0)\approx0.31$, $c(+150)\approx0.16$: Zr moves against its own gradient |
| 8b | hot end at $x=+150$ µm, the Zr-poor side | $c(+150)\approx0.38$; minimum ≈ 0.27 near $x\approx-40$ µm; local maximum ≈ 0.31 near $x\approx-135$ µm; $c(-150)\approx0.28$ |

- Pass: the listed values within ±0.02, extremum positions within ±20 µm.
- Flux check: in 8a, $J_T$ and $J_c$ have opposite signs over most of the domain; in 8b they have the same sign near $x=0$. Mohanty's statement that $J_T$ is 4–5 orders larger than $J_c$ is not used as a criterion, because the two fluxes must approach balance as the profile saturates. Report $\lvert J_T/J_c\rvert$ at $x=0$ versus time instead.

#### Step 1.5: Part 1 closure

- Convergence on Step 1.3 at 30 days: meshes of 150, 300 and 600 elements, and `dtmax` × 1, ½, ¼. $c(0)$ and $c(300)$ must change by less than 0.002.
- Freeze the framework: the kernel block, the material names `M_c`, `M_T` and `f_chem`, the postprocessors, and the scripts. Write `p1_mohanty/README.md` with all filled key-data tables and $s_M$.
- What Part 2 inherits: the `w`-equation block (time derivative, `ADSplitCHWRes`, `ADThermalMigration`) unchanged, and the γ-phase forms of $M_c$ and $M_T$ (Wen's Eqs. 27 and 31 have the same structure as Mohanty's Eqs. 2 and 3).

---

## 4. Part 2: Wen et al. (2022), Multiphase KKS Model with Irradiation-Enhanced Diffusion

Goal: go from one phase to three phases (α, β, γ) with one γ order parameter, in a cylindrical fuel slug (DP-81) with a prescribed temperature profile. Then add irradiation-enhanced diffusion, the paper's main claim.

### 4.1 Model and setup

Phase fractions are $\eta_\gamma=\phi$, $\eta_\beta=h_{\alpha\beta}(1-\phi)$ and $\eta_\alpha=(1-h_{\alpha\beta})(1-\phi)$, with the physical α/β switch (§4.2):

$$
h_{\alpha\beta}(T)=\frac12\left[1+\tanh\frac{T-935}{2}\right],
\qquad
h(\phi)=\phi^3(6\phi^2-15\phi+10),
\qquad
g(\phi)=\phi^2(1-\phi)^2 .
$$

Mobilities (Wen Eqs. 24–31):

$$
M_c=\xi\left[(1-h_{\alpha\beta})(1-h)M_c^{\alpha}+h_{\alpha\beta}(1-h)M_c^{\beta}+h\,M_c^{\gamma}\right],
\qquad
M_T=(1-h_{\alpha\beta})(1-h)M_T^{\alpha}+h_{\alpha\beta}(1-h)M_T^{\beta}+h\,M_T^{\gamma}
$$

$$
M_c^{i}=\frac{D_0^iV_m}{RT}e^{-Q_i/RT},
\qquad
M_T^{i}=\frac{D_0^i\widetilde Q_i^*}{RT}e^{-Q_i/RT}\quad(i=\alpha,\beta)
$$

$$
M_c^{\gamma}=V_m\,c_\gamma(1-c_\gamma)\left[c_\gamma\beta_U+(1-c_\gamma)\beta_{\mathrm{Zr}}\right],
\qquad
M_T^{\gamma}=c_\gamma(1-c_\gamma)\left[\widetilde Q_U^*\beta_U+\widetilde Q_{\mathrm{Zr}}^*\beta_{\mathrm{Zr}}\right],
\qquad
\beta_i=s_W\,\beta_{0,i}e^{-H_i/RT}
$$

Irradiation enhancement (Wen Eqs. 21, 23):

$$
\xi(T)=1+\frac{C_v^r}{C_v^e(T)},
\qquad
C_v^e=\exp\left(-\frac{E_v^f}{k_BT}\right),
\qquad
E_v^f=1.20\ \mathrm{eV},
\quad
C_v^r=7.04\times10^{-7}.
$$

Allen–Cahn mobility: $L=L_0\exp(-Q_\gamma/RT)$ with $Q_\gamma=128000-107000c_\gamma+174000c_\gamma^2$ J/mol.

| Item | Value | Source |
| --- | --- | --- |
| Geometry | 1-D radial, `coord_type = RZ`, $0\le r\le2170$ µm | Wen: "axisymmetric" domain |
| Mesh | 1000 elements (≈ 2 µm); the interface is ≈ 40 µm wide and the $h_{\alpha\beta}$ band ≈ 60 µm | — |
| Temperature | AuxVariable from the digitized orange curve of Wen Fig. 2 (988 K center, 900 K surface; Yacout 1996 profile), as `PiecewiseLinear` | Fig. 2 |
| Initial condition | $c=0.225$; $\phi$, $c_{\alpha\beta}$, $c_\gamma$ from the lever rule at the local $T(r)$ (Step 2.1) | — |
| Boundary condition | zero flux for $c$ and $\phi$ at $r=0$ and $r=r_0$ (natural BC) | Wen |
| Time | $1.3\times10^{5}$ s; Wen used explicit finite differences with 1 s steps, here implicit adaptive steps | Wen |
| Units | m, s, J/m³; $V_m=1.249\times10^{-5}$ m³/mol | Wen Eq. 9 |
| $f_{\alpha\beta}$, $f_\gamma$ | Wen Eqs. 10–11 divided by $V_m$, with the corrections in §4.2 | Wen |
| SGTE end-members | Jung 2025 Appendix A, with $T\ln T$ | §5.2 |
| $D_0^\alpha$, $D_0^\beta$ | $2.0\times10^{-7}$, $5.7\times10^{-5}$ m²/s | Wen Table 1 |
| $Q_\alpha$, $Q_\beta$ | 170, 180 kJ/mol | Wen Table 1 |
| $\widetilde Q_\alpha^*$, $\widetilde Q_\beta^*$ | 200, 450 kJ/mol | Wen Table 1 |
| $\widetilde Q_U^*$, $\widetilde Q_{\mathrm{Zr}}^*$ (γ) | 110, 105 kJ/mol | Wen Table 1 |
| $H_U$, $H_{\mathrm{Zr}}$; $\beta_{0,i}$ | 128, 195 kJ/mol; Part 1 values in molar units (§2.4) | Wen Table 1; Part 1 |
| $\kappa_\gamma$, $\omega$ | $7.8\times10^{-4}$ J/m, $5.0\times10^{5}$ J/m³ | Wen Table 2, selected row |

### 4.2 Paper issues and decisions

| Item | In the paper | Evidence | Decision |
| --- | --- | --- | --- |
| α/β switch, Eqs. 7–8 | the printed formula puts α above 935 K | Wen's text, the Fig. 3 legend (β green in the middle, α blue outside), the Fig. 5 zone labels and the U–Zr phase diagram all put β above 935 K; only the Fig. 2 legend follows the printed formula | physical order: $h_{\alpha\beta}\to1$ (β) above 935 K. This is a one-sign change inside the tanh |
| U end-member in $f_{\alpha\beta}$, Eq. 10 | $G_0^{\gamma\text{-U}}$ | Wen lists $G_0^{\alpha\text{-U}}$ among the functions used, but it appears nowhere in print. With $G_0^{\alpha\text{-U}}$ the lever rule reproduces the Fig. 2 γ curve (single γ above ≈ 972 K, $\phi\approx0.45$ at 950 K, $\phi\approx0.40$ at 900 K); with $G_0^{\gamma\text{-U}}$ the values are 960 K and $\phi\approx0.49$ | $G_0^{\alpha\text{-U}}$ |
| Fixed temperatures in $f_\gamma$, Eq. 11 | 988.4 K in the mixing term, 998.4 K in $L_0$ | — | as printed (about 1 % effect on the excess term) |
| γ-phase $\beta_0$ | not reported | — | Part 1 $\beta_{0,i}$ (molar) times $s_W$, fitted in Step 2.3 |
| $L_0$ | not reported | Jung uses $L=2.5\times10^{-10}$ m³/(J s) | assumed $L_0=3.4\times10^{-4}$ m³/(J s), which gives $L=2.5\times10^{-10}$ at 950 K and $c_\gamma=0.3$ |
| $\kappa_c$ | "of order $10^{-9}$" | $10^{-5}$ of $\kappa_\gamma$; `KKSSplitCHCRes` has no gradient term | omitted |
| α/β mobilities, Eqs. 26 and 30 | no factor $c_{\alpha\beta}$ (Jung's Eqs. 15–16 have one) | — | as printed in Part 2 |
| $\xi$ | multiplies $M_c$ only (Eq. 24) | — | as printed |
| $M_c$ divided by ≈ 5 | not mentioned | Wen's free energies are Hirschhorn 2019's simplified ones, for which Hirschhorn divides $M_c$ by ≈ 5 | as printed (no factor); Part 3 applies it (§5.2) |
| $C_v^r$ | computed at 30 kW/m and ≈ 1000 K (DP-81 ran at 24 kW/m) | — | constant, as printed |
| PIE Zr enrichment at 1700–2000 µm | not explained by the model (Wen's discussion) | — | not a pass criterion |

### 4.3 Steps

#### Step 2.1: Two-phase thermodynamics (no PDE)

- Adds: $f_{\alpha\beta}$ and $f_\gamma$ as `DerivativeParsedMaterial` (J/m³), plus an independent Python implementation (`lever.py`) of the same functions with a common-tangent and lever-rule solver.
- Checks: the MOOSE and Python values of $f$, $\partial f/\partial c$ and $\partial^2f/\partial c^2$ agree at a few $(c,T)$ points (single-element MOOSE run with `MaterialRealAux` output). Then compute the lever-rule $\phi$ for $c=0.225$ at 900–990 K. Compare $\phi$ itself, not $h(\phi)$, with the γ curves.
- Key data (Wen Fig. 2 γ curve; Fig. 2 shows almost no redistribution, so its γ curve is close to local equilibrium): $\phi=1$ for $r\lesssim800$ µm ($T\gtrsim971$ K), $\phi\approx0.45$ near 950 K, $\phi\approx0.40$ at the surface (900 K).
- Expected (precomputed with the §4.2 decisions): single γ above ≈ 972 K; at 950 K $c_{\alpha\beta}\approx0.014$, $c_\gamma\approx0.54$, $\phi\approx0.45$; at 900 K $c_{\alpha\beta}\approx0.015$, $c_\gamma\approx0.68$, $\phi\approx0.40$.
- Pass: the edge of the γ zone within ±10 K, $\phi$ within ±0.05.
- Output: $\phi_0(r)$, $c_{\alpha\beta,0}(r)$, $c_{\gamma,0}(r)$ as a CSV file for the initial conditions of Step 2.2.

#### Step 2.2: KKS multiphase diffusion in the prescribed temperature field

- Adds: variables $\phi$, $c_{\alpha\beta}$, $c_\gamma$; the KKS and Allen–Cahn kernels of §2.2; $L(T,c_\gamma)$ as a `DerivativeParsedMaterial`; $M_c$ as in §4.1 with $\xi=1$; `ADThermalMigration` not active.
- Key data: Wen Fig. 2 γ curve at $1.3\times10^{5}$ s (the same numbers as Step 2.1); $c$ stays ≈ 0.225 away from zone boundaries.
- Pass: $\phi(r)$ within ±0.05 of Fig. 2; the KKS constraints $c=(1-h)c_{\alpha\beta}+hc_\gamma$ and $\partial f_{\alpha\beta}/\partial c_{\alpha\beta}=\partial f_\gamma/\partial c_\gamma$ hold to solver tolerance; $\epsilon_m<10^{-6}$ with the RZ weight.
- Sanity: isothermal 950 K with uniform $c=0.225$ relaxes to $\phi\approx0.45$, $c_{\alpha\beta}\approx0.014$, $c_\gamma\approx0.54$.

#### Step 2.3: Multiphase thermotransport, unirradiated (Wen Fig. 2)

- Adds: `ADThermalMigration` with $M_T$ for all three phases (Eqs. 29–31).
- Fit: $s_W$, once. The γ zone ($T\gtrsim972$ K) is single-phase, so its Zr level is controlled mainly by the γ-phase $\beta_i$, which is exactly the parameter Wen does not report.
- Key data (Fig. 2, $1.3\times10^{5}$ s): $c\approx0.27$ in the γ zone ($r\lesssim800$ µm), $c\approx0.23$ outside, and no Zr well.
- Pass: γ-zone level within ±0.01 (fitted), outer level within ±0.01, minimum $c>0.2$.
- Sanity: Zr enriches at the center; with $M_T=0$ the Step 2.2 result is recovered.
- Plausibility of $s_W$: independent data put the γ mobility within about one order of magnitude of Part 1's. At Wen's γ-zone compositions ($c_\gamma=0.25$–0.47, 970–988 K), Li et al. (2010) gives 0.3–2.8 times the Part 1 value ($s=1$), and Hirschhorn 2019's DP-81 fit multiplies Li's by 10. Hirschhorn et al. (2020) Eq. 16 gives 7–14 times Part 1 at 1023 K and $c_\gamma=0.1$–0.2. An $s_W$ between about 1 and 30 is consistent with these; one of order $10^3$ is not.

#### Step 2.4: Irradiation-enhanced diffusion (Wen Fig. 3)

- Adds: $\xi(T)$ on $M_c$ (Eq. 24). Nothing else changes and $s_W$ stays frozen. This step directly tests Wen's central claim that $\xi$ alone turns Fig. 2 into Fig. 3.
- Key data (Fig. 3, $1.3\times10^{5}$ s):

| Region | $r$ (µm) | $c$ | $\phi$ |
| --- | --- | --- | --- |
| γ zone | 0–800 | ≈ 0.47–0.49, falling to ≈ 0.38 at 800 µm | 1 |
| γ + β | 850–1050 | ≈ 0.27 | ≈ 0.55 |
| Zr well (β-rich) | 1100–1450 | ≈ 0.07 | ≈ 0.25 |
| γ + α | 1500–2170 | ≈ 0.29 at 1500 µm, falling to ≈ 0.22 | ≈ 0.47 falling to ≈ 0.39 |

- DP-81 PIE (Hofman 1996 Fig. 7; the black squares in Wen Figs. 2–4): digitize all points.
- Pass: γ-zone $c$ within ±0.05; well minimum $c<0.12$ with its center within ±150 µm; zone edges within ±100 µm. Report the $L_2$ error against PIE alongside Wen's own curve (Fig. 4).
- Sanity: $\xi(900\ \mathrm K)=4.69$, $\xi(935\ \mathrm K)=3.07$, $\xi(988\ \mathrm K)=1.93$; with $C_v^r=0$ the Step 2.3 result is reproduced exactly.
- If Fig. 3 is not reached with $s_W$ frozen, that is the result of this step. Record it and do not refit.

#### Step 2.5: Interface width (Wen Table 2) and Part 2 closure

- Change: nothing physical. Rerun Step 2.4 with $\omega=5\times10^{5}$, $5\times10^{6}$ and $10^{7}$ J/m³ at $\kappa_\gamma=7.8\times10^{-4}$ J/m, on a 1 µm mesh (2170 elements, as Wen).
- Key data: Table 2 widths $l_{pb}=39.5$, 12.5 and 8.8 µm. These equal $\sqrt{\kappa_\gamma/\omega}$ to four digits, so they are the analytic width scale rather than a measurement. Wen Fig. 5B: the γ profile changes only locally across two orders of magnitude of $\omega$.
- Pass: the width measured with Eq. 35 ($\Delta\phi$ divided by $\max\lvert\mathrm d\phi/\mathrm dr\rvert$ at the γ / (γ + β) boundary) is proportional to $\sqrt{\kappa_\gamma/\omega}$ within ±10 %; zone positions move by less than 50 µm across the three $\omega$ values.
- Closure: convergence on Step 2.4 with 500, 1000 and 2000 elements and `dtmax` × 1, ½, ¼; zone edges change by less than 20 µm and zone concentrations by less than 0.005. Freeze the KKS framework for Part 3.

---

## 5. Part 3: Jung et al. (2025), Heat Conduction, Thermal Conductivity, Porosity and Sodium

Goal: make the temperature an unknown and add, one at a time, the heat source, the composition and temperature dependence of the conductivity, porosity from the phase fractions, and sodium infiltration of the pores.

### 5.1 Model and setup

Heat conduction is solved in its quasi-steady form. The thermal relaxation time $r_0^2\rho c_p/k$ is about 1 s, many orders of magnitude below the redistribution time, and Jung gives no $\rho$ or $c_p$.

$$
0=\nabla\cdot\left(k\nabla T\right)+\dot q,
\qquad
\dot q=\frac{\mathrm{LHGR}\,(1-B)}{\pi r_0^2},
\qquad
B=B_{\mathrm{end}}\frac{t}{t_{\mathrm{end}}}
$$

$$
k=P_c\,(1-P)^{1.5}\,k_0(T,c),
\qquad
k_0=a_1+a_2T+a_3T^2
$$

$$
a_1=17.5\,\frac{1-2.23W_{\mathrm{Zr}}}{1+1.61W_{\mathrm{Zr}}},
\qquad
a_2=1.54\times10^{-2}\,\frac{1+0.061W_{\mathrm{Zr}}}{1+1.61W_{\mathrm{Zr}}},
\qquad
a_3=9.38\times10^{-6},
\qquad
W_{\mathrm{Zr}}=\frac{c\,M_{\mathrm{Zr}}}{c\,M_{\mathrm{Zr}}+(1-c)\,M_U}
$$

$$
P_c=1-\frac{3P_{\mathrm{Na}}}{1-P}\cdot\frac{1-k_{\mathrm{Na}}/k_0}{1.163+1.837\,k_{\mathrm{Na}}/k_0},
\qquad
k_{\mathrm{Na}}=93-0.0581\,\theta+1.173\times10^{-5}\,\theta^2,
\quad
\theta=T-273.15
$$

$$
P=s_P(B)\left[0.4\,\phi+0.13\,(1-\phi)\right],
\qquad
s_P=\min\left(\frac{B}{0.01},1\right),
\qquad
P_{\mathrm{Na}}=P_{\mathrm{Na}}^{\mathrm{target}}\min\left(\max\left(\frac{B-0.01}{0.005},0\right),1\right)
$$

The phase-field part keeps the Part 2 structure, with Jung's Eqs. 13–24: $f_\gamma$ at $T_{\mathrm{ref}}=1010$ K, $f_{\alpha\beta}$ with $G_0^{\alpha\text{-U}}$, α/β mobilities with the factor $c_{\alpha\beta}$, no $\xi$, and $\kappa_c=0$.

| Item | Value | Source |
| --- | --- | --- |
| Geometry | 1-D radial, `coord_type = RZ`, $r_0=2.17$ mm, 1000 elements | Jung §3; see `moose/modules/heat_transfer/test/tests/code_verification/cylindrical_test_no1.i` for a 1-D RZ input |
| Boundary conditions | $T(r_0)=900$ K (DP-81); zero flux for $c$ and $\phi$ | Jung §3 |
| LHGR, $B_{\mathrm{end}}$ | 24 kW/m, 0.036 (DP-81); 25 kW/m, 0.077 (DP-11) | Jung §3 |
| Initial condition | $c=0.23$; $\phi$ linear from 0.5 (center) to 0.4 (surface); $T$ linear from 990 K to 900 K (initial guess only) | Jung §4, Fig. 3a |
| Interface | $\kappa=3l\sigma=7.5\times10^{-4}$ J/m, $w=6\sigma/l=6.0\times10^{3}$ J/m³ ($l=0.5$ mm, $\sigma=0.5$ J/m²); $L=2.5\times10^{-10}$ m³/(J s) | Jung Table 2 |
| α/β kinetics | $D_0$, $Q$, $\widetilde Q^*$ as in Part 2 | Jung Table 3 (= Hirschhorn 2019 Table 1) |
| γ kinetics | $\beta_U(c_\gamma,T)$, $\beta_{\mathrm{Zr}}(c_\gamma,T)$ from Li et al. (2010) (table below); $\widetilde Q_U^*=15.1$, $\widetilde Q_{\mathrm{Zr}}^*=105.4$ kJ/mol, as magnitudes in Eq. 18, so $M_T^\gamma$ has the form of Mohanty's Eq. 3 | Jung Eq. 17 and Hirschhorn 2019 §2.1 (DICTRA mobilities of Li et al.; Mohanty's average heats of transport) |
| Optimized kinetics ("optimized mobilities") | $\widetilde D_0^\beta\times2.722$, $M_c^\gamma\times9.981$, $\widetilde Q_\beta^*\times0.617$, $M_T^\gamma\times1.395$; total $M_c$ divided by 5 | Hirschhorn 2019 Table 2 and §2.1 (§5.2) |
| Sodium level | 0.05 (Baseline), 0.25, 0.50, as fractions | Jung §4 |
| Kernels | `ADHeatConduction` (AD property `k`), `ADBodyForce` with `function` $\dot q(t)$, `DirichletBC`; add `T` to `coupled_variables` of the KKS kernels and free energies | — |
| Burnup ramps | $s_P$ and $P_{\mathrm{Na}}$ as AuxVariables set by `FunctionAux` each step; time steps hit $B=0.01$ and $0.015$ exactly | — |

γ atomic mobilities (Li et al. 2010, Table 2, U–Zr subset; Pu terms dropped). Li's tabulated parameters are the DICTRA quantities $\Phi=-Q^S+RT\ln M^0$ in J/mol, so

$$
\beta_i=\frac{1}{RT}\exp\left(\frac{\Phi_i}{RT}\right),
\qquad
\Phi_i=x_U\Phi_i^{U}+x_{\mathrm{Zr}}\Phi_i^{\mathrm{Zr}}+x_Ux_{\mathrm{Zr}}\sum_{r=0}^{3}{}^{r}\Phi_i^{U,\mathrm{Zr}}\left(x_U-x_{\mathrm{Zr}}\right)^r,
\qquad
x_{\mathrm{Zr}}=c_\gamma,\ x_U=1-c_\gamma .
$$

The text layer of the PDF drops the minus signs; the values below were read from the rendered page.

| Parameter (J/mol) | $i=U$ | $i=\mathrm{Zr}$ |
| --- | --- | --- |
| $\Phi_i^{U}$ | $-107943.2-135.2T$ | $-68383.1-206.3T$ |
| $\Phi_i^{\mathrm{Zr}}$ | $RT\ln\left[3.6\times10^{-5}e^{-242400/RT}+5.3\times10^{-10}e^{-82300/RT}\right]$ | $RT\ln\left[3.2\times10^{-7}e^{-167746/RT}+1.14\times10^{-10}e^{-73741/RT}\right]$ |
| ${}^{0}\Phi_i^{U,\mathrm{Zr}}$ | $-23537.87-54.27T$ | $-277639.59+166T$ |
| ${}^{1}\Phi_i^{U,\mathrm{Zr}}$ | $312070.48-176.99T$ | $-81672.40+23.04T$ |
| ${}^{2}\Phi_i^{U,\mathrm{Zr}}$ | $-252306.31+187.94T$ | $-42148.42+36.49T$ |
| ${}^{3}\Phi_i^{U,\mathrm{Zr}}$ | $6979.57$ | $735973.01$ |

Single-point checks of the implementation:

- Tracer diffusivities $D_i^*=RT\beta_i$ against Li Fig. 2 at 1300 K: $\log_{10}D_U^*=-11.40$ in pure U, and $\log_{10}D_{\mathrm{Zr}}^*=-13.5$ for Zr in U.
- Against Hirschhorn 2019 Fig. 7, with the §5.1 modifiers:

| Point | $\widetilde D_\gamma=RT\,(9.981M_c^\gamma)/(V_mc_\gamma)$ | $\lvert\widetilde Q_\gamma\rvert=(1.395M_T^\gamma)\,V_m/M_c^\gamma$ | Hirschhorn Fig. 7 |
| --- | --- | --- | --- |
| center, $c_\gamma=0.5$, 980 K | $1.4\times10^{-14}$ m²/s | 51 kJ/mol | ≈ $2\times10^{-14}$ m²/s, −52 kJ/mol |
| surface, $c_\gamma=0.67$, 900 K | $7.2\times10^{-16}$ m²/s | 75 kJ/mol | ≈ $8\times10^{-16}$ m²/s, −73 kJ/mol |

Li's fit uses data from 973–1723 K, and Li reports poorer agreement near the γ miscibility gap (1023–1073 K). DP-81's γ region (≈ 940–990 K) sits at the low end of that range.

### 5.2 Paper issues and decisions

| Item | In the paper | Evidence | Decision |
| --- | --- | --- | --- |
| α/β switch, Eqs. 3–4 | the printed formula puts α above 935 K | Fig. 3a: β inside ($T>935$ K), α outside. Jung copies the switch from Hirschhorn 2019 Eq. 3, which is printed as $1.0-\tfrac12\tanh(\ldots)$ and so runs from 0.5 to 1.5; Hirschhorn's text says the fuel goes γ → β+γ → α+γ as temperature falls | physical order, as in Part 2 |
| $f_{\alpha\beta}$ temperature term, Eq. 23 | $+24.0T-23500.0$ (as in Wen) | Hirschhorn 2019 Eq. 18 prints $-29.333333T+24106.6667$. With that term U-22.5Zr never becomes single-phase γ below 1100 K, contradicting Hirschhorn's own reported transitions at 693 °C and 662 °C; Jung's term gives ≈ 970 K (Step 3.2) | Jung's term as printed |
| $g(\phi)$, Eq. 8 | $\phi^2(1-\phi^2)$, which is not a double well on $[0,1]$ | standard form in Wen and Hirschhorn | $\phi^2(1-\phi)^2$ |
| Appendix A | "ln(T)" terms | SGTE uses $T\ln T$. With the printed $\ln T$, U-22.5Zr never becomes single-phase γ between 900 and 1000 K, so Zone I could not form | $T\ln T$ |
| Eq. 31 layout | "3$P_{\mathrm{Na}}$" stacked over "$(1-P)$" | the reading $3P_{\mathrm{Na}}/(1-P)$ gives the ≈ −20 K center shift of Fig. 6b for 50 % versus 5 % sodium; $3P_{\mathrm{Na}}(1-P)$ gives ≈ −13 K (hand estimate) | $3P_{\mathrm{Na}}/(1-P)$, with $P_{\mathrm{Na}}$ the infiltration level as a fraction |
| Eqs. 34–35, $W_{\mathrm{Zr}}(c)$ | circular definition | Jung's own $k$ in Fig. 5c is 15–30 % higher than Eqs. 27–30 give with a standard conversion (case $P=0$, center: 40 versus ≈ 31 W/(m K)) | standard conversion (§5.1). Expected consequence: center temperature ≈ 13 K above Jung's (≈ 956 K instead of 943 K for $P=0$). Temperature criteria in Steps 3.3–3.5 are therefore ratios. Hirschhorn 2019 does not settle this, because it uses a different fresh-fuel conductivity (Kim et al. 2014) with one uniform correction factor |
| Coefficient names | text says $a_0,a_1,a_2$; equations say $a_1,a_2,a_3$ | — | $a_1,a_2,a_3$ |
| Heat equation | transient form, no $\rho$, $c_p$ | thermal time ≈ 1 s | quasi-steady |
| Geometry | "one-dimensional", $\dot q$ normalized by $\pi r_0^2$ | Jung's own $k(r)$ from Fig. 5c gives the Fig. 5b center temperatures in RZ (943 K and 976 K); a slab gives ≈ 1050 K | RZ |
| Total simulated time | not reported | — | $t_{\mathrm{end}}$ fitted in Step 3.2, then frozen |
| Units of $c$ | Table 1 says at%, but plotted values are atom fractions; the abstract's "0.17 at%" and "0.7 at%" | Fig. 5a center: 0.47 → 0.31 ($P=0$); Fig. 6a center: 0.47 → 0.405 (50 % sodium) | use the figure values (atom fraction) |
| Source of Table 2 and the "optimized mobilities" | attributed to "Hirschhorn et al. 2019"; the reference list gives the diffusion-couple paper (J. Nucl. Mater. 529, `Hirschhorn_2020`) | the values ($V_m$, $L$, $l$, $\sigma$, $\widetilde D_0$, $Q$, $\widetilde Q^*$) and the 1-D axisymmetric DP-81/DP-11 setup are those of the redistribution paper, J. Nucl. Mater. 523 (`Hirschhorn_2019.pdf`). The 2020 paper's isothermal polynomial mobilities would grow by about $10^{9}$ at Zone I compositions | Hirschhorn 2019 |
| Optimized modifiers | not printed by Jung | Hirschhorn 2019 Table 2: four modifiers fitted to DP-81 ($\widetilde D_0^\beta$ 2.722, $M_c^\gamma$ 9.981, $\widetilde Q_\beta^*$ 0.617, $M_T^\gamma$ 1.395). They act on Li's DICTRA γ mobilities, which §5.1 also uses | applied as listed in §5.1 |
| $M_c$ divided by ≈ 5 | not mentioned by Jung | Hirschhorn 2019 §2.1: the simplified free energies (Jung's Eqs. 23–24) raise the chemical driving force about fivefold, so Hirschhorn divides the chemical mobility by the same factor, and the Table 2 modifiers were optimized with this division in place | divide the total $M_c$ by 5. Because the factor acts on $M_c$ only, it changes the balance with $M_T$, not just the time scale |
| γ mobilities, Eqs. 19–20 | Arrhenius forms with one $\beta_0$ for U and Zr, value not printed | Jung's Eq. 17 cites Li et al. 2010, and Hirschhorn uses Li's composition- and temperature-dependent mobilities. These reproduce Hirschhorn's Fig. 7a interdiffusivity (§5.1) | Li et al. 2010 (§5.1) |
| γ heats of transport, Eq. 18 | values not printed; Hirschhorn uses "the average values" of Mohanty in a "+" form | entered as magnitudes (15.1 and 105.4 kJ/mol), they reproduce Hirschhorn Fig. 7b (51 → 75 kJ/mol against −52 → −73 kJ/mol; Hirschhorn reports negative values for Zr moving to the hot center). Mohanty's signed values would change the sign of $M_T^\gamma$ for $c_\gamma\gtrsim0.6$, which Fig. 7b does not show | magnitudes, i.e. Mohanty's Eq. 3 |
| Surface Zr drop ($r/r_0>0.95$), over-predicted surface α | acknowledged by Jung | — | reproduce, do not fix |

### 5.3 Steps

#### Step 3.1: Heat conduction alone

- Adds: $T$ as a nonlinear variable with $0=\nabla\cdot(k\nabla T)+\dot q$ and $T(r_0)=900$ K. $k(r)$ is prescribed from digitized Jung Fig. 5c (case $P=0$ and Baseline), with $\dot q$ at $B=0.036$, i.e. $\dot q=1.564\times10^{9}$ W/m³. No phase field.
- Key data (Fig. 5b): center temperature $T_c\approx943$ K ($P=0$) and $T_c\approx977$ K (Baseline).
- Expected (hand calculation with the digitized $k$): 943 K and 976 K. A slab would give ≈ 1050 K for the Baseline.
- Pass: within ±3 K.
- Sanity: with a constant $k=30$ W/(m K), $T(r)=T_s+\dot q\,(r_0^2-r^2)/(4k)$ exactly (center rise 61.4 K); the RZ-weighted integral of $\dot q$ equals $\mathrm{LHGR}(1-B)=23.14$ kW/m.

#### Step 3.2: Jung's phase-field model under a prescribed temperature

- Change (parameter set only, no new physics): the frozen Part 2 framework with Jung's free energies, interface parameters, mobilities and boundary/initial conditions, Hirschhorn's optimized modifiers and the division of $M_c$ by 5 (§5.1, §5.2), and without $\xi$. $T(r)$ is the digitized Baseline curve of Fig. 5b, held fixed in time.
- Thermodynamic pre-check (`lever.py` with Jung's energies): U-22.5Zr is single-phase γ above ≈ 970 K; the γ-side boundary is $c_\gamma\approx0.67$ at 900 K, ≈ 0.57 at 935 K and ≈ 0.51 at 950 K. A Zr-enriched Zone I ($c\approx0.5$) is therefore single γ down to ≈ 950 K, consistent with the Zone I edge at $r/r_0\approx0.40$, where the Baseline $T\approx958$ K. Hirschhorn 2019 reports the γ / β+γ and β+γ / α+γ transitions in DP-81 at about 966 K (693 °C) and 935 K (662 °C).
- Fit: $t_{\mathrm{end}}$, once, on the Zone I center concentration. Neither Jung nor Hirschhorn reports the simulated time. Jung says it is much shorter than DP-81's 284 EFPD ($2.45\times10^{7}$ s), which bounds the search.
- Sanity (before the run): the γ mobility materials reproduce the single-point values of §5.1 (Li Fig. 2 and Hirschhorn Fig. 7) to within 10 %. This checks Li's parameters, the modifiers and the heat-of-transport reading independently of $t_{\mathrm{end}}$.
- Key data (Fig. 3b, the Baseline curve of Fig. 5a): Zone I $c\approx0.47$ at the center, rising to ≈ 0.51 at $r/r_0\approx0.38$; Zone I/II boundary at $r/r_0\approx0.40$; Zr well $c\approx0.02$ at $r/r_0\approx0.53$–0.63; Zone II/III boundary at $r/r_0\approx0.64$; Zone III $c\approx0.22$ at $r/r_0\approx0.65$–0.93; drop to $c\approx0.02$ for $r/r_0>0.95$.
- Pass: zone boundaries within ±0.03 in $r/r_0$; Zone I $c$ within ±0.04; Zone III $c$ within ±0.03; well minimum $c<0.06$.

#### Step 3.3: Two-way coupling through $k_0(T,c)$ ($P=0$)

- Adds: the heat equation of Step 3.1 now uses $k=k_0(T,c)$ (no porosity, $P_{\mathrm{Na}}=0$), and $\dot q$ follows $B(t)$. Composition and phases change $k$, $k$ changes $T$, and $T$ changes $f$, $M_c$, $M_T$. $t_{\mathrm{end}}$ stays frozen. Jung's $P=0$ case keeps 5 % sodium, which is worth about 1 K here.
- Key data (Fig. 5, case $P=0$): $T_c\approx943$ K; center $c\approx0.31$, falling to ≈ 0.25 near $r/r_0\approx0.25$; no Zr well and no single-γ Zone I; $k\approx40$–42 W/(m K).
- Pass: center $c$ within ±0.05; no well. Report $T_c$ (≈ 956 K expected, §5.2) and the ratio of $k(r)$ to Fig. 5c (≈ 0.75–0.8 expected).
- Sanity: freezing $c$ and $\phi$ reproduces a Step 3.1-type steady solution with the same $k(r)$; Newton converges monolithically.

#### Step 3.4: Porosity from phase fractions

- Adds: $(1-P)^{1.5}$ in $k$. First constant $P=0.1$, 0.2, 0.3 from $t=0$ (no ramp), then the phase-dependent $P$ with the burnup ramp $s_P(B)$, which is the Baseline without sodium. $P_{\mathrm{Na}}=0$; Jung's 5 % sodium is worth 2–3 K here.
- Key data (Fig. 5). Temperature criterion: the rise $\Delta T_c=T_c-900$ K divided by the $P=0$ rise of Step 3.3. This ratio is insensitive to the $k_0$ offset in §5.2.

| Case | $T_c$ (K) | $\Delta T_c/\Delta T_c(P=0)$ | Center $c$ | Zone I ($\phi=1$ at center) |
| --- | --- | --- | --- | --- |
| $P=0.1$ | ≈ 951 | 1.19 | ≈ 0.36 | no |
| $P=0.2$ | ≈ 960 | 1.40 | ≈ 0.43 | no ($\phi\approx0.71$) |
| $P=0.3$ | ≈ 973 | 1.70 | ≈ 0.47 | yes |
| $P(\eta)$, Baseline without sodium | ≈ 977 | 1.79 | ≈ 0.47 | yes |

- Pass: ratios within ±0.08; center $c$ within ±0.05; presence or absence of Zone I correct in every case.
- Sanity: $P$ stays in $[0.13,0.4]$ after the ramp; the burnup kink at $B=0.01$ is resolved by the time steps.

#### Step 3.5: Sodium infiltration

- Adds: $P_c$ from Eqs. 31–32 with the sodium ramp between $B=0.01$ and $0.015$, for $P_{\mathrm{Na}}^{\mathrm{target}}=0.05$ (the full Jung Baseline), 0.25 and 0.50.
- Key data (Fig. 6), using the same temperature ratio as Step 3.4, here relative to the 5 % case:

| Sodium | $T_c$ (K) | $\Delta T_c/\Delta T_c(5\,\%)$ | Center $c$ | Center $\phi$ | Center $k$ (W/(m K)) |
| --- | --- | --- | --- | --- | --- |
| 5 % | ≈ 977 | 1 | ≈ 0.47 | 1.0 | ≈ 15 |
| 25 % | ≈ 962 | 0.81 | ≈ 0.47 | ≈ 0.82 | ≈ 21 |
| 50 % | ≈ 956 | 0.73 | ≈ 0.405 | ≈ 0.65 | ≈ 27 |

- Baseline (5 %) against Fig. 3b: the same zone criteria as Step 3.2. Report the $L_2$ error against the DP-81 PIE points (Chang 1994, Fig. II.2, as plotted by Jung).
- Pass: temperature ratios within ±0.05; the center-$c$ drop from 5 % to 50 % equals 0.065 ± 0.03; center $\phi$ within ±0.1; Baseline zone boundaries within ±0.05 in $r/r_0$.
- Sanity: $P_c=1$ when $P_{\mathrm{Na}}=0$, and $P_c>1$ everywhere because $k_{\mathrm{Na}}>k_0$ (at 900 K, $k_{\mathrm{Na}}/k_0\approx2$). This step also confirms the Eq. 31 reading of §5.2.

#### Step 3.6: Second pin DP-11 and Part 3 closure

- Change: operating conditions only: $T_s=910$ K, LHGR 25 kW/m, $B_{\mathrm{end}}=0.077$. $t_{\mathrm{end}}$ is scaled by 7.7/3.6 (the same burnup rate as DP-81, an assumption). Jung gives no separate DP-11 initial conditions, so the DP-81 ones are used.
- Key data (Fig. 4): Zone I $c\approx0.50$–0.53 for $r/r_0<0.5$; zone boundaries at $r/r_0\approx0.52$ and ≈ 0.72; well $c\approx0.02$ at $r/r_0\approx0.58$–0.70; Zone III $c\approx0.20$ at $r/r_0\approx0.73$–0.88; $\eta_\alpha\to1$ at the surface. PIE: Zone I ≈ 0.48, Zone III ≈ 0.20.
- Pass: zone boundaries within ±0.05 in $r/r_0$; Zone I $c$ within ±0.05. This is a prediction with no fitting.
- Closure: convergence on the DP-81 Baseline with 500, 1000 and 2000 elements and `dtmax` × 1, ½, ¼ (the burnup kinks must always be hit); key data change by less than one fifth of their tolerance. Write `p3_jung/README.md`.

---

## 6. Data to Digitize

Store one CSV per curve in `problem/uzr/data/`, with a header recording the paper, figure, axis units and an estimated reading error (WebPlotDigitizer or similar).

| Paper | Figure | Curves | Used in |
| --- | --- | --- | --- |
| Mohanty | 6 | $c(x)$ at 5, 10, 30 days | 1.3 |
| Mohanty | 7 | $c(0,t)$ | 1.3 |
| Mohanty | 8a | isothermal (dashed) and gradient (solid) curves | 1.2, 1.4 |
| Mohanty | 8b | gradient (solid) curve | 1.4 |
| Wen | 2 | $T(r)$, $c(r)$, $\phi(r)$ | 2.1–2.3 (and $T(r)$ for all of Part 2) |
| Wen | 3 | $c(r)$, $\phi(r)$, α and β fractions | 2.4 |
| Wen | 2–4 | DP-81 PIE points (Hofman 1996) | 2.4 |
| Jung | 3b / 5a | Baseline $c(r)$; DP-81 PIE points (Chang 1994) | 3.2, 3.5 |
| Jung | 5a–5d | $c$, $T$, $k$, $\phi$ for $P=0$, 0.1, 0.2, 0.3 and Baseline | 3.1, 3.3, 3.4 |
| Jung | 6a–6d | $c$, $T$, $k$, $\phi$ for 5, 25, 50 % sodium | 3.5 |
| Jung | 4 | DP-11 $c(r)$ and phase fractions; DP-11 PIE points | 3.6 |

`scripts/compare.py` reads an Exodus file and the matching CSV files, samples the simulation at the key points of the step, and prints the key-data table (paper value, simulation value, tolerance, pass/fail) for the step `README.md`.

---

## 7. References and Files

Papers (in `palm/doc/`):

- Mohanty, Bush, Okuniewski, Sohn, J. Nucl. Mater. 414 (2011) 211–216: `Mohanty_2011`
- Wen, Liu, Yun, Qian, Front. Mater. 9 (2022) 950785: `Wen_2022`
- Jung, Kim, Chang, Ann. Nucl. Energy 224 (2025) 111679: `Jung_2025`
- Hirschhorn, Tonks, Aitkaliyeva, Adkins, J. Nucl. Mater. 523 (2019) 143–156: `Hirschhorn_2019.pdf`. The model Jung builds on: source of Jung's Table 2, the heats of transport, the optimized modifiers and the division of $M_c$ by 5 (§5.1, §5.2)
- Hirschhorn, Tonks, Aitkaliyeva, Adkins, J. Nucl. Mater. 529 (2020) 151929, isothermal β/γ diffusion couple at 1023 K: `Hirschhorn_2020`. Used only for the γ-mobility plausibility check in Step 2.3. It is not a reproduction target, because its β and γ free energies are not printed
- Li, Hu, Cui, Zhong, Chang, Li, Zhou, J. Nucl. Mater. 407 (2010) 220–227: `Li_2010.pdf`. DICTRA atomic mobilities of bcc U–Pu–Zr; its U–Zr subset gives Jung's γ kinetics (§5.1)
- Secondary sources named in the plan: Dinsdale, Calphad 15 (1991) 317 (SGTE unary data); Chevalier et al., Calphad 28 (2004) 15 (bcc U–Zr interaction parameters); Hofman et al., J. Nucl. Mater. 227 (1996) 277 (DP-81 PIE, heats of transport)

Code:

- `palm/include/kernels/ADThermalMigration.h`, `palm/src/kernels/ADThermalMigration.C`, `palm/test/tests/kernels/ad_thermal_migration/`
- `palm/problem/uzr/s1_mohanty/` (`s1_mohanty.i`, `ReadMe.md`, `plot_exodus_cx_profile.py`)
- MOOSE references: `moose/modules/phase_field/test/tests/KKS_system/kks_example_split.i` (KKS split form), `moose/modules/heat_transfer/test/tests/code_verification/cylindrical_test_no1.i` (1-D RZ heat conduction), `moose/modules/phase_field/src/kernels/SoretDiffusion.C`

Not covered by this plan: two-dimensional $r$–$z$ geometry, U–Pu–Zr, and any physics not contained in the three papers.

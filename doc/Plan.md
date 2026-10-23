# U–Zr Constituent Redistribution under a Temperature Gradient: Stepwise MOOSE Reproduction Plan

> Code: PALM (MOOSE application; phase-field and heat-transfer modules), split Cahn–Hilliard form with KKS multiphase thermodynamics  
> Part 1: Hirschhorn et al. (2019), the base model, with its data sources (Li 2010, Quaini 2018, Kim 2014, Galloway 2015) and the companion diffusion-couple study (Hirschhorn 2020)  
> Part 2: Jung et al. (2025), porosity and sodium infiltration in the thermal conductivity  
> Part 3: Wen et al. (2022), irradiation-enhanced diffusion  
> Starting point: the PALM kernel `ADThermalMigration` (5 tests passing) and the archived case `problem/uzr/p1_mohanty/step1_1_kernel_swap`

---

## 1. Strategy

### 1.1 Scope and order of the papers

- Hirschhorn 2019 is the base. Jung 2025 copies its model almost verbatim, and Wen 2022 builds on it. It was solved in MARMOT on MOOSE, couples composition, phases and temperature, validates on two EBR-II pins (DP-81 and DP-11), and publishes its optimized parameters.
- Hirschhorn 2019 has no isothermal result. Li 2010, the source of its γ mobilities, and Hirschhorn 2020, an isothermal β/γ diffusion couple solved with the same model by the same group, supply the single-phase and two-phase diffusion checks that come before the temperature gradient. Quaini 2018, Kim 2014 and Galloway 2015 supply inputs that Hirschhorn 2019 cites but does not print.
- Jung 2025 is Hirschhorn 2019 with one block replaced: the uniform, burnup-dependent conductivity factor becomes a porosity- and sodium-dependent conductivity. Part 2 starts from the end state of Part 1 and changes only that block and Jung's setup.
- Wen 2022 adds one physical process, the irradiation-enhancement factor $\xi$ of the chemical mobility. Its own setup is under-specified (no γ mobility prefactor, no $L_0$, a 1.5-day run with explicit finite differences), so a literal reproduction would mostly test Wen's unreported choices. Part 3 asks Wen's actual question on the common model: can $\xi$ replace Hirschhorn's fitted, roughly tenfold mobility enhancement over the real irradiation time?
- Mohanty 2011 is not reproduced as a part. Its printed mobilities make its own Figs. 6–7 evolve 24–38 times too fast; this was established with the archived `step1_1_kernel_swap` case, whose steady state matches an independent ODE solution. Its γ thermal-mobility form, Eq. (3), is still used: it is the form Hirschhorn 2019 needs to reproduce its own Fig. 7b (§3.3).

Each part starts from the end state of the previous part, so the final model contains every process of the three papers.

### 1.2 Roadmap

Each row is one small, runnable step that adds one physical process (rule 1 of §1.3).

| Step | Physics added (or change made)                                                                    | Key data used for validation                       | Fitted                             |
| ---- | ------------------------------------------------------------------------------------------------- | -------------------------------------------------- | ---------------------------------- |
| 1.1  | Thermodynamics: Hirschhorn's simplified energies and Quaini's full energies (no PDE)              | Hirschhorn 2019 Fig. 2; Hirschhorn 2020 tie-line   | —                                  |
| 1.2  | γ kinetics: Li atomic mobilities, γ heat of transport (no PDE)                                    | Li Figs. 2–3; Hirschhorn 2019 Fig. 7               | —                                  |
| 1.3  | Single-phase γ interdiffusion (isothermal couple)                                                 | Li Fig. 5b (Muller data)                           | —                                  |
| 1.4  | Two-phase β/γ KKS diffusion with a moving interface (isothermal couple)                           | Hirschhorn 2020 Table 1, Fig. 2                    | —                                  |
| 1.5  | Three-phase KKS model in a prescribed radial temperature field                                    | Hirschhorn 2019 Fig. 3, phases                     | —                                  |
| 1.6  | + thermotransport in α, β, γ                                                                      | Hirschhorn 2019 Fig. 3, Zr                         | —                                  |
| 1.7  | + heat conduction: Kim conductivity, burnup factor $k_c$, composition-dependent power             | Hirschhorn 2019 Fig. 3; Fig. 6 center temperatures | —                                  |
| 1.8  | + Hirschhorn's four optimized modifiers (DP-81)                                                   | Hirschhorn 2019 Figs. 4, 6, 7; DP-81 PIE           | $t_{\mathrm{end}}$, only if needed |
| 1.9  | DP-11 (no new physics); closure                                                                   | Hirschhorn 2019 Fig. 5; DP-11 PIE                  | —                                  |
| 2.1  | Jung's heat source and setup, conductivity prescribed from Jung's own figure                      | Jung Fig. 5b                                       | —                                  |
| 2.2  | Jung's $k_0(T,c)$ replaces Kim's conductivity times $k_c$ ($P=0$)                                 | Jung Fig. 5, case $P=0$                            | $t_{\mathrm{end}}$                 |
| 2.3  | + porosity from phase fractions                                                                   | Jung Fig. 5                                        | —                                  |
| 2.4  | + sodium infiltration (full Jung Baseline)                                                        | Jung Figs. 6, 3b; DP-81 PIE                        | —                                  |
| 2.5  | DP-11 (no new physics); closure                                                                   | Jung Fig. 4; DP-11 PIE                             | —                                  |
| 3.1  | Irradiation-enhanced diffusion $\xi(T)$ replaces the fitted mobility modifiers (DP-81, real time) | Wen Fig. 3; DP-81 PIE; Step 2.4                    | —                                  |
| 3.2  | DP-11 with $\xi(T)$ (no new physics)                                                              | DP-11 PIE; Step 2.5                                | —                                  |

### 1.3 Rules for every step

1. One change per step. Each step adds one physical process. Steps 1.3–1.5 and 2.1 also switch to the setup of the experiment they reproduce, and Steps 1.9, 2.5 and 3.2 change only the pin. From Step 1.5 on, consecutive inputs within a part differ only by the block that implements the new process.
2. Validation by key data. Each step is judged against one to four numbers or features read from one figure or table, each with a stated tolerance, plus at most two sanity checks (a limiting case or a conservation check). Mesh and time-step convergence is checked once per part, in its closing step.
3. Printed equations come first. A printed equation is replaced only when it contradicts the paper's own figures or is mathematically impossible. Each replacement is one row in the part's "Paper issues and decisions" table, with the evidence.
4. One frozen fit per part. Part 1 uses Hirschhorn's published modifiers and needs no fit; its slot is reserved for the simulated time $t_{\mathrm{end}}$, which Hirschhorn does not state, and is used only if Step 1.8 misses its key data. Part 2 fits Jung's unreported $t_{\mathrm{end}}$ once in Step 2.2. Part 3 has no fit, because its purpose is to test whether $\xi$ removes the need for one.
5. No physics beyond the papers. There is no defect rate-theory PDE, no pore nucleation or growth, and no mechanics, FCCI or fission gas. Porosity and sodium act only through the thermal conductivity, and irradiation only through $\xi$.
6. Units are the same in all parts (§2.4).
7. Values marked "≈" in this plan were read by eye from the PDF figures. They are good enough to plan with, but must be replaced by digitized values (§6) before any error is computed. Values marked "expected" were precomputed from the equations of this plan and indicate what a correct implementation should give.

### 1.4 Standard checks and outputs

- Mass conservation: $\epsilon_m=\lvert\int_\Omega c\,\mathrm dV-\int_\Omega c_0\,\mathrm dV\rvert/\int_\Omega c_0\,\mathrm dV<10^{-6}$. With `coord_type = RZ`, MOOSE integral postprocessors include the $2\pi r$ weight automatically.
- Bounds: $0<c<1$, $0\le\phi\le1$ (small overshoot allowed and recorded), $T>0$, $k>0$.
- KKS consistency (from Step 1.4 on): $c=(1-h)c_{\alpha\beta}+hc_\gamma$ and equal phase chemical potentials to solver tolerance.
- Every step directory contains one input file and a `README.md` of at most one page: what changed from the previous step, the key-data table (paper value, this run, tolerance, pass/fail), the sanity checks, and the run command.

---

## 2. Common Framework

### 2.1 Governing equations and sign convention

The conserved variable is the Zr atom fraction $c$; from Step 1.4 on, the γ order parameter $\phi$ and the KKS phase concentrations are added; from Step 1.7 on, the temperature is an unknown.

$$
\frac{\partial c}{\partial t}=\nabla\cdot\left(M_c\nabla w\right)-\nabla\cdot\left(M_T\frac{\nabla T}{T}\right),
\qquad
w=\frac{\delta F}{\delta c},
\qquad
\frac{\partial\phi}{\partial t}=-L\frac{\delta F}{\delta\phi},
\qquad
0=\nabla\cdot\left(k\nabla T\right)+\dot q .
$$

PALM's convention is the flux $\mathbf J=-M_c\nabla w+M_T\nabla T/T$: a positive $M_T$ moves Zr toward the hot side. The heat equation is solved in its quasi-steady form, because the thermal relaxation time $r_0^2\rho c_p/k\approx1$ s is many orders of magnitude below the redistribution time and none of the papers gives $\rho$ and $c_p$.

The papers do not share a sign convention, and Hirschhorn and Jung print one that their own results contradict:

| Source | Printed | What the paper's own results require | `M_T` in PALM |
| --- | --- | --- | --- |
| Mohanty 2011, Eqs. (1), (3) | $J_T=+M_Q\nabla T/T$ with $M_Q=c(1-c)[\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*]$ | consistent | γ phase: this expression |
| Hirschhorn 2019 Eq. (1); Jung 2025 Eq. (1) | $-M_T\nabla T/T$ inside the divergence, which reads as $M_T>0$ → hot side | standard convention, positive heat of transport → cold side (evidence in §3.3, row 1) | α, β: the printed $M_T^i$ with a minus sign; γ: Mohanty's expression |
| Wen 2022 Eq. (4) | same as Hirschhorn | not needed: Part 3 takes only $\xi$ from Wen | — |

### 2.2 MOOSE objects

Variable `w` carries the mass balance and variable `c` carries the definition of `w`, as in the standard MOOSE split form.

| Term | Step 1.3 (single phase) | Steps 1.4 onward (KKS) |
| --- | --- | --- |
| $\partial c/\partial t$ | `ADCoupledTimeDerivative` (`variable = w`, `v = c`) | same |
| $\nabla\cdot(M_c\nabla w)$ | `ADSplitCHWRes` | same |
| $-\nabla\cdot(M_T\nabla T/T)$ | — | `ADThermalMigration` (from Step 1.6) |
| $w=\delta F/\delta c$ | `ADSplitCHParsed` with $\kappa_c=0$ | `KKSSplitCHCRes`, `KKSPhaseConcentration`, `KKSPhaseChemicalPotential` |
| $\phi$ equation | — | `TimeDerivative`, `KKSACBulkF`, `KKSACBulkC`, `ACInterface` |
| Temperature | uniform constant | uniform constant (1.4); AuxVariable (1.5–1.6); from 1.7: `ADHeatConduction`, `ADMatHeatSource`, `DirichletBC` at $r_0$ |

The KKS kernels are not AD, so the phase free energies, $h(\phi)$, $g(\phi)$ and $L$ are `DerivativeParsedMaterial` objects, with `T` added to their `coupled_variables` once temperature is a variable. `M_c`, `M_T`, `k` and $\dot q$ are `ADParsedMaterial` objects of $(c,\phi,c_{\alpha\beta},c_\gamma,T)$ that write $h(\phi)$ and $h_{\alpha\beta}(T)$ inline, which avoids a name clash between AD and non-AD properties. Burnup-dependent factors enter as AuxVariables set by `FunctionAux`. Solver: monolithic `NEWTON` with LU (the problems are 1-D), `IterationAdaptiveDT`.

### 2.3 Why `ADThermalMigration` replaces `SoretDiffusion`

MOOSE's `SoretDiffusion` discretizes $\int_\Omega D\,Q\,c/(k_BT^2)\,\nabla T\cdot\nabla\psi\,\mathrm dV$ with $k_B$ hard-coded in eV/K. Expressing the papers' flux $M_T\nabla T/T$ with it needs $D=M_ck_BT/c$ and $Q=-M_T/M_c$. That mapping has three problems that grow in the multiphase model:

1. Incomplete Jacobian. The kernel differentiates only its explicit dependence on $c$ and $T$, not the derivatives of $D$ and $Q$. Here $M_T$ depends on $\phi$, $c_{\alpha\beta}$, $c_\gamma$ and $T$, and those Jacobian blocks would be missing entirely.
2. Fixed units. A hard-coded $k_B$ in eV/K forces eV-based units, while all three parts use J/mol and J/m³.
3. Artificial singularity. Splitting $M_T$ into $D\,Q\,c$ puts $1/c$ into $D$, which is singular as $c\to0$, a limit the Zr-depleted β-rich bathtub ($c\approx0.01$) actually reaches.

`ADThermalMigration` (in `include/kernels/` and `src/kernels/`) discretizes

$$
R_{\mathrm{TM}}=-\int_\Omega\frac{M_T}{T}\,\nabla T\cdot\nabla\psi\,\mathrm dV .
$$

It reads $M_T$ directly as an AD material property, accepts either sign (negative for α and β here), contains no physical constants, and obtains the full Jacobian through automatic differentiation. Its tests in `test/tests/kernels/ad_thermal_migration/` cover an exact steady state, a reversed gradient, a `PetscJacobianTester` run and equivalence with `SoretDiffusion`. The archived case `p1_mohanty/step1_1_kernel_swap` reproduces Mohanty's printed equations exactly, which was confirmed against an independent ODE steady state. MOOSE's own `ADCHSoretMobility` has the same weak form with the opposite sign; the PALM kernel is kept because its sign matches the convention of §2.1.

### 2.4 Units and conversions

| Quantity                                                                                    | Unit in PALM             | Note                                                                                     |
| ------------------------------------------------------------------------------------------- | ------------------------ | ---------------------------------------------------------------------------------------- |
| Length, time, temperature                                                                   | m, s, K                  | °C + 273.15                                                                              |
| Molar free energy $f_{\mathrm m}$, heats of transport $\widetilde Q^*$, activation energies | J/mol                    | the code uses $f=f_{\mathrm m}/V_m$ in J/m³                                              |
| Molar volume $V_m$                                                                          | m³/mol                   | $1.249\times10^{-5}$ (Hirschhorn 2019, Jung); $1.309\times10^{-5}$ (Hirschhorn 2020)     |
| Atomic mobility $\beta_i$                                                                   | m² mol/(J s)             | Li 2010: $\beta_i=\exp(Q_i/RT)/(RT)$, tracer $D_i^*=RT\beta_i$                           |
| $M_c$, $M_T$, $L$                                                                           | m⁵/(J s), m²/s, m³/(J s) | —                                                                                        |
| Conductivity, power density                                                                 | W/(m K), W/m³            | —                                                                                        |
| Burnup $B$                                                                                  | fraction                 | 3.6 at.% is written 0.036                                                                |
| Mohanty's per-atom heats                                                                    | → J/mol                  | $+2.5\times10^{-20}$ J = +15.06 kJ/mol (U); $-17.5\times10^{-20}$ J = −105.4 kJ/mol (Zr) |

### 2.5 Directory layout

```text
palm/problem/uzr/
├── data/                        digitized paper curves, one CSV per curve (§6)
├── scripts/                     plot_exodus_cx_profile.py, compare.py, lever.py, mobility.py, matano.py
├── p0_mohanty/                  archived SoretDiffusion runs (SD_Fitted, SD_non_Fitted, ADTM_non_Fitted)
├── p1_hirschhorn/
│   ├── step1_1_thermo/          Python only: phase boundaries and lever rule
│   ├── step1_2_gamma_kinetics/  Python only: Li mobilities, D*, interdiffusivity, heat of transport
│   ├── step1_3_gamma_couple/
│   ├── step1_4_beta_gamma_couple/
│   ├── step1_5_kks_prescribed_T/
│   ├── step1_6_thermotransport/
│   ├── step1_7_heat/
│   ├── step1_8_dp81/
│   ├── step1_9_dp11/
│   └── README.md                filled key-data tables of Part 1
├── p2_jung/
│   ├── step2_1_heat/  step2_2_k0/  step2_3_porosity/  step2_4_sodium/  step2_5_dp11/
│   └── README.md
└── p3_wen/
    ├── step3_1_red_dp81/  step3_2_red_dp11/
    └── README.md
```

The empty placeholder directories `p1_mohanty/step1_2…step1_4`, `p2_wen/` and `p3_jung/` are replaced by this layout. No C++ beyond `ADThermalMigration` is planned. A C++ material is written only if a parsed expression becomes unreadable (the likely candidates are the Li mobilities and Jung's conductivity). MOOSE sources are not modified.

---

## 3. Part 1: Hirschhorn et al. (2019), the Base Model

Goal: build the complete three-phase model step by step, from thermodynamics and isothermal diffusion to thermotransport and heat conduction, and reproduce Hirschhorn's DP-81 and DP-11 results.

### 3.1 Model

Free energies (J/mol; Hirschhorn Eqs. 18–19 with the corrections of §3.3; SGTE end-members in Appendix A):

$$
f_{\alpha\beta}=(1-c)\,G^{\alpha\text{-U}}+c\left(38000+G^{\alpha\text{-Zr}}\right)+RT\left[(1-c)\ln(1-c)+c\ln c\right]-7000\,c(1-c)+24.0\,T-23500.0
$$

$$
f_{\gamma}=(1-c)\,G^{\gamma\text{-U}}+c\,G^{\beta\text{-Zr}}+RT_{\mathrm{ref}}\left[(1-c)\ln(1-c)+c\ln c\right]+c(1-c)\left[28369.0639-14.741714\,T_{\mathrm{ref}}+4992.02639(1-2c)+2609.10711(1-2c)^2\right],
\qquad T_{\mathrm{ref}}=1010\ \mathrm K
$$

Phase field (Hirschhorn Eqs. 2–12): $h=\phi^3(6\phi^2-15\phi+10)$, $g=\phi^2(1-\phi)^2$, $F=\int\{[(1-h)f_{\alpha\beta}+hf_\gamma]/V_m+wg+\tfrac{\kappa}{2}\lvert\nabla\phi\rvert^2\}\,\mathrm dV$, $\kappa=3l\sigma$, $w=6\sigma/l$ with $l=500$ µm and $\sigma=0.5$ J/m², $L=2.5\times10^{-10}$ m³/(J s), and the KKS constraints. The α/β switch puts β above 662 °C:

$$
h_{\alpha\beta}(T)=\frac12\left[1+\tanh\frac{T-935.15}{2}\right],
\qquad
\eta_\gamma=\phi,\quad\eta_\beta=h_{\alpha\beta}(1-\phi),\quad\eta_\alpha=(1-h_{\alpha\beta})(1-\phi).
$$

Mobilities (Hirschhorn Eqs. 13–17 in PALM's sign convention). The modifiers $s$ are 1 until Step 1.8:

$$
M_c=\frac15\left[(1-h_{\alpha\beta})(1-h)M_c^{\alpha}+h_{\alpha\beta}(1-h)M_c^{\beta}+h\,M_c^{\gamma}\right],
\qquad
M_T=(1-h_{\alpha\beta})(1-h)M_T^{\alpha}+h_{\alpha\beta}(1-h)M_T^{\beta}+h\,M_T^{\gamma}
$$

$$
M_c^{i}=\frac{c_{\alpha\beta}V_m}{RT}\,s_{D,i}\widetilde D_0^{i}e^{-Q_i/RT},
\qquad
M_T^{i}=-\,s_{Q,i}\widetilde Q_i^{*}\,\frac{c_{\alpha\beta}}{RT}\,s_{D,i}\widetilde D_0^{i}e^{-Q_i/RT}
\qquad(i=\alpha,\beta)
$$

$$
M_c^{\gamma}=V_m\,c_\gamma(1-c_\gamma)\left[c_\gamma\beta_U+(1-c_\gamma)\beta_{\mathrm{Zr}}\right],
\qquad
M_T^{\gamma}=s_{Q,\gamma}\,c_\gamma(1-c_\gamma)\left[\beta_U\widetilde Q_U^{*}-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^{*}\right],
\qquad
\beta_i=\frac{s_\beta}{RT}\exp\left(\frac{Q_i(c_\gamma,T)}{RT}\right)
$$

$Q_i$ follows Li 2010 Eq. (2) and Table 2 for the U–Zr binary:

$$
Q_i=(1-c)\,Q_i^{U}+c\,Q_i^{\mathrm{Zr}}+c(1-c)\sum_{r=0}^{3}{}^{r}Q_i^{U,\mathrm{Zr}}(1-2c)^r .
$$

| Li Table 2 (J/mol) | U | Zr |
| --- | --- | --- |
| $Q_i^{U}$ | $-107943.2-135.2T$ | $-68383.1-206.3T$ |
| $Q_i^{\mathrm{Zr}}$ | $RT\ln[3.6\times10^{-5}e^{-242400/RT}+5.3\times10^{-10}e^{-82300/RT}]$ | $RT\ln[3.2\times10^{-7}e^{-167746/RT}+1.14\times10^{-10}e^{-73741/RT}]$ |
| ${}^{0}Q_i^{U,\mathrm{Zr}}$ | $-23537.87-54.27T$ | $-277639.59+166T$ |
| ${}^{1}Q_i^{U,\mathrm{Zr}}$ | $312070.48-176.99T$ | $-81672.40+23.04T$ |
| ${}^{2}Q_i^{U,\mathrm{Zr}}$ | $-252306.31+187.94T$ | $-42148.42+36.49T$ |
| ${}^{3}Q_i^{U,\mathrm{Zr}}$ | $6979.57$ | $735973.01$ |

Heat conduction (from Step 1.7; Kim 2014 Eqs. 1, 2, 13, 14; Galloway 2015 Fig. 12):

$$
k=k_c(B)\,k_{\mathrm{U\text{-}Zr}}(T,x),
\qquad
k_{\mathrm{U\text{-}Zr}}=\left(1-\sqrt{1-x}\right)k_{\mathrm{Zr}}+\sqrt{1-x}\left[(1-x)k_U+x\,k_{c,\mathrm{Zr}}\right],
\qquad
x=W_{\mathrm{Zr}}=\frac{c\,M_{\mathrm{Zr}}}{c\,M_{\mathrm{Zr}}+(1-c)M_U}
$$

$$
k_U=21.73+1.591\times10^{-2}T+5.907\times10^{-6}T^2,
\qquad
k_{\mathrm{Zr}}=8.853+7.082\times10^{-3}T+2.533\times10^{-6}T^2+2992/T
$$

$$
k_{c,\mathrm{Zr}}=-102.0+200.1x-109.2x^2+9.435\times10^{-3}T+3.459\times10^{-5}T^2-0.02093\,xT
$$

$$
\dot q=\frac{q'(1-B)}{\pi r_0^2}\left(1.248345-1.128840\,c\right),
\qquad
B=B_{\mathrm{end}}\frac{t}{t_{\mathrm{end}}}
$$

$k_c(B)$ is piecewise linear: 1 at $B=0$, 0.50 at $B=0.0075$, 0.70 at $B=0.015$, then constant. The power factor equals 1 at $c=0.220$, so with the conserved mean composition the total power stays $q'(1-B)$.

### 3.2 Setup and parameters

| Item | Value | Source |
| --- | --- | --- |
| Geometry (Steps 1.5–1.9) | 1-D radial, `coord_type = RZ`, $0\le r\le r_0=2170$ µm, 1000 elements | Hirschhorn §3 ("1D, axisymmetric") |
| Boundary conditions | zero flux for all composition and phase variables; $T(r_0)=T_s$ | Hirschhorn §3 |
| DP-81 | $q'=24$ kW/m, $T_s=627$ °C = 900.15 K, $B_{\mathrm{end}}=0.036$, $t_{\mathrm{end}}=284$ EFPD | Hirschhorn §2.2; EFPD from Hofman 1996 as quoted by Jung |
| DP-11 | $q'=25$ kW/m, $T_s=637$ °C = 910.15 K, $B_{\mathrm{end}}=0.077$, $t_{\mathrm{end}}\approx583$ EFPD | Hirschhorn §2.2; time scaled at the same burnup per unit power (derived) |
| Initial state | $c_0=0.22$; $\phi$, $c_{\alpha\beta}$, $c_\gamma$ from the lever rule (Step 1.1) at the initial $T$: the prescribed $T$ in Steps 1.5–1.6, the beginning-of-life steady solution ($k_c=1$) from Step 1.7 on | Hirschhorn §4 (22 at.%), Fig. 6a |
| Time step | adaptive, 1 s to $10^5$ s; steps hit $B=0.0075$ and $0.015$ exactly | Hirschhorn §3 |
| α, β kinetics | $\widetilde D_0^{\alpha}=2.0\times10^{-7}$, $\widetilde D_0^{\beta}=5.7\times10^{-5}$ m²/s; $Q_\alpha=1.7\times10^{5}$, $Q_\beta=1.8\times10^{5}$ J/mol | Hirschhorn §2.1 (Hofman) |
| α, β heats of transport | $\widetilde Q_\alpha^*=2.0\times10^{5}$, $\widetilde Q_\beta^*=4.5\times10^{5}$ J/mol | Hirschhorn §2.1 (Kim 2006, Galloway 2015) |
| γ heats of transport | $\widetilde Q_U^*=+15.06$, $\widetilde Q_{\mathrm{Zr}}^*=-105.4$ kJ/mol | "average values" of Mohanty 2011, confirmed by Fig. 7b (§3.3) |
| Optimized modifiers (Step 1.8) | $s_{D,\beta}=2.722$, $s_\beta=9.981$, $s_{Q,\beta}=0.617$, $s_{Q,\gamma}=1.395$ | Hirschhorn Table 2 |

### 3.3 Paper issues and decisions

| Item | Printed | Evidence | Decision |
| --- | --- | --- | --- |
| Sign of the thermal term, Eq. (1) | $M_T>0$ moves Zr to the hot side | Fig. 7b gives a negative $\widetilde Q_\gamma^*=M_T^\gamma V_m/M_c^\gamma$ (−50 to −73 kJ/mol) while Zr enriches at the center; the text says the γ and β+γ regions "possess heats of transport with opposite signs" and that larger $\widetilde Q_\beta^*$ moves Zr outward (Fig. 11); Galloway 2015, the source of $\widetilde Q_\alpha^*$ and $\widetilde Q_\beta^*$, states that positive heats of transport drive Zr to the cold side | standard convention, implemented as the minus sign in $M_T^{\alpha,\beta}$ (§3.1) |
| γ thermal mobility, Eq. (17) | $c_\gamma(1-c_\gamma)(\widetilde Q_U^*\beta_U+\widetilde Q_{\mathrm{Zr}}^*\beta_{\mathrm{Zr}})$ with Mohanty's averages | With Li's mobilities, only Mohanty's form $\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*$ times 1.395 reproduces Fig. 7b: expected −52, −55, −54, −75 kJ/mol at the center, 950 µm, bathtub and surface against ≈ −52, −56, −52, −73; the printed form gives +41 kJ/mol at the center, the wrong sign | Mohanty's form (§3.1) |
| α/β switch, Eq. (3) | $h_{\alpha\beta}=1.0-\tfrac12\tanh[(T-T_{\alpha\beta})/2]$, which runs from 0.5 to 1.5 | Figs. 3–6: β inside where $T>662$ °C, α outside | $h_{\alpha\beta}=\tfrac12[1+\tanh((T-935.15)/2)]$ |
| α/β correction term, Eq. (18) | $-29.333333T+24106.6667$ | with it U-22Zr never becomes single-phase γ below 1060 K, contradicting the stated 693 °C transition and Figs. 2–4; Wen Eq. (10) and Jung Eq. (23) print $+24T-23500$, which puts the transition at ≈ 970 K | $+24.0T-23500.0$ |
| $f_\gamma$, Eq. (19) | "$c+G_0^{\beta\text{-Zr}}$"; $4992.02639(1-c)c$; $2609.10711[(1-c)c]^2$ | Quaini's bcc parameters, which these coefficients are, multiply $(1-2c)^k$; with $(1-c)c$ the transition drops to ≈ 942 K | $c\,G^{\beta\text{-Zr}}$, $(1-2c)$, $(1-2c)^2$ |
| Free-energy simplification | chemical mobility divided by "approximately five" | — | 5, in Steps 1.5–1.9 only; Steps 1.3–1.4 use Quaini's unsimplified energies and no factor |
| Modifier meaning, Table 2 | "$M_c^\gamma$" 9.981 and "$\widetilde D_0^\beta$" 2.722 | Fig. 7a and 7b are reproduced only if 9.981 scales $\beta_U,\beta_{\mathrm{Zr}}$ (so both $M_c^\gamma$ and $M_T^\gamma$) and the 1.395 heat-of-transport factor appears alone in $\widetilde Q_\gamma^*$ | $s_\beta$ multiplies $\beta_i$; $s_{D,\beta}$ multiplies $\widetilde D_0^\beta$ in both $M_c^\beta$ and $M_T^\beta$; $s_{Q}$ multiply the heats of transport |
| Linear power | "240 W/m" | Hofman, Jung: 24 kW/m | 24 kW/m (DP-81), 25 kW/m (DP-11) |
| Simulated time | not stated | DP-81 needed 284 EFPD (Hofman, quoted by Jung) | real irradiation time (§3.2); the Part 1 fit slot is reserved for it |
| Fresh conductivity, radial power | cited only (Kim 2014; Galloway 2015) | Kim Eqs. 1, 2, 13, 14 (the leading −102.0 of Eq. 14 and +2992/T of Eq. 2 checked on the PDF image and against Kim Fig. 1); Galloway Fig. 12 fit | as in §3.1 |
| SGTE end-members | cited (Dinsdale 1991) | Jung's Appendix A prints α-U and γ-U, with "ln(T)" in place of $T\ln T$ | Appendix A of this plan (checked against the pure-element transition temperatures) |
| Initial $\phi$, $T$; heat-equation form | not given | Fig. 6a (beginning of life) | §3.2; quasi-steady heat equation |
| PIE data | captured from Hofman 1996 | Figs. 3–5 | digitize from Figs. 4 (DP-81) and 5 (DP-11) |

### 3.4 Steps

#### Step 1.1: Thermodynamics (no PDE)

- Adds: $f_{\alpha\beta}$ and $f_\gamma$ of §3.1, and Quaini's full energies for α-U, β-U and bcc (Quaini Appendix D: α-U $G_{\mathrm{Zr}}=G_{\mathrm{HSER}}^{\mathrm{Zr}}+38000$, ${}^0L=-7000$; β-U $G_U=G^{\mathrm{TET}}_{U}$, $G_{\mathrm{Zr}}=G_{\mathrm{HSER}}^{\mathrm{Zr}}+3400$, ${}^0L=+27000$; bcc ${}^0L=28369.0639-14.741714T$, ${}^1L=4992.02639$, ${}^2L=2609.10711$), as `DerivativeParsedMaterial` objects and in `scripts/lever.py` (common tangent and lever rule).
- Checks: MOOSE and Python agree on $f$, $\partial f/\partial c$ and $\partial^2f/\partial c^2$ at a few $(c,T)$ points (single-element run with `MaterialRealAux` output).
- Key data (Hirschhorn Fig. 2 and §4): for U-22Zr, γ above ≈ 693 °C, β+γ between 662 and 693 °C, α+γ below 662 °C; the simplified and Quaini diagrams agree except in the hatched region ($c<0.15$ above 693 °C). Hirschhorn 2020: the β/γ gap at 750 °C is 2–4 at.% Zr.
- Expected: U-22Zr single-phase γ above ≈ 970 K; $c_\gamma\approx0.67$, 0.57, 0.51 and $\phi\approx0.40$, 0.44, 0.46 at 900, 935, 950 K; $c_{\alpha\beta}\approx0.014$. Quaini at 1023 K: β/γ tie-line 0.0066/0.0225 with ${}^0L_\beta=27000$ and 0.020/0.040 with Hirschhorn 2020's ${}^0L_\beta=21500$.
- Pass: transition temperature within ±10 K, $c_\gamma$ within ±0.03, $\phi$ within ±0.05, the 1023 K tie-line within ±0.005.
- Output: $\phi$, $c_{\alpha\beta}$, $c_\gamma$ as functions of $T$ for the initial conditions of Steps 1.5–1.9.

#### Step 1.2: γ kinetics (no PDE)

- Adds: Li's atomic mobilities (§3.1 table) and the γ thermal mobility in `scripts/mobility.py` and as AD materials.
- Key data and expected values:

| Quantity | Paper | Expected |
| --- | --- | --- |
| Tracer $D^*_U$ in bcc U, 1100 K | Li Fig. 2a | $6.5\times10^{-13}$ m²/s |
| Tracer $D^*_{\mathrm{Zr}}$ in bcc U, 1100 K | Li Fig. 2c | $9.5\times10^{-15}$ m²/s |
| Tracer $D^*_{\mathrm{Zr}}$ in bcc Zr, 1300 K | Li Fig. 2b | $1.8\times10^{-13}$ m²/s |
| $\log_{10}\widetilde D$ at 1223 K, $c=0.2/0.4/0.6/0.8$ (Quaini bcc at local $T$) | Li Fig. 3c,d: ≈ −12.95/−13.3/−13.55/−13.6 | −13.25/−13.56/−13.74/−13.78 |
| $\widetilde D_\gamma=RT M_c^\gamma/(V_mc_\gamma)$ with $s_\beta=9.981$ (Hirschhorn Eq. 20): center, surface | Hirschhorn Fig. 7a: ≈ $2.5\times10^{-14}$, ≈ $8\times10^{-16}$ m²/s | $2.4\times10^{-14}$, $7\times10^{-16}$ m²/s |
| $\widetilde Q_\gamma^*$ with $s_{Q,\gamma}=1.395$ at (974 K, 0.44), (955 K, 0.40), (940 K, 0.55), (900 K, 0.67) | Hirschhorn Fig. 7b: ≈ −52, −56, −52, −73 kJ/mol | −52, −55, −54, −75 kJ/mol |

- Pass: tracer diffusivities within 0.1 decade of Li's curves (they test the mobility parameters alone); $\widetilde Q_\gamma^*$ within ±5 kJ/mol; $\widetilde D_\gamma$ within 0.15 decade. The interdiffusivity sits 0.2–0.3 decades below Li's curve because Li used Chevalier's bcc energy, whose thermodynamic factor is larger; record the offset, which carries into Step 1.3.
- This step fixes the γ thermal-mobility form and its sign before any PDE is run.

#### Step 1.3: Single-phase γ interdiffusion (isothermal couple)

- Adds: the first PDE: split Cahn–Hilliard with `ADSplitCHParsed`, $f_\gamma$ of Quaini's bcc at local $T$ (Hirschhorn's $f_\gamma$ with $T$ in place of $T_{\mathrm{ref}}$), $M_c^\gamma$ with Li's mobilities, no thermotransport.
- Setup: Li Fig. 5b, Zr against U-13 at.% Zr at 1223 K for 120 h; 1-D Cartesian, $-1000\le x\le1000$ µm, 1000 elements; $c=0.999$ for $x<0$ and 0.13 for $x>0$ with a 2 µm tanh step.
- Key data (Muller's points in Li Fig. 5b): $x_U=0.1$ at ≈ −180 µm, $x_U=0.5$ at the Matano plane, $x_U=0.8$ at ≈ +280 µm.
- Expected: ≈ −145, +8 and +220 µm, about 20 % narrower than the data because of the thermodynamic factor (Step 1.2).
- Pass: the two outer positions within ±25 %, the Matano plane within ±10 µm.
- Sanity: $\epsilon_m<10^{-6}$; the Boltzmann–Matano analysis of the simulated profile (`scripts/matano.py`) returns the input $\widetilde D(c)$ of Step 1.2 within 10 %.

#### Step 1.4: Two-phase β/γ KKS diffusion (isothermal couple)

- Adds: KKS with one order parameter and a moving β/γ interface, without a temperature gradient. Thermodynamics: Quaini's β-U (with Hirschhorn 2020's ${}^0L_\beta=21500$) and bcc at local $T$; $V_m=1.309\times10^{-5}$ m³/mol.
- Setup (Hirschhorn 2020 §3–4): β-U against U-18.6 at.% Zr (γ) at 750 °C for 16.5 h; 1-D Cartesian, $-161\le x\le839$ µm, 1000 elements, initial interface at 339 µm (tanh profile of Hirschhorn 2020 Eq. 14 for $c$ from 0 to 0.186 and $\phi$ from 0 to 1). Interface $l=10$ µm, $\sigma=1$ J/m², $\alpha=2.94$, so $\kappa=3l\sigma/\alpha$ and $w=6\alpha\sigma/l$; $L=10^{-9}$ m³/(J s). Mobility $M=(1-h)M_\beta+hM_\gamma$ with the fitted forms of Hirschhorn 2020 Eqs. 15–16, $M_\beta=c_\beta(1-c_\beta)(8\times10^{-19}c_\beta^2+10^{-23})$ and $M_\gamma=c_\gamma(1-c_\gamma)[0.11(0.186-c_\gamma)^{24}+10^{-22}]$ m⁵/(J s).
- Key data (Hirschhorn 2020 Table 1, experiment; the paper's own simulation in parentheses): phase plane $x_P=263$ (260) µm; peak $\lvert\widetilde J/\rho\rvert=40$ (41) $\times10^{-6}$ µm/s; $\widetilde D^{\mathrm{eff}}_\beta=60$ (67) and $\widetilde D^{\mathrm{eff}}_\gamma=30$ (33) $\times10^{-15}$ m²/s; isoconcentration step at the interface from $c\approx0.02$ to $\approx0.04$ (Fig. 2).
- Pass: $x_P$ within ±10 µm; peak flux within ±15 %; $\widetilde D^{\mathrm{eff}}$ within ±25 % (Boltzmann–Matano as in Hirschhorn 2020 Eqs. 17–23, excluding $0.01\le\phi\le0.99$).
- Sanity: KKS constraints; replacing $M_\gamma$ by the Li-based form (Hirschhorn 2020 Eq. 25) reproduces the "DICTRA γ" curve of Hirschhorn 2020 Fig. 7.

#### Step 1.5: Three-phase KKS model in a prescribed temperature field

- Adds: Hirschhorn's α/β/γ model (§3.1) in RZ: simplified energies, $h_{\alpha\beta}(T)$, interface $l=500$ µm, $L=2.5\times10^{-10}$, $M_c$ divided by 5, all modifiers 1; no thermotransport.
- Setup: $T(r)$ prescribed as an AuxVariable from the digitized dashed curve of Hirschhorn Fig. 3 (≈ 707 °C at the center, 627 °C at the surface), held fixed; run for 284 EFPD.
- Key data (Hirschhorn Fig. 3, phase curves): $\phi=1$ for $r\lesssim960$ µm; β+γ from ≈ 960 to ≈ 1640 µm with $\eta_\beta\approx0.5$–0.6; α+γ beyond ≈ 1650 µm with $\eta_\alpha\approx0.56$–0.6.
- Pass: zone boundaries within ±60 µm; $\phi$ within ±0.05.
- Sanity: away from the zone boundaries, $\phi$ equals the lever-rule value of Step 1.1 at the local $T$ within ±0.02; $\epsilon_m<10^{-6}$ with the RZ weight.

#### Step 1.6: Thermotransport in α, β and γ

- Adds: `ADThermalMigration` with $M_T$ of §3.1 (negative in α and β, positive in γ), still under the prescribed temperature of Step 1.5 and with all modifiers 1.
- Key data (Hirschhorn Fig. 3, Zr curve): ≈ 0.26 in the γ zone, a dip to ≈ 0.09 at ≈ 990 µm, ≈ 0.20 in β+γ, ≈ 0.22 in α+γ.
- Pass: the sign of the change in each zone (gain in γ, loss in β+γ, ≈ 0 in α+γ) and zone levels within ±0.02. The prescribed temperature does not include the swelling peak of the real history (Step 1.7), so a closer match is not expected yet.
- Sanity: $M_T=0$ recovers Step 1.5.

#### Step 1.7: Heat conduction

- Adds: $T$ as a nonlinear variable with $k=k_c(B)k_{\mathrm{U\text{-}Zr}}$, the burnup-dependent, composition-dependent heat source of §3.1, and $T(r_0)=900.15$ K; DP-81 conditions of §3.2; modifiers still 1. This step reproduces Hirschhorn's pre-optimization run.
- Key data: Hirschhorn Fig. 3 (end of life): center ≈ 707 °C, Zr levels as in Step 1.6. Hirschhorn Fig. 6, center temperature versus burnup: ≈ 685, 738, 704, 700 °C at 0, 0.75, 1.5, 3.6 at.% (Hirschhorn notes that the temperature hardly changes with the modifiers).
- Expected (composition frozen at 0.22): 686, 741, 709, 707 °C.
- Pass: center temperatures within ±8 K at all four burnups; Fig. 3 Zr levels within ±0.02.
- Sanity: with constant $k$ and uniform $\dot q$, $T(r)=T_s+\dot q(r_0^2-r^2)/(4k)$; the RZ integral of $\dot q$ equals $q'(1-B)$ at all times.

#### Step 1.8: Optimized modifiers (DP-81)

- Adds: Hirschhorn's Table 2 modifiers (§3.2). Nothing else changes.
- Key data (Hirschhorn Fig. 4 at 3.6 at.%): Zr ≈ 0.44 for $r\lesssim850$ µm (PIE ≈ 0.45–0.50); γ edge ≈ 930 µm; β-rich bathtub ($\eta_\beta\approx1$) at ≈ 1100–1350 µm with Zr ≈ 0.01–0.02; Zr ≈ 0.18–0.20 at ≈ 1370–1550 µm, peak ≈ 0.25 at ≈ 1620 µm, ≈ 0.21 in α+γ; α/β switch ≈ 1650 µm; center ≈ 701 °C. Fig. 6: center Zr ≈ 0.30, 0.37, 0.44 at 0.75, 1.5, 3.6 at.%. Fig. 7a/b along the radius as in Step 1.2.
- Pass: zone edges within ±60 µm; γ-zone Zr within ±0.03; bathtub Zr < 0.05; center temperatures within ±8 K. Report the $L_2$ error against the DP-81 PIE next to Hirschhorn's own curve.
- If the key data are missed, the single Part 1 fit is $t_{\mathrm{end}}$, fitted on the γ-zone Zr and then frozen.

#### Step 1.9: DP-11 and Part 1 closure

- Change: DP-11 conditions of §3.2 only; no new physics and no fitting.
- Key data (Hirschhorn Fig. 5): Zr ≈ 0.53 for $r\lesssim950$ µm (PIE ≈ 0.48); γ edge ≈ 1040 µm; bathtub ≈ 0.01 at ≈ 1050–1400 µm; Zr ≈ 0.13–0.20 at ≈ 1400–1700 µm with a peak ≈ 0.30 at ≈ 1750 µm; ≈ 0.20 at ≈ 1800–2100 µm; center ≈ 706 °C.
- Pass: zone edges within ±100 µm; γ-zone Zr within ±0.04; center temperature within ±10 K.
- Closure: convergence on Step 1.8 with 500, 1000 and 2000 elements and the maximum time step × 1, ½, ¼ (burnup kinks always hit); key data change by less than one fifth of their tolerance. Write `p1_hirschhorn/README.md`.

---

## 4. Part 2: Jung et al. (2025), Porosity and Sodium Infiltration

Goal: starting from the Step 1.8 model, replace Hirschhorn's uniform conductivity factor by Jung's porosity- and sodium-dependent conductivity and reproduce Jung's DP-81 parameter study.

### 4.1 Model changes

Jung keeps Hirschhorn's thermodynamics, interface parameters and "optimized mobilities" (Jung §2.1 and Tables 2–3), so everything in §3.1 stays except the conductivity, the heat source and the setup:

$$
k=P_c\,(1-P)^{1.5}\,k_0(T,c),
\qquad
k_0=a_1+a_2T+a_3T^2,
\qquad
a_1=17.5\,\frac{1-2.23W_{\mathrm{Zr}}}{1+1.61W_{\mathrm{Zr}}},
\qquad
a_2=1.54\times10^{-2}\,\frac{1+0.061W_{\mathrm{Zr}}}{1+1.61W_{\mathrm{Zr}}},
\qquad
a_3=9.38\times10^{-6}
$$

$$
P_c=1-\frac{3P_{\mathrm{Na}}}{1-P}\cdot\frac{1-k_{\mathrm{Na}}/k_0}{1.163+1.837\,k_{\mathrm{Na}}/k_0},
\qquad
k_{\mathrm{Na}}=93-0.0581\,\theta+1.173\times10^{-5}\theta^2,
\quad
\theta=T-273.15
$$

$$
P=s_P(B)\left[0.4\,\phi+0.13\,(1-\phi)\right],
\qquad
s_P=\min\left(\frac{B}{0.01},1\right),
\qquad
P_{\mathrm{Na}}=P_{\mathrm{Na}}^{\mathrm{target}}\min\left(\max\left(\frac{B-0.01}{0.005},0\right),1\right),
\qquad
\dot q=\frac{\mathrm{LHGR}\,(1-B)}{\pi r_0^2}
$$

| Item | Value | Source |
| --- | --- | --- |
| DP-81 | LHGR 24 kW/m, $T_s=900$ K, $B_{\mathrm{end}}=0.036$ | Jung §3 |
| DP-11 | LHGR 25 kW/m, $T_s=910$ K, $B_{\mathrm{end}}=0.077$ | Jung §3 |
| Initial state | $c=0.23$; $\phi$ linear from 0.5 (center) to 0.4 (surface); $T$ linear from 990 K to 900 K (initial guess only) | Jung §4, Fig. 3a |
| Sodium level | 0.05 (Baseline), 0.25, 0.50, as fractions | Jung §4 |
| Simulated time | not reported; start from 284 EFPD | fitted once in Step 2.2 if needed |

### 4.2 Paper issues and decisions

| Item | Printed | Evidence | Decision |
| --- | --- | --- | --- |
| α/β switch, Eqs. 3–4 | same inversion as Hirschhorn | Fig. 3a: β inside ($T>935$ K), α outside | as Part 1 |
| $g(\phi)$, Eq. 8 | $\phi^2(1-\phi^2)$, not a double well on $[0,1]$ | Hirschhorn Eq. 7 prints $\phi^2(1-\phi)^2$ | $\phi^2(1-\phi)^2$ |
| Appendix A | "ln(T)" terms; A.3/A.4 set $G^{\alpha\text{-Zr}}=G^{\beta\text{-Zr}}$ equal to SGTE bcc Zr plus 5000 J/mol | SGTE uses $T\ln T$; with the printed $\ln T$ U-22.5Zr never becomes single-phase γ between 900 and 1000 K; the A.3/A.4 variant moves the phase boundaries by less than 0.01 | Appendix A of this plan, as in Part 1 |
| Sign convention, Eq. 1; γ heats of transport | as Hirschhorn; γ values not given | the Fig. 3b bathtub needs opposite signs in γ and β; Jung follows Hirschhorn's "equations and optimized values" | as Part 1 |
| Eq. 31 layout | "3$P_{\mathrm{Na}}$" stacked over "$(1-P)$" | the reading $3P_{\mathrm{Na}}/(1-P)$ gives the ≈ −20 K center shift of Fig. 6b for 50 % versus 5 % sodium; $3P_{\mathrm{Na}}(1-P)$ gives ≈ −13 K (hand estimate) | $3P_{\mathrm{Na}}/(1-P)$, with $P_{\mathrm{Na}}$ the infiltration level as a fraction |
| Eqs. 34–35, $W_{\mathrm{Zr}}(c)$ | circular definition | Jung's own $k$ in Fig. 5c is 15–30 % higher than Eqs. 27–30 with the standard conversion (and than Kim's correlation): $P=0$ center 40 against ≈ 31 W/(m K) | standard conversion; expect the center ≈ 13 K above Jung's, so temperature criteria are stated as ratios |
| Coefficient names | text $a_0,a_1,a_2$, equations $a_1,a_2,a_3$ | — | $a_1,a_2,a_3$ |
| Heat equation | transient form, no $\rho$, $c_p$ | thermal time ≈ 1 s | quasi-steady |
| Geometry | "one-dimensional", $\dot q$ normalized by $\pi r_0^2$ | Jung's own $k(r)$ from Fig. 5c gives the Fig. 5b center temperatures only in RZ (943 and 976 K); a slab gives ≈ 1050 K | RZ |
| Heat source | uniform, no composition factor | differs from Hirschhorn's | as printed by Jung |
| Simulated time | not reported ("considerably shorter" than the experiment) | — | the single Part 2 fit |
| Units of $c$ in the abstract | "0.17 at%" and "0.7 at%" | Fig. 5a center 0.47 → 0.31 ($P=0$); Fig. 6a center 0.47 → 0.405 (50 % sodium) | use the figure values (atom fraction) |
| Surface Zr drop ($r/r_0>0.95$), over-predicted surface α | acknowledged by Jung | — | reproduce, do not fix |

### 4.3 Steps

#### Step 2.1: Jung's heat source with a prescribed conductivity

- Change: Jung's setup (§4.1) and heat source; $k(r)$ prescribed from digitized Jung Fig. 5c for the $P=0$ and Baseline cases; no phase field. $\dot q$ at $B=0.036$ is $1.564\times10^{9}$ W/m³.
- Key data (Jung Fig. 5b): center $T_c\approx943$ K ($P=0$) and ≈ 977 K (Baseline).
- Expected: 943 K and 976 K; a slab would give ≈ 1050 K for the Baseline.
- Pass: within ±3 K.
- Sanity: the RZ integral of $\dot q$ equals $\mathrm{LHGR}(1-B)=23.14$ kW/m.

#### Step 2.2: Jung's $k_0(T,c)$ without porosity

- Adds: the coupled model of Step 1.8 with $k=k_0(T,c)$ in place of $k_c(B)k_{\mathrm{U\text{-}Zr}}$ ($P=0$, $P_{\mathrm{Na}}=0$; Jung's $P=0$ case keeps 5 % sodium, worth ≈ 1 K).
- Fit: run with $t_{\mathrm{end}}=284$ EFPD first; if the center Zr misses the key data, fit $t_{\mathrm{end}}$ once on it and freeze it for Steps 2.3–2.5.
- Key data (Jung Fig. 5, $P=0$): $T_c\approx943$ K; center Zr ≈ 0.31, falling to ≈ 0.25 near $r/r_0\approx0.25$; no Zr bathtub and no single-γ zone; $k\approx40$–42 W/(m K).
- Pass: center Zr within ±0.05; no bathtub. Report $T_c$ (≈ 956 K expected, §4.2) and the ratio of $k(r)$ to Fig. 5c (≈ 0.75–0.8 expected).

#### Step 2.3: Porosity from phase fractions

- Adds: $(1-P)^{1.5}$ in $k$. First constant $P=0.1$, 0.2, 0.3 from $t=0$ (no ramp), then the phase-dependent $P$ with the ramp $s_P(B)$, which is the Baseline without sodium ($P_{\mathrm{Na}}=0$; Jung's 5 % is worth 2–3 K).
- Key data (Jung Fig. 5). The temperature criterion is the center rise $\Delta T_c=T_c-900$ K divided by the $P=0$ rise of Step 2.2, which is insensitive to the conductivity offset of §4.2.

| Case | $T_c$ (K) | $\Delta T_c/\Delta T_c(P=0)$ | Center Zr | γ-only zone at the center |
| --- | --- | --- | --- | --- |
| $P=0.1$ | ≈ 951 | 1.19 | ≈ 0.36 | no |
| $P=0.2$ | ≈ 960 | 1.40 | ≈ 0.43 | no ($\phi\approx0.71$) |
| $P=0.3$ | ≈ 973 | 1.70 | ≈ 0.47 | yes |
| $P(\eta)$, Baseline without sodium | ≈ 977 | 1.79 | ≈ 0.47 | yes |

- Pass: ratios within ±0.08; center Zr within ±0.05; presence or absence of the γ-only zone correct in every case.
- Sanity: after the ramp $P$ stays in $[0.13,0.4]$; the kink at $B=0.01$ is resolved.

#### Step 2.4: Sodium infiltration (full Jung Baseline)

- Adds: $P_c$ (Jung Eqs. 31–32) with the sodium ramp between $B=0.01$ and 0.015, for $P_{\mathrm{Na}}^{\mathrm{target}}=0.05$ (Baseline), 0.25 and 0.50.
- Key data (Jung Fig. 6), with the same temperature ratio as Step 2.3 relative to the 5 % case:

| Sodium | $T_c$ (K) | $\Delta T_c/\Delta T_c(5\,\%)$ | Center Zr | Center $\phi$ | Center $k$ (W/(m K)) |
| --- | --- | --- | --- | --- | --- |
| 5 % | ≈ 977 | 1 | ≈ 0.47 | 1.0 | ≈ 15 |
| 25 % | ≈ 962 | 0.81 | ≈ 0.47 | ≈ 0.82 | ≈ 21 |
| 50 % | ≈ 956 | 0.73 | ≈ 0.405 | ≈ 0.65 | ≈ 27 |

- Baseline against Jung Fig. 3b: center Zr ≈ 0.47 rising to ≈ 0.51 at $r/r_0\approx0.38$; Zone I/II boundary $r/r_0\approx0.40$; bathtub Zr ≈ 0.02 at $r/r_0\approx0.53$–0.63; Zone II/III boundary ≈ 0.64; Zone III Zr ≈ 0.22; surface drop beyond $r/r_0\approx0.95$. Report the $L_2$ error against the DP-81 PIE (Chang 1994, as plotted by Jung).
- Pass: temperature ratios within ±0.05; center-Zr drop from 5 % to 50 % equal to 0.065 ± 0.03; center $\phi$ within ±0.1; Baseline zone boundaries within ±0.05 in $r/r_0$.
- Sanity: $P_c=1$ when $P_{\mathrm{Na}}=0$, and $P_c>1$ everywhere because $k_{\mathrm{Na}}>k_0$ ($k_{\mathrm{Na}}/k_0\approx2$ at 900 K). This step also confirms the Eq. 31 reading of §4.2.

#### Step 2.5: DP-11 and Part 2 closure

- Change: DP-11 conditions of §4.1; $t_{\mathrm{end}}$ scaled by the burnup ratio at the same burnup rate per unit power (an assumption); Jung gives no separate DP-11 initial state, so the DP-81 one is used.
- Key data (Jung Fig. 4): Zone I Zr ≈ 0.50–0.53 for $r/r_0<0.5$; zone boundaries at $r/r_0\approx0.52$ and ≈ 0.72; bathtub Zr ≈ 0.02 at $r/r_0\approx0.58$–0.70; Zone III Zr ≈ 0.20 at $r/r_0\approx0.73$–0.88; $\eta_\alpha\to1$ at the surface. PIE: Zone I ≈ 0.48, Zone III ≈ 0.20.
- Pass: zone boundaries within ±0.05 in $r/r_0$; Zone I Zr within ±0.05.
- Closure: convergence on the Baseline as in Step 1.9; write `p2_jung/README.md`.

---

## 5. Part 3: Wen et al. (2022), Irradiation-Enhanced Diffusion

Goal: test Wen's central claim on the common model. Hirschhorn and Jung enhance the γ and β mobilities by fitted factors of 9.981 and 2.722. Wen proposes that irradiation-produced vacancies supply this enhancement physically.

### 5.1 Model change

Starting from the Step 2.4 Baseline, the two mobility modifiers are set back to 1 ($s_\beta=s_{D,\beta}=1$; the heat-of-transport modifiers stay), and the chemical mobility is multiplied by Wen's factor (Wen Eqs. 21, 23, 24):

$$
M_c^{I}=\xi(T)\,M_c,
\qquad
\xi(T)=1+\frac{C_v^r}{C_v^e(T)},
\qquad
C_v^e=A_v\exp\left(-\frac{E_v^f}{k_BT}\right),
\qquad
E_v^f=1.20\ \mathrm{eV},\quad A_v=1,\quad C_v^r=7.04\times10^{-7}.
$$

The run uses the real irradiation time (284 EFPD for DP-81), not Part 2's fitted time, because the purpose is to see whether $\xi$ explains the kinetics without a fitted factor.

### 5.2 Paper issues and decisions

| Item | Printed | Note | Decision |
| --- | --- | --- | --- |
| What $\xi$ multiplies | $M_c$ only (Eq. 24) | a vacancy mechanism would scale the atomic mobilities and hence $M_T$ too; applying $\xi$ to $M_c$ alone lowers the ratio $M_T/M_c$ by $\xi$ | as printed |
| $C_v^r$ | $7.04\times10^{-7}$, computed for 30 kW/m and ≈ 1000 K | DP-81 ran at 24 kW/m, DP-11 at 25 kW/m | constant, as printed |
| Wen's other ingredients | α/β mobilities without $c_{\alpha\beta}$, unreported γ prefactor, all-positive heats of transport, $L(T,c_\gamma)$, $\kappa_c\sim10^{-9}$, 1.3×10⁵ s explicit run | Hirschhorn 2019 and 2020 show the results are insensitive to $L$ once the process is diffusion-controlled | not adopted: Part 3 changes only the enhancement |

### 5.3 Steps

#### Step 3.1: $\xi(T)$ in place of the fitted enhancement (DP-81)

- Adds: §5.1. Nothing else changes.
- Key data: the DP-81 PIE (as in Step 2.4); Wen Fig. 3 as the paper's own result with irradiation enhancement (γ zone Zr ≈ 0.47–0.49 up to ≈ 800 µm, Zr well ≈ 0.07 at ≈ 1100–1450 µm, γ+α zone ≈ 0.22–0.29); the Step 2.4 Baseline as the fitted-enhancement reference.
- Expected: $\xi\approx4.69$ at 900 K, 3.07 at 935 K and 1.93 at 988 K, so the enhancement in the γ zone (≈ 2) is well below the fitted 9.981 and in the β zone (≈ 3) close to the fitted 2.722. Hirschhorn's sensitivity study (Fig. 10) shows that γ mobility at the low end of its range gives an almost unredistributed profile, so a weaker Zone I is expected.
- Outcome criterion: Wen's claim is supported if the $L_2$ error against the PIE is at most 1.2 times that of the Step 2.4 Baseline and the zone boundaries are within ±0.05 in $r/r_0$ of Step 2.4. Otherwise record that $\xi$ alone does not reproduce DP-81, with the γ-zone Zr and bathtub depth that it gives.
- Sanity: $C_v^r=0$ gives $\xi=1$ and reproduces the unenhanced model (the analogue of Wen Fig. 2 and Hirschhorn Fig. 3).

#### Step 3.2: DP-11 with $\xi(T)$ and Part 3 closure

- Change: DP-11 conditions (real time ≈ 583 EFPD); no new physics.
- Key data: DP-11 PIE (as in Step 2.5); Step 2.5 as the fitted-enhancement reference.
- Outcome criterion: as Step 3.1. Write `p3_wen/README.md` with the conclusion on Wen's claim for both pins.

---

## 6. Data to Digitize

Store one CSV per curve in `problem/uzr/data/`, with a header recording the paper, figure, axis units and an estimated reading error (WebPlotDigitizer or similar).

| Paper | Figure or table | Curves | Used in |
| --- | --- | --- | --- |
| Hirschhorn 2019 | Fig. 2 | phase boundaries (simplified and Quaini) | 1.1 |
| Li 2010 | Figs. 2, 3c–d, 5b | tracer diffusivities; $\widetilde D(c)$ at 1223 K; Muller's points and Li's curve | 1.2, 1.3 |
| Hirschhorn 2020 | Fig. 1 or 2, Fig. 7; Table 1 | experimental $c(x)$; DICTRA-γ curve; response values (no digitizing) | 1.4 |
| Hirschhorn 2019 | Fig. 3 | $T(r)$, $\phi$, $\eta_\alpha$, $\eta_\beta$, Zr | 1.5–1.7 |
| Hirschhorn 2019 | Figs. 4, 6, 7 | DP-81 final and intermediate profiles; center $T$ and Zr per burnup; $\widetilde D_\gamma$, $\widetilde Q_\gamma^*$ | 1.2, 1.7, 1.8 |
| Hirschhorn 2019 | Figs. 4, 5 | DP-81 and DP-11 PIE points | 1.8, 1.9, 3.1, 3.2 |
| Jung 2025 | Figs. 3b, 5a–d, 6a–d, 4 | $c$, $T$, $k$, $\phi$ for all cases; DP-81 and DP-11 PIE points | 2.1–2.5, 3.1 |
| Wen 2022 | Fig. 3 | $c(r)$, $\phi(r)$ with irradiation enhancement | 3.1 |

`scripts/compare.py` reads an Exodus file and the matching CSV files, samples the simulation at the key points of the step and prints the key-data table (paper value, simulation value, tolerance, pass/fail) for the step `README.md`. `scripts/matano.py` performs the Boltzmann–Matano analysis of Steps 1.3–1.4.

---

## 7. References and Files

Papers (in `palm/doc/`):

| File | Reference | Role |
| --- | --- | --- |
| `paper_Hirschhorn_2019` | Hirschhorn, Tonks, Aitkaliyeva, Adkins, J. Nucl. Mater. 523 (2019) 143–156 | Part 1 |
| `paper_Li_2010` | Li et al., J. Nucl. Mater. 407 (2010) 220–227 | γ atomic mobilities; Step 1.3 |
| `paper_Hirschhorn_2020` | Hirschhorn et al., J. Nucl. Mater. 529 (2020) 151929 | Step 1.4 |
| `paper_Quaini_2018` | Quaini et al., J. Nucl. Mater. 501 (2018) 104–131 | full U–Zr CALPHAD energies (Steps 1.1, 1.3, 1.4) |
| `paper_Kim_2014` | Kim, Cho, Sohn, J. Nucl. Mater. 445 (2014) 272–280 | fresh-fuel conductivity (Step 1.7) |
| `paper_Galloway_2015` | Galloway et al., Nucl. Eng. Des. 286 (2015) 1–17 | power profile; sign convention of the heats of transport |
| `paper_Jung_2025` | Jung, Kim, Chang, Ann. Nucl. Energy 224 (2025) 111679 | Part 2 |
| `paper_Wen_2022` | Wen, Liu, Yun, Qian, Front. Mater. 9 (2022) 950785 | Part 3 |
| `paper_Mohanty_2011` | Mohanty, Bush, Okuniewski, Sohn, J. Nucl. Mater. 414 (2011) 211–216 | source of the γ thermal-mobility form and heats of transport; not reproduced |

Not in `palm/doc/`: Dinsdale, Calphad 15 (1991) 317 (SGTE unary data; the functions used are listed in Appendix A) and Hofman et al., J. Nucl. Mater. 227 (1996) 277 (DP-81 and DP-11 PIE, cited through Hirschhorn and Jung).

Code:

- `palm/include/kernels/ADThermalMigration.h`, `palm/src/kernels/ADThermalMigration.C`, `palm/test/tests/kernels/ad_thermal_migration/`
- `palm/problem/uzr/p1_mohanty/step1_1_kernel_swap/` and `palm/problem/uzr/s1_mohanty/` (archived)
- MOOSE references: `moose/modules/phase_field/test/tests/KKS_system/kks_example_split.i` (KKS split form), `moose/modules/heat_transfer/test/tests/code_verification/cylindrical_test_no1.i` (1-D RZ heat conduction)

Not covered by this plan: two-dimensional $r$–$z$ geometry, U–Pu–Zr, and any physics not contained in the papers.

---

## Appendix A: SGTE Unary Functions Used (J/mol)

| Function | Phase | Range (K) | Expression |
| --- | --- | --- | --- |
| $G^{\alpha\text{-U}}$ (GHSERUU) | orthorhombic α-U | 298.15–955 | $-8407.734+130.955151T-26.9182\,T\ln T+1.25156\times10^{-3}T^2-4.42605\times10^{-6}T^3+38568/T$ |
| | | 955–3000 | $-22521.8+292.121093T-48.66\,T\ln T$ |
| $G^{\beta\text{-U}}$ (GTETUU) | tetragonal β-U | 298.15–941.5 | $-5156.136+106.976316T-22.841\,T\ln T-1.084475\times10^{-2}T^2+2.7889\times10^{-8}T^3+81944/T$ |
| | | 941.5–3000 | $-14327.309+244.16802T-42.9278\,T\ln T$ |
| $G^{\gamma\text{-U}}$ (GBCCUU) | bcc γ-U | 298.15–1049 | $-752.767+131.5381T-27.5152\,T\ln T-8.35595\times10^{-3}T^2+9.67907\times10^{-7}T^3+204611/T$ |
| | | 1049–3000 | $-4698.365+202.685635T-38.2836\,T\ln T$ |
| $G^{\alpha\text{-Zr}}$ (GHSERZR) | hcp α-Zr | 130–2128 | $-7827.595+125.64905T-24.1618\,T\ln T-4.37791\times10^{-3}T^2+34971/T$ |
| $G^{\beta\text{-Zr}}$ (GBCCZR) | bcc β-Zr | 130–2128 | $-525.539+124.9457T-25.607406\,T\ln T-3.40084\times10^{-4}T^2-9.729\times10^{-9}T^3+25233/T-7.6143\times10^{-11}T^4$ |

Checks that a transcription must pass, and that these expressions pass: $G^{\beta\text{-U}}=G^{\alpha\text{-U}}$ at 941.5 K and $G^{\beta\text{-U}}=G^{\gamma\text{-U}}$ at 1049 K (pure-U transitions); $G^{\beta\text{-Zr}}=G^{\alpha\text{-Zr}}$ at 1139 K (pure-Zr transition); each two-range function is continuous at its break. Within 900–1020 K the two α-U ranges differ by less than 1 J/mol, so the low range printed in Jung's appendix can be used throughout.

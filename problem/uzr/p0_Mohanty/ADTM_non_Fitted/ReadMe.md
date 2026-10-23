# Fig. 6 with `ADThermalMigration`, no fitted parameters

Single-phase $\gamma$(bcc) U–Zr thermotransport on the 300 $\mu\mathrm{m}$ bar of Mohanty et al., J. Nucl. Mater. 414 (2011) 211–216, Fig. 6.

Nothing is adjusted to bring the profile closer to Fig. 6. `ADThermalMigration` reads $M_T$ directly, and the split Cahn–Hilliard kernels are the automatic-differentiation versions. `SoretDiffusion` is not used.

Both atomic mobilities depend on temperature through Eq. (4). That is the combination in Fig. 5c, whose 30-day curve has the same shape as Fig. 6. Figs. 5a and 5b hold $\beta_{\mathrm{Zr}}$ constant and are not this calculation. $\widetilde Q^*$ is the constant average from the text.

Input: `step1_1.i`. Papers: `palm/doc/Mohanty_2011`, `palm/doc/Jung_2025`.

## Equation decomposition

Mohanty Eq. (1), with the gradient-energy coefficient set to zero:

$$
\frac{\partial c}{\partial t}
= \nabla\cdot\left(M\nabla w\right)
- \nabla\cdot\left(M_T\frac{\nabla T}{T}\right),
\qquad
w=\frac{\partial f}{\partial c}.
$$

Eqs. (2)–(4):

$$
M=c(1-c)\left[c\,\beta_U+(1-c)\beta_{\mathrm{Zr}}\right],
\qquad
M_T=c(1-c)\left[\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right],
$$

$$
\beta_i=\beta_{0,i}\exp\left(-\frac{H_i}{RT}\right).
$$

$M_T>0$ sends Zr toward the hot end. The bracket in $M_T$ has no extra factor of $c$ or $(1-c)$.

Integration by parts, with zero boundary flux, gives the weak form actually solved. $\psi$ is the test function of the $w$-equation:

$$
\left(\psi,\frac{\partial c}{\partial t}\right)
+\left(\nabla\psi,\,M\nabla w\right)
-\left(\nabla\psi,\,\frac{M_T}{T}\nabla T\right)
=0,
$$

$$
\left(\psi_c,\frac{\partial f}{\partial c}-w\right)
+\left(\nabla\psi_c,\,\kappa_c\nabla c\right)
=0.
$$

$\kappa_c=0$, so the second line is $w=\partial f/\partial c$.

## Kernel syntax

Each line of the weak form is one kernel in `step1_1.i`. All four are AD kernels, and every material they read is an AD material, so `solve_type = NEWTON` differentiates $M(c,T)$, $M_T(c,T)$ and $\partial f/\partial c$.

| Weak-form term | Block | Type | Arguments |
| --- | --- | --- | --- |
| $(\psi,\partial c/\partial t)$ | `c_dot` | `ADCoupledTimeDerivative` | `variable = w`, `v = c` |
| $(\nabla\psi, M\nabla w)$ | `w_res` | `ADSplitCHWRes` | `variable = w`, `mob_name = M` |
| $-(\nabla\psi, M_T\nabla T/T)$ | `w_res_tm` | `ADThermalMigration` | `variable = w`, `T = T`, `thermal_mobility = M_T` |
| $w=\partial f/\partial c$ | `c_res` | `ADSplitCHParsed` | `variable = c`, `f_name = f_bulk`, `kappa_name = kappa_c`, `w = w`, `coupled_variables = T` |

`ADThermalMigration` implements

$$
R=-\int_\Omega\frac{M_T}{T}\,\nabla T\cdot\nabla\psi\,\mathrm dV.
$$

The kernel contains no $k_B$, no diffusivity and no heat of transport. The sign and the units sit entirely in the property `M_T`. That property is Eq. (3).

`variable = w` on the first three kernels means their test function is the shape function of $w$. `v = c` means the time derivative is taken of $c$. `coupled_variables = T` tells `ADSplitCHParsed` that $f$ depends on the auxiliary temperature; $T$ is not a nonlinear unknown.

Temperature is an `AuxVariable` set by `FunctionAux` on `INITIAL TIMESTEP_BEGIN`:

```text
T = 1405 - (1405 - 1050)/300 * x
```

## Parameters taken from Mohanty as printed

| Quantity | Value | Where it is set | Source |
| --- | --- | --- | --- |
| Length | 300 $\mu\mathrm{m}$ | `Mesh/xmax` | §3 |
| Elements | 300 linear, $\Delta x=1\,\mu\mathrm{m}$ | `Mesh/nx` | not given in the paper; see below |
| Initial $c$ | 0.39 | `ICs/c_IC` | Fig. 6, U–39 at.% Zr |
| Duration | $2.592\times10^{6}\,\mathrm{s}$ = 30 day | `Executioner/end_time` | Figs. 6 and 7 |
| Composition BC | no `[BCs]` block | natural BC of the weak form | Eq. (6); see below |
| $H_U$ | 128000 J/mol | `beta_U` | text, Hofman 1996 |
| $H_{\mathrm{Zr}}$ | 195000 J/mol | `beta_Zr` | text, Hofman 1996 |
| $\widetilde Q_U^*$ | $+2.5\times10^{-20}\,\mathrm{J}$ = $+0.156037728\,\mathrm{eV}$ | `M_T`, constant `QU` | text; average of Fig. 3 |
| $\widetilde Q_{\mathrm{Zr}}^*$ | $-17.5\times10^{-20}\,\mathrm{J}$ = $-1.092264098\,\mathrm{eV}$ | `M_T`, constant `QZ` | text; average of Fig. 3 |
| $\kappa_c$ | 0 | `constants` | “vanishes in a single-phase alloy” |
| Hot end | $x=0$, 1405 K | `T_profile` | Figs. 5 and 6 put the hot end on the left; the text states 1405 K |
| Cold end | $x=300\,\mu\mathrm{m}$, 1050 K | `T_profile` | same figures; text states 1050 K |
| $R$ | 8.314462618 J/(mol·K) | `beta_U`, `beta_Zr` | gas constant in Eq. (4) |

Eq. (8) prints $T(0)=T_{\min}$ and $T(L)=T_{\max}$. Figs. 5 and 6 draw the opposite assignment, and those are the figures being reproduced, so the hot end is at $x=0$. The end values are the stated 1405 K and 1050 K. The 1400 K mark on the figure axis is a tick, not a second temperature.

### Why the mesh has 300 elements

Section 3 states the domain length, $L=300\,\mu\mathrm{m}$, and that Eq. (5) was solved by a control-volume method. It does not state how many control volumes were used. `nx = 300` is a numerical choice: one linear element per micrometre, so $\Delta x=1\,\mu\mathrm{m}$. It is not a parameter of the paper.

### Why there is no boundary-condition block

Eq. (6) sets $\partial c/\partial x=0$ and $\partial^3 c/\partial x^3=0$ at $x=0$ and $x=L$.

The kernels integrate the divergences by parts and do not add a boundary integral. The resulting natural condition is zero total flux,

$$
M\frac{\partial w}{\partial x}-M_T\frac{1}{T}\frac{\partial T}{\partial x}=0,
$$

at both ends. That is the closed-system condition. Writing a `NeumannBC` of value 0 would impose the same flux a second time, so it is not added.

$\partial^3 c/\partial x^3=0$ is the natural condition that belongs to the gradient-energy term $\kappa_c\nabla^2 c$. The paper says $\kappa_c$ vanishes in a single-phase alloy, and the input sets it to 0, so that boundary condition has nothing left to constrain.

Temperature is assigned by `FunctionAux`. It is not obtained by solving Eq. (8). Eq. (8) as printed asks for zero heat flux and for both end temperatures at once, which a constant-conductivity bar cannot satisfy unless the two end temperatures are equal. The linear field used here is the one the text says those conditions were meant to produce.

The heats are joule per atom. Division by $1.60217662\times10^{-19}\,\mathrm{J/eV}$ converts them to eV per atom. Division by 96485 would treat $10^{-20}\,\mathrm{J}$ as J/mol, which the text does not say. $\widetilde Q_{\mathrm{Zr}}^*$ is stored with its printed negative sign, and the minus in Eq. (3) is written once.

## Quantities Mohanty does not tabulate

$\beta_{0,i}$ is not tabulated. The values below were fitted to Fig. 2 and are used directly. They are not adjusted to Fig. 6. Eq. (4) is applied to both U and Zr:

$$
\beta_{0,\mathrm{U}}=119945.49892\times10^{8}=1.1994549892\times10^{13},
$$

$$
\beta_{0,\mathrm{Zr}}=5194864.7\times10^{8}=5.1948647\times10^{14}.
$$

### Unit printed on Fig. 2

Fig. 2 labels the axis $(10^{8}\,\mathrm{m}^{2}\cdot\mathrm{sec}/\mathrm{J})$. Read as written, that is square metres multiplied by time and divided by joules, $\mathrm{m}^{2}\,\mathrm{s}/\mathrm{J}$.

That label does not close the other statements in the same paper:

- Fig. 4 plots $\beta_i\widetilde Q_i^*$ in $\mathrm{m}^{2}/\mathrm{s}$. With $\widetilde Q^*$ in joules, the product has those units only if $\beta$ is $\mathrm{m}^{2}/(\mathrm{J\,s})$, not $\mathrm{m}^{2}\,\mathrm{s}/\mathrm{J}$. The two labels differ by $\mathrm{s}^{2}$.
- In Eq. (1), $\beta\nabla(\partial f/\partial c)$ has to be a velocity when $\partial f/\partial c$ is an energy per atom. That requires $\beta$ in $\mathrm{m}^{2}/(\mathrm{J\,s})$.
- The reduced time defined in the paper is proportional to $\beta t/L^{2}$. That ratio is dimensionless only for $\mathrm{m}^{2}/(\mathrm{J\,s})$.

The input therefore keeps the Fig. 2 magnitudes and assigns them the dimension $\mathrm{m}^{2}/(\mathrm{J\,s})$. At 1400 K this produces $\beta_U\widetilde Q_U^*=5.03\times10^{-12}\,\mathrm{m}^{2}/\mathrm{s}$ and $\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*=-4.82\times10^{-12}\,\mathrm{m}^{2}/\mathrm{s}$, which sit on the scale of Fig. 4. Treating the Fig. 2 label literally as $\mathrm{m}^{2}\,\mathrm{s}/\mathrm{J}$ would make $\beta\widetilde Q^*$ come out in $\mathrm{m}^{2}\,\mathrm{s}$ and would make $\partial c/\partial t$ not have dimension $1/\mathrm{s}$.

Conversion into the input, where energy is eV per atom and length is $\mu\mathrm{m}$:

$$
1\,\mathrm{m}^{2}/(\mathrm{J\,s})=1.60217662\times10^{-7}\,\mu\mathrm{m}^{2}/(\mathrm{eV\,s}).
$$

| | $\beta_{0,i}$ from Fig. 2 | value stored in `step1_1.i` |
| --- | --- | --- |
| U | $1.1994549892\times10^{13}\,\mathrm{m}^{2}/(\mathrm{J\,s})$ | $1.92173874\times10^{6}\,\mu\mathrm{m}^{2}/(\mathrm{eV\,s})$ |
| Zr | $5.1948647\times10^{14}\,\mathrm{m}^{2}/(\mathrm{J\,s})$ | $8.32309077\times10^{7}\,\mu\mathrm{m}^{2}/(\mathrm{eV\,s})$ |

Mohanty does not print $f$. The $\gamma$ molar Gibbs energy is Jung et al. (2025) Eq. (24), with $T_{\mathrm{ref}}=1010\,\mathrm{K}$:

$$
L_0=28369.0639-14.741714\,T_{\mathrm{ref}},\quad
L_1=4992.02639,\quad
L_2=2609.10711.
$$

End-members are Appendix A. $G_0^{\beta\text{-Zr}}=G_0^{\alpha\text{-Zr}}$ as (A.4). The numerical coefficients in `free_energy` are those appendix entries:

| Symbol in `step1_1.i` | Appendix term |
| --- | --- |
| `aU` … `fU` | $G_0^{\gamma\text{-U}}$: $-752.767$, $131.5381\,T$, $-27.5152\,T\ln T$, $-8.35595\times10^{-3}T^2$, $0.967907\times10^{-6}T^3$, $204611/T$ |
| `aZ` … `gZ` | $G_0^{\alpha\text{-Zr}}$: $4474.461$, $124.9457\,T$, $-25.607406\,T\ln T$, $-3.40084\times10^{-4}T^2$, $-9.792\times10^{-9}T^3$, $25233/T$, $-7.6143\times10^{-11}T^4$ |

Appendix A prints $\ln T$. The numbers are Dinsdale’s coefficients of $T\ln T$, and the appendix cites that SGTE source, so the input evaluates $T\ln T$. Jung prints one $G_0^{\gamma\text{-U}}$ polynomial, the branch valid below 1049 K. That polynomial is used on the whole 1050–1405 K bar, because no other branch is given.

$G$ is in J/mol. Division by `JeV = 96485.3321` puts $f$ in eV/atom, which matches $\beta$ in $\mu\mathrm{m}^2/(\mathrm{eV\,s})$. $R/F=k_B$, so the ideal-mixing term is `kB*Tref*(c*log(c)+(1-c)*log(1-c))`. In the parsed expression, `log` is the natural logarithm.

## What this input does not contain

- No multiplier on $\beta_U$ or $\beta_{\mathrm{Zr}}$.
- No composition weight inside the bracket of Eq. (3).
- No `D_soret` and no `Qeff`. Those two properties exist only to feed `SoretDiffusion`.
- $\kappa_c$ is zero, not $10^{-5}$.

## Run

```bash
conda activate moose
cd /home/mzw/opt/moose_projects/palm/problem/uzr/p1_mohanty/step1_1_kernel_swap
../../../../palm-opt -i step1_1.i
conda activate analysis
python plot_exodus_cx_profile.py step1_1_out.e
```

`Executioner` is BDF2, Newton with LU, and `IterationAdaptiveDT` from $\mathrm{d}t=1\,\mathrm{s}$ up to $\mathrm{d}t_{\max}=3600\,\mathrm{s}$. Output times are 5, 10, 20 and 30 days.

## Result

The run finished in 114 s. The element average stays at $c_{\mathrm{avg}}=0.390000$.

| $t$ (day) | $c(0)$, hot | $c(300\,\mu\mathrm{m})$, cold |
| --- | --- | --- |
| 0.2 | 0.507 | 0.333 |
| 1 | 0.550 | 0.144 |
| 5 | 0.585 | 0.025 |
| 10 | 0.591 | 0.002 |
| 30 | 0.591 | 0.002 |

Fig. 6 at 30 days is still evolving: about $+11$ at.% at the hot end ($c\approx0.50$) and about $-7$ at.% at the cold end ($c\approx0.32$), with a shoulder near the cold end. This run is already close to steady state by 10 days, and the steady segregation is stronger. That is the result of the printed mobilities and Eq. (3). No factor was added to slow the kinetics or to reduce $M_T$.

![c(x)](step1_1.png)

# Mohanty et al. (2011) Fig. 6, paper parameters

Single-phase $\gamma$(bcc) U–Zr thermotransport on a 300 $\mu\mathrm{m}$ bar. The input uses `SoretDiffusion`. No parameter is adjusted to improve the match with Fig. 6.

Paper: `palm/doc/Mohanty_2011`. Free energy: `palm/doc/Jung_2025`.

## What is taken from Mohanty as printed

$$
\frac{\partial c}{\partial t}
= \nabla\cdot\left(M\nabla w\right)
- \nabla\cdot\left(M_T\frac{\nabla T}{T}\right),
\qquad
w=\frac{\partial f}{\partial c}
$$

$$
M=c(1-c)\left[c\,\beta_U+(1-c)\beta_{\mathrm{Zr}}\right],
\qquad
M_T=c(1-c)\left[\beta_U\widetilde Q_U^*-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right]
$$

$$
\beta_i=\beta_{0,i}\exp(-H_i/RT)
$$

| Quantity | Value | Source |
| --- | --- | --- |
| Length, mesh | 300 $\mu\mathrm{m}$, 300 elements | §3 |
| Initial $c$ | 0.39 | Fig. 6 |
| Time | 30 day | Figs. 6 and 7 |
| Composition BC | zero flux | Eq. (6) |
| $H_U$, $H_{\mathrm{Zr}}$ | 128, 195 kJ/mol | text, Hofman 1996 |
| $\widetilde Q_U^*$, $\widetilde Q_{\mathrm{Zr}}^*$ | $+2.5\times10^{-20}$, $-17.5\times10^{-20}$ J | text; heat per atom; Fig. 3 |
| $\kappa_c$ | 0 | “vanishes in a single-phase alloy” |
| Hot / cold | $x=0$ at 1405 K, $x=300\,\mu\mathrm{m}$ at 1050 K | Figs. 5–6 put the hot end on the left; the text states 1405 K and 1050 K |

Eq. (8) prints $T(0)=T_{\min}$ and $T(L)=T_{\max}$. Figs. 5 and 6 do the opposite and are the figure being reproduced, so the hot end is at $x=0$. The end values are the stated 1405 K and 1050 K, not the 1400 K axis tick.

$M_T$ is written as Eq. (3). `SoretDiffusion` does not accept $M_T$; its residual is $\int D\,Q\,c/(k_B T^2)\,\nabla T\cdot\nabla\psi\,\mathrm dV$ with $k_B$ in eV/K. Setting $D=M k_B T/c$ and $Q=-M_T/M$ makes that residual equal to the printed flux. This $Q$ is the kernel argument, not $\widetilde Q^*$. Because $\beta$ is stored in $\mu\mathrm{m}^2/(\mathrm{eV\,s})$, the joule values are divided by $1.60217662\times10^{-19}\,\mathrm{J/eV}$ and kept as eV per atom. Dividing by 96485 would treat $10^{-20}\,\mathrm{J}$ as J/mol.

## Quantities the paper does not tabulate

$\beta_{0,i}$ is only Eq. (4) plus Figs. 2 and 4. Read at 1400 K, with the printed $\widetilde Q^*$ so that $\beta_i\widetilde Q_i^*$ lies on Fig. 4:

| | $\beta_i(1400\,\mathrm{K})$ | $\beta_{0,i}$ |
| --- | --- | --- |
| U | $1.98\times10^{8}\,\mathrm{m}^2/(\mathrm{J\,s})$ | $1.892434\times10^{6}\,\mu\mathrm{m}^2/(\mathrm{eV\,s})$ |
| Zr | $2.65\times10^{7}\,\mathrm{m}^2/(\mathrm{J\,s})$ | $8.004853\times10^{7}\,\mu\mathrm{m}^2/(\mathrm{eV\,s})$ |

Reading error is about 10%. These prefactors are not adjusted to Fig. 6.

Mohanty does not print $f$. The $\gamma$ molar Gibbs energy is Jung et al. (2025) Eq. (24), with $T_{\mathrm{ref}}=1010\,\mathrm{K}$ and

$$
L_0=28369.0639-14.741714\,T_{\mathrm{ref}},\quad
L_1=4992.02639,\quad
L_2=2609.10711.
$$

End-members are Appendix A. $G_0^{\beta\text{-Zr}}=G_0^{\alpha\text{-Zr}}$ as (A.4). Jung prints one $G_0^{\gamma\text{-U}}$ polynomial, Dinsdale’s branch below 1049 K; that polynomial is used on the whole 1050–1405 K bar.

Appendix A prints $\ln T$. The numbers are Dinsdale’s coefficients of $T\ln T$, and the appendix cites that source, so the input evaluates $T\ln T$. Energies are divided by 96485.3321 J/mol so that $f$ is in eV/atom, consistent with $k_B$ inside `SoretDiffusion` and with $\beta$ in $\mu\mathrm{m}^2/(\mathrm{eV\,s})$.

## Run

```bash
conda activate moose
cd /home/mzw/opt/moose_projects/palm/problem/uzr/s1_mohanty/2
../../../../palm-opt -i s1_mohanty.i
conda activate analysis
python plot_exodus_cx_profile.py s1_mohanty_out.e
```

## Result

Mass is conserved: $c_{\mathrm{avg}}=0.390000$. The run finished in 152 s. The profile is already near steady state by 10 days.

| $t$ (day) | $c(0)$, hot | $c(300\,\mu\mathrm{m})$, cold |
| --- | --- | --- |
| 0.2 | 0.506 | 0.333 |
| 1 | 0.548 | 0.145 |
| 5 | 0.583 | 0.027 |
| 10 | 0.589 | 0.002 |
| 30 | 0.589 | 0.002 |

Fig. 6 at 30 days is still evolving, with about $+11$ at.% at the hot end ($c\approx0.50$) and about $-7$ at.% at the cold end ($c\approx0.32$), and a shoulder near the cold end. This calculation does not reproduce that figure. The printed mobilities and Eq. (3) reach a much stronger steady segregation in a few days. No scale factor was introduced to slow it down.

![c(x)](s1_mohanty.png)

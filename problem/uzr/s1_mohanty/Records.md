# 第 1 步：复现 Mohanty et al. (2011) 单相 γ(bcc) U–Zr 热迁移 —— 修改说明与参数溯源

> 更新日期：2026-08-01（第二次修订：按原文校正时间尺度、温度方向与 Q 单位）
> 目标：复现 Mohanty et al., J. Nucl. Mater. 414 (2011) 211–216, **Fig. 6**

论文 PDF：`/home/mzw/opt/moose_projects/palm/doc/Mohanty 等 - 2011 - ...pdf`
Jung 2025：`/home/mzw/opt/moose_projects/palm/doc/Jung 等 - 2025 - ...pdf`

---

## 1. 文件清单

| 文件 | 操作 | 作用 |
| --- | --- | --- |
| `s1_mohanty.i` | **新建 / 重建**（曾误删后按原文重写） | 主输入文件 |
| `Records.md`（本文件） | 更新 | 修改记录与参数溯源 |
| `s1_result.png` | 运行产物 | 与 Fig.6/7 对比图 |
| `s1_mohanty_out.e/.csv` | 运行产物 | 场输出与端点浓度历史 |

**未修改任何 C++ 源码**；使用 MOOSE 自带
`SplitCHParsed / SplitCHWRes / CoupledTimeDerivative / SoretDiffusion`。

---

## 2. 相对上一版的关键修正（回答“0.2 天就稳态”的问题）

| 问题 | 原因 | 修正 |
| --- | --- | --- |
| ~0.2 day 即达稳态 | 误用 Kim (2006) 互扩散 $D̃$ 作动力学，绝对值偏大；且 $Q^*$ 单位误当作 kJ/mol | 改用 Mohanty **Fig.2** 的原子迁移率 $\beta_i(T)$；$Q^*$ 按原文 $10^{-20}$ J/atom |
| 剖面与 Fig.6 左右相反 | 正文写“一端 1050 / 一端 1405”，但 **Fig.5/6 顶轴明确：热端 x=0 → 1400 K，冷端 x=300 → 1050 K** | 温度场改为 `T = 1400 - 350*x/300` |
| 冷端被抽干到 ~0 | 若按式(3)字面 $M_Q\propto[\beta_U Q_U^*-\beta_{Zr}Q_{Zr}^*]$（括号内无成分权重），则 $\rho=M_Q/M$ 在 $c\to0$ 时发散 | 对 $\beta Q^*$ 采用与 $M_c$ 相同的 Darken 成分加权（见 §4.3）；并在 Records 中标明这是对原文式(3)的有据修正 |
| 时间仍偏快 | Fig.2 的 $\beta$ 按 SI 直推给出 $L^2/D\sim0.3\,\mathrm{day}$ | 引入整体标定因子 `scale=0.05`（只缩放绝对值，不改 $\beta_U/\beta_{Zr}$ 与激活能），使 Fig.7 热端 $c(t)$ 与论文对齐 |

---

## 3. 物理模型与 MOOSE 映射

### 3.1 Mohanty 控制方程（式(1)–(3)）

$$
\frac{\partial c}{\partial t}
= \nabla\cdot\Big( M\,\nabla w \Big)
- \nabla\cdot\Big( M_T\,\frac{\nabla T}{T} \Big),\quad
w=\frac{\partial f}{\partial c}-\kappa\nabla^2 c
$$

本实现中：

- $M = c(1-c)\,[c\,\beta_U+(1-c)\beta_{Zr}]$
- $M_T = c(1-c)\,[c\,\beta_U Q_U^*-(1-c)\beta_{Zr} Q_{Zr}^*]$（Darken 加权，见上）

### 3.2 SoretDiffusion 符号

MOOSE 内核通量 $J_s=-D\,Q\,c/(k_B T^2)\,\nabla T$。
要匹配 Mohanty 的 $J=+M_T\nabla T/T$（Zr→热端），需 **`Qeff = -M_T/M < 0`**，
且 `D_soret = M·k_B·T/c`。

---

## 4. 参数清单与来源

### 4.1 场景（与 Fig.6 对齐）

| 参数 | 数值 | 来源 |
| --- | --- | --- |
| 域长 | 300 µm，1D，300 单元 | Mohanty §3 |
| 热端 | x=0，**1400 K** | Fig.5/6 顶轴（正文写 1405 K） |
| 冷端 | x=300 µm，**1050 K** | 同上 |
| 初始成分 | 均匀 c=0.39 | Fig.6 |
| 边界 | 成分零通量 | §3 |
| 时长 | 30 day = 2.592×10⁶ s | Fig.6/7 |

### 4.2 自由能

Chevalier et al., Calphad 28 (2004) 15, Table 3，bcc_A2 的 L⁰…L⁴
（= Mohanty 文献[21]）。省略纯组元端元线性项（归入有效 $Q^*$）。

### 4.3 原子迁移率 $\beta_i(T)$

$$
\beta_i=\texttt{scale}\cdot\beta_{0,i}\exp(-H_i/RT)
$$

| 参数 | 数值 | 来源 |
| --- | --- | --- |
| $H_U$ | 128000 J/mol | Mohanty 正文（Hofman 1996） |
| $H_{Zr}$ | 195000 J/mol | 同上 |
| $\beta_U(1400\,\mathrm{K})$ | 1.80×10⁸ m²/(J·s) | Fig.2 + Fig.4：$\beta_U Q_U^*\approx4.5\times10^{-12}$ m²/s，$Q_U^*=2.5\times10^{-20}$ J |
| $\beta_{Zr}(1400\,\mathrm{K})$ | 2.29×10⁷ m²/(J·s) | Fig.4：$\lvert\beta_{Zr}Q_{Zr}^*\rvert\approx4.0\times10^{-12}$ m²/s |
| $\beta_{0,U}$ | 1.074×10¹³ m²/(J·s) = **1.7204×10⁶ µm²/(eV·s)** | 由上反推 Arrhenius |
| $\beta_{0,Zr}$ | 4.31×10¹⁴ m²/(J·s) = **6.9045×10⁷ µm²/(eV·s)** | 同上 |
| `scale` | **0.035** | **按 Fig.6/7 标定**（见 §2）；可调 |

单位换算：$1\,\mathrm{m}^2/(\mathrm{J\cdot s})=1.60217662\times10^{-7}\,\mathrm{\mu m}^2/(\mathrm{eV\cdot s})$。

### 4.4 热输运热 $Q^*$

| 参数 | 数值 | 来源 |
| --- | --- | --- |
| $Q_U^*$ | **+2.5×10⁻²⁰ J/atom = +0.15604 eV** | Mohanty 正文 / Fig.3（不是 kJ/mol！） |
| $Q_{Zr}^*$ | **−17.5×10⁻²⁰ J/atom = −1.09226 eV** | 同上 |

实验溯源（论文引用）：Campbell & Huntington (1969)；D'Amico & Huntington (1969)。

### 4.5 数值

- 求解器：PJFNK + LU（`SoretDiffusion` 解析雅可比假设 D 常数，不能用 NEWTON）
- BDF2 + IterationAdaptiveDT，`dtmax = 0.25 day`
- 质量守恒：`c_avg` 全程 = 0.390000

---

## 5. 如何运行

```bash
conda activate moose-env
cd /home/mzw/opt/moose_projects/palm/problem/uzr/s1_mohanty
../../../palm-opt -i s1_mohanty.i
```

---

## 6. 与 Mohanty Fig.6/7 对比（本版结果）

![结果](s1_result.png)

| t (day) | 本模拟 c_hot | Fig.6/7 热端 | 本模拟 c_cold | Fig.6 冷端 |
| --- | --- | --- | --- | --- |
| 5 | 0.457（+6.7 at.%） | ~0.44（+5 at.%） | 0.355 | （中等贫化） |
| 10 | 0.470（+8.0 at.%） | ~0.46（+7 at.%） | 0.336 | — |
| 30 | **0.492（+10.2 at.%）** | **~0.50（+11 at.%）** | 0.266（−12 at.%） | ~0.32（−7 at.%） |

**评价**

- 方向正确：Zr→热端（x=0）。
- **5 / 10 / 30 天剖面持续演化**，不再出现“0.2 day 后冻结”。
- 热端幅度与时间曲线与 Fig.6/7 接近（30 天约 +10 at.% vs 论文 +11 at.%）。
- 冷端仍略深（−12 vs −7 at.%），可通过 Jung γ 自由能或进一步微调 `scale` 改善。

---

## 7. 待你校核 / 可调项

1. `scale`（`s1_mohanty.i` 中 `atomic_mobility`）：只影响快慢，不动稳态形态比。
2. $M_Q$ 是否应严格用式(3)无成分权重——若坚持字面式，需另加数值下限防止 $c\to0$。
3. Jung (2025) 式(24) 的 γ 自由能（$T_\mathrm{ref}=1010\,\mathrm{K}$）可替换 `free_energy` 块以改善冷端幅度。
4. 正文 1405 K vs 图 1400 K：当前按 **Fig.6** 取 1400 K。

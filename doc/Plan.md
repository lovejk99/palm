# U–Zr 燃料温度梯度组分重分布：MOOSE 逐步复现与扩展规划

> 研究对象：U–Zr 金属燃料在温度梯度下的组分重分布  
> 数值框架：MOOSE，相场模块，split 形式 Cahn–Hilliard 方程  
> 规划范围：严格限制在 Mohanty et al. (2011)、Wen et al. (2022) 和 Jung et al. (2025) 明确使用的物理模型内  
> 推荐实施顺序：一维分项验证 $\rightarrow$ 一维多物理全耦合 $\rightarrow$ 二维几何扩展

## 1. 规划原则

### 1.1 复现目标

本项目的主要目标不是立即建立一个包含所有可能物理过程的“最大模型”，而是建立一个可追溯、可验证、可逐层扩展的 U–Zr 组分重分布计算框架。每一层只增加一个或一组紧密关联的物理过程，并在进入下一层之前完成：

1. 方程和符号检查；
2. 单位与量纲检查；
3. 守恒检查；
4. 网格和时间步收敛性检查；
5. 与论文曲线或实验数据的定量比较；
6. 参数来源、拟合参数和无法确定参数的分类记录。

最终框架应能回答三个层次的问题：

- Mohanty 2011：单相 $\gamma$-U–Zr 中，化学势梯度与温度梯度如何共同决定 Zr 迁移？
- Wen 2022：多相热力学、KKS 相平衡和辐照增强扩散如何改变宏观径向组分分布？
- Jung 2025：成分、相分数、孔隙率、钠渗入和热导率之间的反馈如何改变温度场与组分重分布？

### 1.2 不超出论文的耦合边界

本规划不在论文模型之外引入以下物理：

- 不显式求解空位和间隙原子速率理论方程；
- 不使用相场显式解析孔洞形核、长大、并合或迁移；
- 不令孔隙率直接改变扩散系数、自由能或力学性能；
- 不增加应力、裂纹、燃料–包壳相互作用或裂变气体输运；
- 不将二维扩展解释为新的物理模型；二维阶段仅检查几何和数值效应。

Wen 2022 中的辐照作用只通过辐照增强因子 $\xi$ 修改化学迁移率。Jung 2025 中的孔隙率只按相分数计算，并只通过有效热导率影响温度场。

### 1.3 单位策略

不同论文和复现目录允许采用不同单位制，不强制全项目统一：

- `s1_mohanty` 可继续使用 $\mathrm{\mu m}$–s–eV/atom–K；
- Wen 2022 和 Jung 2025 建议按论文使用 m–s–J/mol 或 J/m³–K；
- 不允许在同一个输入文件中未经显式换算混用 eV/atom、J/atom、J/mol 和 J/m³；
- 每个复现目录必须包含独立的 `UNITS.md` 或 README 单位章节；
- 每个 Material 属性名建议携带物理含义，必要时在注释中注明单位。

每个目录的单位审计至少回答：

1. 自由能是每原子、每摩尔还是每体积？
2. 化学势变量的单位是什么？
3. 迁移率与化学势梯度相乘后是否得到正确通量？
4. 热迁移率与 $\nabla T/T$ 相乘后是否与化学通量同量纲？
5. 坐标和梯度使用 m 还是 $\mathrm{\mu m}$？
6. 热导率、体积热源、密度和比热是否构成一致的热方程？

## 2. 三篇论文的模型关系

### 2.1 Mohanty et al. (2011)：单相热迁移基准

论文只考虑单相 bcc-$\gamma$ U–Zr。守恒变量为 Zr 原子分数 $c$，控制方程为

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

化学迁移率和热迁移率为

$$
M_c=\frac{1}{V_m}c(1-c)
\left[c\beta_U+(1-c)\beta_{\mathrm{Zr}}\right],
$$

$$
M_Q=\frac{1}{V_m}c(1-c)
\left[\beta_U\widetilde Q_U^*
-\beta_{\mathrm{Zr}}\widetilde Q_{\mathrm{Zr}}^*\right].
$$

原子迁移率服从

$$
\beta_i=\beta_{0,i}\exp\left(-\frac{H_i}{RT}\right).
$$

论文明确说明单相合金中的成分梯度能系数消失。因此严格论文分支应使用

$$
\kappa_c=0,
\qquad
\mu=\frac{\partial f}{\partial c}.
$$

当前 `s1_mohanty.i` 使用非零 `kappa_c = 1e-5` 作为数值正则化。该设置只能保留在校准/正则化分支中，并应通过 $\kappa_c\to0$ 收敛检查确认它没有控制 Fig. 6 的边界层或内部最低点。

关键基准：

- 一维长度 300 $\mathrm{\mu m}$；
- 初始成分 39 at.% Zr；
- 两端温度约 1050 K 和 1405 K；
- 5、10、30 day 的浓度剖面；
- 30 day 后热端约富集 11 at.% Zr，冷端约贫化 7 at.% Zr；
- 扩散偶初始成分为 39 at.% 与 22.5 at.% Zr；
- 热迁移通量在论文扩散偶算例中比化学势梯度通量大约 4–5 个数量级。

论文温度场由 $\nabla^2T=0$ 得到，并假定热导率与温度、成分无关。因此，这一阶段不应引入成分相关热导率。

### 2.2 Wen et al. (2022)：多相 KKS 与辐照增强扩散

Wen 2022 在 Mohanty 型守恒方程基础上加入：

- $\alpha$、$\beta$、$\gamma$ 三相的宏观表示；
- 一个非守恒 $\gamma$ 相序参量；
- KKS 相内浓度约束；
- 温度相关相界面迁移率；
- 辐照增强化学迁移率；
- 宏观相界面宽度与界面参数分析。

组分演化方程为

$$
\frac{\partial c}{\partial t}
=
\nabla\cdot\left(
M_c^I\nabla\frac{\delta F}{\delta c}
-M_T\frac{\nabla T}{T}
\right),
$$

其中

$$
M_c^I=\xi M_c,
\qquad
\xi=\frac{C_v^e+C_v^r}{C_v^e},
$$

$$
C_v^e=A_v\exp\left(-\frac{E_v^f}{k_BT}\right).
$$

论文采用 $E_v^f=1.20$ eV，并在约 1000 K、30 kW/m 条件下引用 $C_v^r=7.04\times10^{-7}$。论文没有在相场计算中显式求解缺陷浓度演化，因此复现时应把 $C_v^r$ 作为给定值或参数函数，而不是新增缺陷 PDE。

相序参量满足 Allen–Cahn 方程：

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

总自由能包含相体自由能、双阱势、成分梯度能和相场梯度能：

$$
F=\int_V
\left\{
\frac{(1-h)f_{\alpha\beta}+hf_\gamma}{V_m}
+\omega g(\gamma)
+\frac{\kappa_c}{2}\lvert\nabla c\rvert^2
+\frac{\kappa_\gamma}{2}\lvert\nabla\gamma\rvert^2
\right\}\,\mathrm dV.
$$

KKS 约束为

$$
c=(1-h)c_{\alpha\beta}+hc_\gamma,
$$

$$
\frac{\partial f_{\alpha\beta}}{\partial c_{\alpha\beta}}
=
\frac{\partial f_\gamma}{\partial c_\gamma}.
$$

Wen 式 (7)–(8) 的 $\alpha/\beta$ 权重与紧邻文字及 U–Zr 相图相反：按印刷公式，高温区会被标成 $\alpha$、低温区会被标成 $\beta$。因此必须同时保留：

- `literal`：逐字采用论文印刷公式；
- `phase_consistent`：交换 $\alpha$ 和 $\beta$ 权重，使 $T<935\,\mathrm K$ 为 $\alpha$、$T>935\,\mathrm K$ 为 $\beta$。

论文主图的相区顺序应作为两者取舍依据，不能在代码中无记录地交换。

Wen 2022 使用给定的径向温度剖面，没有求解成分依赖热导率的全耦合热传导。因此，辐照增强阶段应保持温度场为给定场，避免同时引入 Jung 2025 的热反馈。

### 2.3 Jung et al. (2025)：热传导、孔隙率与钠渗入反馈

Jung 2025 延续 KKS 多相相场模型，并求解热传导：

$$
\rho c_p\frac{\partial T}{\partial t}
=
\nabla\cdot\left[
k(P,P_{\mathrm{Na}},T,c)\nabla T
\right]
+\dot q(B).
$$

有效热导率为

$$
k(P,P_{\mathrm{Na}},T,c)
=P_c(P,P_{\mathrm{Na}},T,c)(1-P)^{1.5}k_0(T,c).
$$

基体热导率写成

$$
k_0(T,c)=a_1(c)+a_2(c)T+a_3T^2.
$$

对 U–Zr，$W_{\mathrm{Pu}}=0$，论文式 (28)–(30)应整理为

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

论文正文先称系数为 $a_0,a_1,a_2$，随后的公式却编号为 $a_1,a_2,a_3$。代码应采用一个一致命名，并在 README 中保留与论文编号的映射。

论文的孔隙率不是独立演化变量，而是相分数的代数函数：

$$
P=0.4\eta_\gamma+0.13(\eta_\alpha+\eta_\beta).
$$

钠热导率为

$$
k_{\mathrm{Na}}(T)
=93-0.0581(T-273.15)
+1.173\times10^{-5}(T-273.15)^2.
$$

孔隙修正因子按论文式 (31) 实现。孔隙率随燃耗的处理为：

- 0–1% burnup：孔隙率从 0 线性增大到目标值；
- 大于 1% burnup：孔隙率保持目标值；
- 孔隙稳定后钠渗入开始增加；
- 钠渗入在 1.5% burnup 达到目标值并保持不变。

体积热源按

$$
\dot q(B)=\frac{\mathrm{LHGR}(1-B)}{\pi r_0^2}
$$

近似。燃耗在总模拟时间内线性增大到 DP-81 的 3.6 at.% 或 DP-11 的 7.7 at.%。

这一模型中的反馈闭环为：

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

注意：论文明确指出孔隙率只影响热导率，不直接影响扩散动力学或相平衡。本复现应保持这一限制。

## 3. 当前第一步实现的定位与需要修正的事项

当前目录 `palm/problem/uzr/s1_mohanty` 已实现：

- 一维 300 $\mathrm{\mu m}$ 网格；
- Zr 原子分数 $c$ 和化学势 $w$ 的 split CH；
- 线性给定温度场；
- CALPHAD 型 $\gamma$ 相自由能；
- 温度相关原子迁移率；
- MOOSE 内置 `SoretDiffusion`；
- 5、10、30 day 输出；
- 质量平均值、端点浓度和空间剖面后处理。

这一工作可作为“现象复现分支”，但在进入后续阶段前，应建立“论文原式分支”，原因如下。

### 3.1 当前模型包含两个人为修正

当前输入文件使用：

- 原子迁移率整体缩放因子 `scale = 0.035`；
- 对 Mohanty 式 (3) 增加了成分权重。

当前 README 的修订说明中还残留 `scale = 0.05`，而参数表和实际输入文件使用 `0.035`。后续整理时应以输入文件为准，并修正文档中的旧值。当前 `kappa_c` 注释写成“热导率”，实际应为成分梯度能系数。

这两个处理改善了 Fig. 6/7 的拟合，但它们不是 Mohanty 2011 原式。建议保留两个并行算例：

- `paper_equation.i`：严格使用论文式 (2)–(4)，不加成分权重，不先拟合时间尺度；
- `calibrated.i`：保留当前稳定化和拟合版本。

所有图中必须明确标记使用的是“原式”还是“校准式”，避免将拟合结果误称为无参数复现。

### 3.2 当前 `SoretDiffusion` 不适合后续强耦合

MOOSE 内置 `SoretDiffusion` 的残差为

$$
R_T=
\int_\Omega
\frac{DQc}{k_BT^2}\nabla T\cdot\nabla\psi\,\mathrm dV.
$$

当前通过

$$
D_{\mathrm{Soret}}=\frac{M_ck_BT}{c},
\qquad
Q_{\mathrm{eff}}=-\frac{M_T}{M_c}
$$

将其映射为目标热迁移通量。代数上映射可行，但内置 Kernel 的 Jacobian 只显式处理 $c$ 和 $T$，没有包含 $D(c,T)$ 与 $Q(c,T)$ 的材料导数。当前使用 PJFNK 可以规避解析 Jacobian 不完整带来的部分问题，但在 Jung 2025 的 $c$–$\phi$–$T$ 强耦合问题中会导致：

- Newton 收敛速度下降；
- 预条件矩阵缺少关键耦合块；
- 难以区分物理刚性与 Jacobian 不完整；
- 后续材料属性依赖相内浓度时更难维护。

建议在 PALM 中新增一个直接使用 $M_T$ 的 AD Kernel，例如 `ADThermalMigration`，其弱式残差为

$$
R_{\mathrm{TM}}
=-\int_\Omega
\frac{M_T}{T}\nabla T\cdot\nabla\psi\,\mathrm dV.
$$

该对象应：

- 继承 `ADKernel`；
- 耦合温度变量 `T`；
- 读取 AD 材料属性 `thermal_mobility`；
- 由自动微分生成对 $c$、$\phi$、$T$、$c_{\alpha\beta}$ 和 $c_\gamma$ 的完整 Jacobian；
- 不在 Kernel 内硬编码 Boltzmann 常数或单位制；
- 不再把 $M_T$ 人为拆成 `D_soret` 和 `Qeff`。

`SoretDiffusion` 可以继续用于第一步对照，但不建议修改 MOOSE 上游源码。新增对象应放在 PALM：

- `include/kernels/ADThermalMigration.h`
- `src/kernels/ADThermalMigration.C`
- `test/tests/kernels/ad_thermal_migration/`

PALM 已通过 `Registry::registerObjectsTo(f, {"palmApp"})` 注册本应用对象，因此在新类源文件中使用

```cpp
registerMooseObject("palmApp", ADThermalMigration);
```

即可纳入应用。

### 3.3 第一步完成标准

在进入 Wen 2022 前，应完成：

- 原式与校准式的并列比较；
- 在恒温条件下关闭热迁移，确认均匀成分保持不变；
- 令 $Q_U^*=Q_{\mathrm{Zr}}^*=0$，确认热迁移项严格为零；
- 反转温度梯度，确认迁移方向反转；
- 检查 $\int_\Omega c\,\mathrm dV$ 的相对漂移；
- 至少使用 150、300、600 单元进行网格收敛；
- 至少使用三组 `dtmax` 进行时间步收敛；
- 对热端、冷端及剖面误差分别量化；
- 增加化学通量、热迁移通量和总通量的输出。

## 4. 推荐目录结构

建议在 `palm/problem/uzr` 下建立独立、可运行的阶段目录：

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

每个目录至少包含：

- 主输入文件；
- 最小验证输入文件；
- `README.md`；
- `UNITS.md`；
- 参数来源文件；
- 自动提取和绘图脚本；
- 网格/时间步收敛脚本；
- MOOSE tests 入口或指向 PALM `test/tests` 的测试说明。

## 5. 分阶段复现路线

## 阶段 0：建立可重复验证基础设施

### 目标

在继续增加物理之前，将运行、采样、论文曲线数字化和误差计算自动化。

### 工作内容

1. 将每篇论文的目标曲线数字化为 CSV，并记录图号、坐标轴单位和人工读取误差。
2. 统一输出以下指标：
   - 总 Zr 量；
   - 最小和最大 $c$；
   - 中心与表面 $c$；
   - 中心与表面温度；
   - 各相体积分数；
   - 相区边界位置；
   - $L_1$、$L_2$ 和最大相对误差。
3. 为每个算例保存 MOOSE 版本、PALM commit、输入文件哈希和命令行。
4. 将绘图脚本改为读取指定文件，而不是默认选择目录中最新的 Exodus 文件。
5. 增加非物理值检查：
   - $0\le c\le1$；
   - $0\le\phi\le1$；
   - $0\le P<1$；
   - $k>0$；
   - $T>0$；
   - 相分数之和为 1。

### 通过标准

- 同一输入重复运行得到一致结果；
- 后处理脚本不依赖人工修改；
- 误差指标可由单条命令生成；
- 质量守恒误差建议先控制在 $10^{-8}$ 相对量级以内；若由于求解容差不能达到，应记录误差随容差的变化。

## 阶段 1：完成 Mohanty 2011 单相均匀合金基准

### 子任务 1A：方程级验证

建立小尺寸制造解或简化算例，分别验证：

- `SplitCHParsed`/`SplitCHWRes` 的化学扩散部分；
- `ADThermalMigration` 的热迁移部分；
- 两项叠加后的总通量；
- 温度反转、热迁移热符号反转和零温度梯度极限。

### 子任务 1B：论文 Fig. 5–7

依次运行：

1. 常数 $\beta_i$、常数 $Q_i^*$；
2. 常数 $\beta_i$、温度相关 $Q_i^*$；
3. 温度相关 $\beta_i$、常数 $Q_i^*$；
4. 5、10、30 day 瞬态；
5. 热端浓度随时间变化。

如果温度相关 $Q_i^*$ 的完整数据无法从论文获得，应把该子项标为“不可唯一复现”，不要自行构造函数。

### 验证量

- 迁移方向：Zr 向热端迁移；
- 30 day 热端与冷端成分变化；
- 5 day 初期快速迁移与后期减缓；
- 端点值误差和全剖面误差；
- 质量守恒；
- 150/300/600 单元网格收敛；
- 三组最大时间步收敛。

## 阶段 2：Mohanty 2011 扩散偶与通量分解

### 目标

复现论文 Fig. 8，并确认热迁移与化学扩散的竞争关系。

### 算例

1. 1225 K 等温扩散偶；
2. 温度梯度与初始 Zr 浓度梯度同向；
3. 温度梯度与初始 Zr 浓度梯度反向。

初始成分为 39 at.% 和 22.5 at.% Zr。阶跃初值建议先用窄的平滑双曲正切过渡，避免有限元初始不连续导致的网格依赖；过渡宽度必须小于论文尺度且进行敏感性分析。

### 必需输出

定义并输出

$$
\mathbf J_c=-M_c\nabla w,
\qquad
\mathbf J_T=M_T\frac{\nabla T}{T},
\qquad
\mathbf J_{\mathrm{total}}=\mathbf J_c+\mathbf J_T.
$$

验证论文所述 $\lvert J_T/J_c\rvert$ 的数量级，并检查其空间与时间变化。不要只比较最终浓度曲线。

## 阶段 3：无辐照 KKS 多相骨架

### 目标

在加入辐照增强之前，独立验证 $\alpha/\beta/\gamma$ 热力学、KKS 约束和 Allen–Cahn 相演化。

### 推荐变量

- `c`：总 Zr 原子分数，守恒变量；
- `w`：split CH 化学势；
- `phi`：$\gamma$ 相序参量；
- `c_ab`：$\alpha\beta$ 混合相内 Zr 浓度；
- `c_g`：$\gamma$ 相内 Zr 浓度；
- `T`：本阶段为给定 AuxVariable；
- `eta_a`、`eta_b`、`eta_g`：输出用 AuxVariable 或材料属性。

### MOOSE 对象映射

优先复用 phase-field 模块的 KKS 对象：

- `KKSPhaseConcentration`
- `KKSPhaseChemicalPotential`
- `KKSSplitCHCRes`
- `KKSACBulkF`
- `KKSACBulkC`
- `KKSGlobalFreeEnergy`
- `TimeDerivative`
- `ACInterface`

可参考：

- `moose/modules/phase_field/test/tests/KKS_system/two_phase.i`
- `moose/modules/phase_field/test/tests/KKS_system/kks_example_split.i`
- `moose/modules/phase_field/examples/kim-kim-suzuki/`

自由能可先用 `DerivativeParsedMaterial` 实现，但应把 $\ln c$ 和 $\ln(1-c)$ 的定义域问题作为显式测试。不要通过过大的数值截断掩盖相内浓度求解失败。

`KKSSplitCHCRes` 负责 KKS 化学势关系，本身不自动加入 Wen 的 $\kappa_c\nabla^2c$。Wen 分支若使用非零 $\kappa_c$，需用 `MatDiffusion` 或等价 Kernel 明确加入梯度项；Jung 严格分支不加入该项。若按 Mohanty 式 (1) 保留 $2\kappa_c$ 的定义，则传给 split Kernel 的系数应是论文 $\kappa_c$ 的两倍；在严格 $\kappa_c=0$ 分支中该因子不影响结果。

### 分步验证

1. 单点自由能和一阶、二阶导数与独立 Python/SymPy 结果比较；
2. 给定 $c,T,\phi$ 时，检查 KKS 两个约束残差；
3. 等温两相平衡，检查化学势相等；
4. 无热迁移下的界面平衡；
5. 给定温度场下的三相区位置；
6. 检查 $\eta_\alpha+\eta_\beta+\eta_\gamma=1$；
7. 检查界面宽度随 $\kappa_\phi/\omega$ 的变化。

### 重要歧义

Wen 2022 与 Jung 2025 的插值函数、势垒函数排印可能存在差异。特别是 Jung 式 (8) 印为 $g(\phi)=\phi^2(1-\phi^2)$，而常见双阱形式为 $\phi^2(1-\phi)^2$。实施前应：

- 按原文逐字建立一个版本；
- 按其引用模型建立一个参考版本；
- 检查两个版本是否都在 $\phi=0,1$ 具有预期极小值；
- 若原文形式不能产生合理双阱，应记录为论文排印或定义歧义。

## 阶段 4：Wen 2022 辐照增强扩散

### 目标

只按 Wen 论文的有效迁移率方法加入辐照，不引入额外缺陷方程。

### 实现顺序

1. 计算热平衡空位浓度 $C_v^e(T)$；
2. 使用论文给定的 $C_v^r$；
3. 计算 $\xi(T)$；
4. 令 `chemical_mobility_irradiated = xi * chemical_mobility`；
5. 热迁移率仍按论文式 (29)–(31)计算，不擅自乘以 $\xi$，除非经论文方程和文字再次确认；
6. 对比有辐照与无辐照结果。

### 建议材料属性

- `thermal_vacancy_concentration`
- `irradiation_vacancy_concentration`
- `irradiation_enhancement_factor`
- `chemical_mobility_unirradiated`
- `chemical_mobility_irradiated`

这些关系可先用 `ADParsedMaterial` 或 AD 解析材料实现，不必立即编写 C++ Material。只有在分段函数、数值保护或性能成为问题时再新增 `ADU_ZrIrradiationMobilityMaterial`。

### 论文复现设置

- 一维长度 2170 $\mathrm{\mu m}$；
- 2170 个有限差分单元为原论文设置，MOOSE 中应做独立网格收敛，不机械要求同样单元数；
- 中心约 988 K，表面 900 K；
- 初始 Zr 原子分数约 0.225；
- 零组分通量和零相场通量；
- 论文显式时间步 1 s，总时间 $1.3\times10^5$ s；
- MOOSE 可使用隐式自适应时间步，但必须证明结果不依赖更大的时间步。

### 目标图

- Wen Fig. 2：无辐照增强；
- Wen Fig. 3：有辐照增强；
- Wen Fig. 4：与其他模型和实验数据比较；
- Wen Fig. 5：$\kappa_\gamma$ 和 $\omega$ 敏感性；
- Wen Table 2：界面能和宏观界面宽度。

### 通过标准

- 有辐照模型明显增强低温区动力学；
- 能形成论文描述的中心富 Zr 区、中间 “Zr well” 和外侧贫化区；
- 三相区顺序与论文一致；
- KKS 约束、质量守恒和相分数归一性通过；
- 对 $\xi$、$\kappa_\gamma$、$\omega$ 分别做单参数敏感性分析；
- 明确区分文献给定参数、引用参数和试配参数。

## 阶段 5：Jung 2025 热传导基线

### 目标

先建立不含孔隙和钠渗入的 $c$–$\phi$–$T$ 全耦合基线，再添加热导率修正。

### 热方程对象

建议使用：

- `ADHeatConductionTimeDerivative`
- `ADHeatConduction`
- `ADBodyForce` 或等价 AD 热源 Kernel
- 表面 `ADDirichletBC` 或论文一致的 Dirichlet BC

材料属性至少包括：

- `density`
- `specific_heat`
- `thermal_conductivity`
- `volumetric_heat_source`

若论文没有给出 $\rho$ 和 $c_p$，而最终结果实质为准稳态热场，应先做稳态热方程验证，并把瞬态热参数标为缺失。不能仅为运行瞬态热方程而无来源地指定 $\rho c_p$。

### 成分到重量分数换算

Jung 热导率关系使用 Zr 重量分数，而相场变量是 Zr 原子分数。应使用

$$
W_{\mathrm{Zr}}
=
\frac{cM_{\mathrm{Zr}}}
{cM_{\mathrm{Zr}}+(1-c)M_U},
$$

$$
W_U=1-W_{\mathrm{Zr}}.
$$

论文式 (34)–(35) 的文本提取存在循环定义或排印不清，代码中必须使用可独立验证的原子分数–重量分数换算，并在 README 中说明。

### Jung 严格模型中的成分梯度项

Jung 式 (5)–(11)只给出相序参量梯度能

$$
\frac{\kappa_\phi}{2}\lvert\nabla\phi\rvert^2,
$$

没有报告 $\kappa_c\lvert\nabla c\rvert^2/2$。因此：

- `jung_literal` 中应设置 $\kappa_c=0$；
- 若沿用 Wen 的非零 $\kappa_c$，必须命名为 `wen_inherited` 或敏感性分支；
- 两个分支都可使用 split CH 软件结构，但 $\kappa_c=0$ 时方程物理上退化为二阶非线性扩散。

### 子步骤

1. 常数 $k$ 的稳态圆柱/平板解析解验证；
2. 只令 $k=k(T)$；
3. 令 $k=k(T,c)$；
4. 打开 $c$–$\phi$–$T$ 单向耦合：温度影响相场，但相场不反馈热导率；
5. 打开双向耦合：$c$ 影响 $k$，$k$ 影响 $T$；
6. 检查单体求解与分裂求解的差异。

论文没有报告总仿真时间，却令燃耗在总仿真时间内线性增长，并使用 Hirschhorn 优化迁移率压缩约 284 EFPD 的真实辐照过程。因此 Jung 结果可用于终态剖面复现，但在获得作者输入文件或补充数据前，不能用于验证真实时间动力学。

### 几何注意事项

论文把一维区间左端视为燃料中心、右端视为表面，但一维有限元实现可能使用平板算子。应建立两个明确区分的输入：

- `paper_1d.i`：尽量复现论文实际一维实现；
- `cylindrical_check.i`：使用轴对称径向热传导或等价径向权重，检查几何项影响。

几何修正结果不能替代论文复现结果，应作为单独敏感性研究。

## 阶段 6：Jung 2025 孔隙率耦合

### 目标

按论文代数模型把相分数映射为孔隙率，并只修正热导率。

### 实现

孔隙率建议作为 AD 材料属性，而不是非线性变量：

$$
P_{\mathrm{phase}}
=0.4\eta_\gamma+0.13(\eta_\alpha+\eta_\beta).
$$

再乘以燃耗增长函数：

$$
P(B)=s_P(B)P_{\mathrm{phase}},
$$

其中

$$
s_P(B)=
\begin{cases}
B/0.01, & 0\le B<0.01,\\
1, & B\ge0.01.
\end{cases}
$$

此处 $B$ 应统一使用无量纲分数，避免把 1% 写成 1 而不是 0.01。

### 对比算例

- 常数 $P=0$；
- 常数 $P=0.1$；
- 常数 $P=0.2$；
- 常数 $P=0.3$；
- 相分数相关 Baseline 孔隙率。

### 目标图与指标

复现 Jung Fig. 5 的：

- Zr 浓度；
- 温度；
- 热导率；
- $\gamma$ 相分数。

额外输出：

- 孔隙率剖面；
- Zone I/II/III 边界位置；
- 中心温度随孔隙率的变化；
- 中心 Zr 浓度随孔隙率的变化。

论文摘要给出的量级检查包括：忽略孔隙率时，中心 Zr 浓度约降低 0.17 at.%，中心温度约降低 36 K。具体比较必须以复现基线定义为准。

## 阶段 7：Jung 2025 钠渗入

### 目标

按论文式 (31)–(32)加入钠对孔隙有效热导率的修正。

### 实现顺序

1. 独立验证 $k_{\mathrm{Na}}(T)$；
2. 独立验证 $P_c$ 在 $P_{\mathrm{Na}}=0$ 时的极限；
3. 检查 $k_{\mathrm{Na}}/k_0$ 在全部温度和成分范围内的数值；
4. 加入 5%、25%、50% 钠渗入；
5. 加入随燃耗增长的 $P_{\mathrm{Na}}(B)$；
6. 打开全部反馈。

钠渗入百分数必须明确是 0–1 的分数还是 0–100 的百分数。建议代码内部统一用 0–1，并在输入层转换。

论文还没有明确 $P_{\mathrm{Na}}$ 表示“被钠填充的孔隙比例”还是“钠占燃料总体积的绝对体积分数”。应建立两个显式分支：

$$
P_{\mathrm{Na}}^{\mathrm{literal}}=f_{\mathrm{Na}},
$$

$$
P_{\mathrm{Na}}^{\mathrm{pore}}=f_{\mathrm{Na}}P.
$$

主复现先采用论文式 (31) 的字面定义，另一分支用于不确定性分析。不得在同一输入文件中根据局部孔隙率自动切换定义。

### 目标图

复现 Jung Fig. 6 的：

- Zr 浓度；
- 温度；
- 热导率；
- $\gamma$ 相分数。

论文给出的量级检查为：50% 钠渗入相对 Baseline 使中心 Zr 浓度约低 0.7 at.%，中心温度约降低 20 K。应同时比较整条剖面，不能只拟合中心点。

## 阶段 8：DP-81 与 DP-11 验证

### DP-81

- 半径/一维长度：2.17 mm；
- 表面温度：900 K；
- LHGR：24 kW/m；
- 目标燃耗：3.6 at.%；
- 初始 Zr：约 23 at.%；
- 初始 $\gamma$ 相分数：中心 0.5，表面 0.4；
- 初始温度：中心 990 K，表面 900 K；
- Baseline 钠渗入：5%。

### DP-11

- 表面温度：910 K；
- LHGR：25 kW/m；
- 目标燃耗：7.7 at.%；
- 其余缺失初始条件应从论文图、引用报告或补充资料中追溯，不直接沿用 DP-81。

### 验证顺序

1. 温度场；
2. 热导率场；
3. 三相分数与相区边界；
4. Zr 浓度；
5. 孔隙率；
6. 中心量随时间/燃耗的演化。

应对 PIE 数据进行数字化，并将实验散点与模拟插值到相同归一化半径 $r/r_0$ 后计算误差。

## 阶段 9：二维扩展

### 目标

在不增加论文外新物理的前提下，检查一维假设、几何离散和非均匀边界对结果的影响。

### 顺序

1. 二维矩形域重现一维解，用于代码回归；
2. 二维轴对称 $r$–$z$ 域，轴向条件均匀；
3. 允许轴向功率或表面温度变化；
4. 比较中心线、不同轴向位置和体积平均结果；
5. 检查二维网格各向异性和并行可扩展性。

二维结果的第一项验收标准是：当边界和初值沿第二维均匀时，二维解应退化为一维解。

## 6. 建议新增的 PALM 代码

## 6.1 必需自定义对象

### `ADThermalMigration`

用途：直接离散

$$
-\nabla\cdot\left(M_T\frac{\nabla T}{T}\right)
$$

并保留所有材料依赖的自动微分 Jacobian。

最小参数：

- `variable`：split 质量平衡方程所在变量；
- `T`：温度变量；
- `thermal_mobility`：AD 材料属性名。

测试：

- 常数 $M_T$ 制造解；
- $M_T(c,T)$ Jacobian 测试；
- 温度反转；
- 与 `SoretDiffusion` 在常系数映射下对比。

## 6.2 优先用解析材料、必要时再自定义

以下关系优先由 `ADParsedMaterial`、`DerivativeParsedMaterial`、`ParsedFunction` 或 `PiecewiseLinear` 实现：

- $\alpha\beta/\gamma$ 自由能；
- 各相化学迁移率；
- 各相热迁移率；
- $\xi(T)$；
- 相分数；
- 孔隙率；
- $k_0(T,c)$；
- $k_{\mathrm{Na}}(T)$；
- $P_c$；
- burnup 和 LHGR 历史。

当出现以下情况时再编写 C++ Material：

- 复杂分段表达式难以审查；
- 同一公式在多个算例重复；
- 需要集中执行单位换算和范围检查；
- 解析表达式导致明显性能问题；
- 需要清晰的错误消息和参数合法性检查。

届时可考虑：

- `ADUZrKKSFreeEnergyMaterial`
- `ADUZrMobilityMaterial`
- `ADUZrThermalConductivityMaterial`

但不建议在方程尚未通过单点验证前一次性实现这些类。

## 6.3 不建议修改的 MOOSE 上游代码

不建议直接修改：

- `moose/modules/phase_field/src/kernels/SoretDiffusion.C`
- MOOSE KKS Kernels；
- MOOSE HeatConduction Kernels。

原因是上游修改会增加升级、回归和结果追溯成本。PALM 自定义对象可以更准确地表达论文方程，同时保持 MOOSE checkout 可更新。

PALM 当前没有自定义 C++ 物理对象，`HEAT_TRANSFER` 和 `PHASE_FIELD` 已启用，`POROUS_FLOW` 未启用。Jung 的孔隙率只是代数材料属性，不需要为此启用 PorousFlow；只有未来真正求解孔隙流动时才需要重新评估模块依赖。

## 7. 求解策略

### 7.1 第一阶段

当前 `PJFNK + LU` 可继续用于 Mohanty 对照。新增 AD Kernel 后，应尝试完整 Newton：

- 小问题：`NEWTON + LU`；
- 较大问题：PJFNK 或 JFNK 配合字段分裂预条件；
- 对 $c,w,\phi,c_{\alpha\beta},c_\gamma,T$ 建立变量缩放。

### 7.2 强耦合阶段

Jung 模型建议先采用单体求解，因为温度、相分数和组分构成明显反馈。若单体系统难以收敛，再尝试：

1. 每个时间步先求热稳态，再求相场；
2. 热–相场 Picard 外迭代；
3. 最后再比较分裂方案与单体方案的差异。

不能仅因分裂方案容易收敛就默认其等价。需要检查时间步缩小时两种方案是否趋于同一解。

### 7.3 数值保护

允许的数值保护包括：

- 对初始阶跃做有记录的平滑；
- 对时间步进行自适应截断；
- 对非线性变量设置合理缩放；
- 使用变量边界或线搜索避免 $c$ 离开 $(0,1)$。

不建议：

- 无记录地截断自由能；
- 使用很大的梯度能系数掩盖网格不足；
- 用迁移率缩放同时拟合曲线形状和时间尺度；
- 在同一次拟合中自由调整多个高度相关参数。

## 8. 验证与不确定性量化

### 8.1 守恒

对无通量边界，应计算

$$
\epsilon_m(t)
=
\frac{\left|\int_\Omega c(t)\,\mathrm dV
-\int_\Omega c(0)\,\mathrm dV\right|}
{\int_\Omega c(0)\,\mathrm dV}.
$$

一维平板、圆柱和二维模型的积分测度必须与几何一致。

### 8.2 网格收敛

每个关键阶段至少使用三组网格。比较：

- 中心和表面 $c$；
- 最大/最小 $c$；
- Zone 边界；
- 中心温度；
- 全剖面 $L_2$ 差异。

界面模型还必须保证每个弥散界面内有足够单元。网格收敛不能只看全局平均量。

### 8.3 时间步收敛

使用相同空间网格和求解容差，至少将最大时间步连续减半两次。对燃耗驱动的分段函数，时间步必须能解析 1% 和 1.5% burnup 附近的转折点。

### 8.4 论文曲线误差

对数字化论文曲线 $y_i^{\mathrm{ref}}$ 和模拟曲线 $y_i$，建议报告：

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

同时报告：

- 数字化点数；
- 插值方法；
- 坐标归一化方法；
- 图像读取不确定度；
- 是否对参数进行了拟合。

### 8.5 参数分类

每个参数标为以下之一：

- `direct`：论文直接给出；
- `cited`：来自论文引用文献；
- `digitized`：从图中读取；
- `derived`：由已知量换算；
- `calibrated`：为拟合结果而调整；
- `assumed`：缺失时采用的假设。

只有 `direct`、`cited` 和可验证的 `derived` 参数能够被称为无拟合输入。

## 9. 论文中的关键缺口与风险

### Mohanty 2011

- 式 (3) 按排印不含成分权重，可能导致与直觉或数值表现不同的极限行为；
- 论文说明单相中 $\kappa_c=0$，当前算例的非零梯度能是额外数值正则化；
- 原子迁移率前因子未在正文中完整列表，需要从图中读取或追溯引用；
- 自由能来自商业数据库，不能完全由论文正文唯一恢复；
- 文本写 1405 K，而图轴可能显示 1400 K；
- 论文式 (6) 的 $c_x=0,c_{xxx}=0$ 不严格等价于含热迁移项的总无通量边界；
- 无法仅根据最终曲线唯一确定时间尺度修正。

### Wen 2022

- $C_v^r$ 来自外部速率理论程序，而非论文内可重复推导；
- 辐照增强因子缺少完整辐照条件依赖；
- $\alpha/\beta$ 相分数和迁移率插值权重与论文相区文字相反；
- $L_0$ 未给出，$\kappa_c$ 只给数量级；
- 部分自由能表达使用固定温度或排印不清；
- 摘要称讨论辐照对热扩散的影响，但控制方程只令 $M_c^I=\xi M_c$；
- 显式有限差分和 MOOSE 隐式有限元的时间离散不同；
- 相界面参数含明显试配过程；
- 燃料表面附近 Zr 富集无法由该模型解释。

### Jung 2025

- 采用 Hirschhorn 优化迁移率以压缩真实 284 EFPD 时间，严格时间尺度不可直接解释；
- 总仿真时间未报告，因此燃耗随时间的斜率不可唯一恢复；
- $\rho$ 和 $c_p$ 等瞬态热参数在正文参数表中不完整；
- $\beta_0$、$\widetilde Q_U^*$、$\widetilde Q_{\mathrm{Zr}}^*$ 和优化迁移率缩放不完整；
- 热导率系数编号和排印可能存在错位；
- 原子分数到重量分数的公式文本存在歧义；
- $\alpha/\beta$ 切换方向、双阱函数和钠渗入量定义均存在歧义；
- 论文未报告成分梯度能，严格分支应取 $\kappa_c=0$；
- 孔隙率模型为相区平均代数模型，不是孔结构演化；
- LHGR 与燃耗均使用简化线性历史；
- 论文承认表面 $\alpha$ 相和表面浓度存在偏差；
- 论文中的“一维”是否包含严格圆柱几何项需要通过输入文件或作者代码进一步确认。

这些缺口应成为复现报告的一部分，而不是通过未说明的参数调整隐藏。

## 10. 里程碑与完成判据

### M1：单相热迁移可信

- Mohanty Fig. 5–7 与扩散偶趋势复现；
- 原式与校准式明确分离；
- 质量、网格、时间步检查通过；
- `ADThermalMigration` 测试通过。

### M2：多相 KKS 可信

- 单点自由能和导数通过；
- KKS 约束残差收敛；
- 相分数归一；
- 界面宽度网格无关；
- 无辐照基线稳定。

### M3：Wen 辐照增强复现

- 无辐照/有辐照差异复现；
- 三个相区和 Zr well 出现；
- $\xi$ 及界面参数敏感性完成；
- 不引入论文外缺陷方程。

### M4：Jung 热反馈基线可信

- 热传导解析基准通过；
- $k(T,c)$ 单向和双向耦合逐层通过；
- DP-81 Baseline 的温度、相分数、成分达到可量化一致。

### M5：孔隙率和钠渗入复现

- Jung Fig. 5 和 Fig. 6 的四类剖面均完成；
- 中心温度与中心 Zr 变化量级一致；
- 燃耗分段函数经过时间步收敛检查；
- 孔隙率只通过热导率耦合。

### M6：二维框架完成

- 均匀第二维条件下退化到一维解；
- 轴对称几何影响被独立量化；
- 二维计算保持质量守恒和网格收敛；
- 一维论文复现结果仍作为主验证基准。

## 11. 推荐立即执行的下一步

1. 冻结当前 `s1_mohanty.i` 为校准版本。
2. 建立 Mohanty 论文原式版本，去掉迁移率缩放和额外成分权重。
3. 在 PALM 中实现并测试 `ADThermalMigration`。
4. 为 Fig. 6/7 建立数字化参考 CSV 和自动误差脚本。
5. 完成 Mohanty 网格、时间步和通量分解验证。
6. 复现 Mohanty 扩散偶后，再建立 KKS 多相骨架。
7. KKS 无辐照模型通过后才加入 Wen 的 $\xi$。
8. Wen 阶段通过后才把温度从给定场升级为 Jung 的热传导变量。
9. 按“常数热导率 $\rightarrow k(T,c)$ $\rightarrow$ 孔隙率 $\rightarrow$ 钠渗入”的顺序打开 Jung 耦合。
10. 完成全部一维验证后再进入二维。

## 12. 主要参考文件

- Mohanty et al. (2011)：`palm/doc/Mohanty 等 - 2011 - Thermotransport in γ(bcc) U–Zr alloys A phase-field model study.pdf`
- Wen et al. (2022)：`palm/doc/Wen 等 - 2022 - A phase-field model with irradiation-enhanced diffusion for constituent redistribution in U-10wt%Zr.pdf`
- Jung et al. (2025)：`palm/doc/Jung 等 - 2025 - Investigating Constituent Redistribution in U-Zr Metallic Fuels A Phase-Field Approach Incorporatin.pdf`
- 当前第一步说明：`palm/problem/uzr/s1_mohanty/ReadMe.md`
- 当前第一步输入：`palm/problem/uzr/s1_mohanty/s1_mohanty.i`
- MOOSE Soret Kernel：`moose/modules/phase_field/src/kernels/SoretDiffusion.C`
- MOOSE KKS 示例：`moose/modules/phase_field/test/tests/KKS_system/`
- MOOSE 热传导示例：`moose/modules/heat_transfer/test/tests/`

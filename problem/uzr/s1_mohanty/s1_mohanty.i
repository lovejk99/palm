# ==============================================================================
# Step 1: 复现 Mohanty et al., J. Nucl. Mater. 414 (2011) 211-216, Fig.6
# 单相 gamma(bcc) U-Zr 热迁移
#
# 控制方程（Mohanty 式(1)–(3)）:
#   dc/dt = div[ M * grad(w) ] - div[ M_T * grad(T)/T ]
#   w     = df/dc - kappa * lap(c)
#   M   = c(1-c)[ c*beta_U + (1-c)*beta_Zr ]
#   M_T = c(1-c)[ beta_U*Qstar_U - beta_Zr*Qstar_Zr ]
#
# 单位: 长度 um, 时间 s, 能量 eV/atom, 温度 K
# ==============================================================================

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 300
  xmin = 0
  xmax = 300 # um（Mohanty 2011 §3 / Fig.6）
  elem_type = EDGE
[]

[Variables]
  [c] # Zr 原子分数
    order = FIRST
    family = LAGRANGE
  []
  [w] # 化学势 df/dc (eV/atom)
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [c_IC]
    # U-39 at.% Zr（Mohanty Fig.6）
    type = ConstantIC
    variable = c
    value = 0.39
  []
[]

[AuxVariables]
  [T]
  []
[]

[AuxKernels]
  [T_profile]
    # Mohanty Fig.5/6 顶轴: 热端 x=0 → 1400 K, 冷端 x=300 → 1050 K
    # （正文写 1405 K，图中刻度为 1400 K；此处与 Fig.6 对齐取 1400 K）
    type = FunctionAux
    variable = T
    function = '1400 - (1400 - 1050)/300 * x'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Kernels]
  [c_dot]
    type = CoupledTimeDerivative
    variable = w
    v = c
  []
  [w_res]
    type = SplitCHWRes
    variable = w
    mob_name = M
  []
  [w_res_soret]
    # D_soret*Qeff*c/(kB*T^2)*grad(T)  映射为  M_T*grad(T)/T
    # MOOSE 中 Soret 通量 J_s = -D*Q*c/(kB*T^2)*grad(T)
    # Mohanty 中 J     = +M_T*grad(T)/T  （Zr→热端）
    # 故需 Qeff < 0（见 Materials/Qeff_mat）
    type = SoretDiffusion
    variable = w
    c = c
    T = T
    diff_name = D_soret
    Q_name = Qeff
  []
  [c_res]
    type = SplitCHParsed
    variable = c
    f_name = f_bulk
    kappa_name = kappa_c
    w = w
  []
[]

[Materials]
  # --------------------------------------------------------------------------
  # gamma(bcc) 自由能 (eV/atom)
  # Redlich–Kister: Chevalier et al., Calphad 28 (2004) 15, Table 3
  # = Mohanty 2011 文献 [21]
  # --------------------------------------------------------------------------
  [free_energy]
    type = DerivativeParsedMaterial
    property_name = f_bulk
    coupled_variables = 'c T'
    constant_names = 'kB JeV L0a L0b L0c L1 L2 L3 L4'
    constant_expressions = '8.617343e-5 96485.33
                            60574.02 -221.45371 24.779079
                            8418.51 512.70 3700.10 5860.56'
    expression = 'kB*T*(c*log(c) + (1-c)*log(1-c))
                  + c*(1-c)*( (L0a + L0b*T + L0c*T*log(T))
                            + L1*(1-2*c)
                            + L2*(1-2*c)^2
                            + L3*(1-2*c)^3
                            + L4*(1-2*c)^4 )/JeV'
    derivative_order = 3
  []

  # --------------------------------------------------------------------------
  # 原子迁移率 beta_i = scale * beta0_i * exp(-H_i/(R*T))   [um^2/(eV*s)]
  #
  # H_U=128 kJ/mol, H_Zr=195 kJ/mol（Mohanty 正文, 引自 Hofman 1996）
  #
  # beta0 由 Mohanty Fig.2 + Fig.4 在 ~1400 K 标定（SI）:
  #   beta_U(1400K)=1.80e8 m^2/(J*s),  beta_Zr(1400K)=2.29e7 m^2/(J*s)
  #   ⇒ beta0_U=1.074e13, beta0_Zr=4.31e14  m^2/(J*s)
  # 换算 1 m^2/(J*s) = 1.60217662e-7 um^2/(eV*s):
  #   beta0_U = 1.7204e6,  beta0_Zr = 6.9045e7  um^2/(eV*s)
  #
  # scale: 绝对时间标定因子（保持 beta_U/beta_Zr 比值与激活能不变）。
  #   仅用 Fig.2/4 的 beta 绝对值时动力学偏快（见 Records.md）；
  #   按 Fig.7 热端 c(t) 标定。当前值使 5/10/30 day 剖面与 Fig.6 同量级。
  # --------------------------------------------------------------------------
  [atomic_mobility]
    type = DerivativeParsedMaterial
    property_name = beta_eff
    coupled_variables = 'c T'
    constant_names = 'scale b0U b0Zr HU HZr R'
    constant_expressions = '0.035 1.7204e6 6.9045e7 128000 195000 8.314462'
    expression = 'scale*(c*b0U*exp(-HU/(R*T)) + (1-c)*b0Zr*exp(-HZr/(R*T)))'
    derivative_order = 2
  []

  # 化学迁移率 M = c(1-c)*beta_eff
  [mobility]
    type = DerivativeParsedMaterial
    property_name = M
    coupled_variables = 'c T'
    material_property_names = 'beta_eff(c,T)'
    expression = 'c*(1-c)*beta_eff'
    derivative_order = 2
    outputs = exodus
  []

  # --------------------------------------------------------------------------
  # Qeff（SoretDiffusion 的 Q_name）[eV]
  #
  # 热迁移迁移率采用与 Mc 相同的 Darken 成分加权（避免 c→0 时 ρ 发散、
  # 冷端被抽干；若按 Mohanty 式(3)字面 [bU*QU - bZr*QZr] 不含 c 权重,
  # 则 ρ=(...)/(c*bU+(1-c)*bZr) 在 c→0 时 → (bU/bZr)*QU+|QZr| 过大）:
  #   M_T = c(1-c)[ c*beta_U*QU - (1-c)*beta_Zr*QZr ]
  #   rho = M_T/M = [c*bU*QU - (1-c)*bZr*QZr]/[c*bU+(1-c)*bZr]
  #
  # MOOSE: J_s = -D*Qeff*c/(kB*T^2)*grad(T), 要 J = +M_T*grad(T)/T
  # ⇒ Qeff = -rho < 0（Zr→热端）
  #
  # Q*_U = +2.5e-20 J/atom = +0.15604 eV（Mohanty 正文 / Fig.3）
  # Q*_Zr= -17.5e-20 J/atom = -1.09226 eV
  # ※ 原文单位是 J/atom，不是 kJ/mol
  # --------------------------------------------------------------------------
  [Qeff_mat]
    type = DerivativeParsedMaterial
    property_name = Qeff
    coupled_variables = 'c T'
    constant_names = 'b0U b0Zr HU HZr R QU QZr'
    constant_expressions = '1.7204e6 6.9045e7 128000 195000 8.314462
                            0.156037728 -1.092264098'
    expression = '-( c*b0U*exp(-HU/(R*T))*QU - (1-c)*b0Zr*exp(-HZr/(R*T))*QZr )
                  / ( c*b0U*exp(-HU/(R*T)) + (1-c)*b0Zr*exp(-HZr/(R*T)) )'
    derivative_order = 1
    outputs = exodus
  []

  # D_soret = M*kB*T/c  ⇒  D_soret*Qeff*c/(kB*T^2) = M*Qeff/T
  [soret_diffusivity]
    type = DerivativeParsedMaterial
    property_name = D_soret
    coupled_variables = 'c T'
    material_property_names = 'M(c,T)'
    constant_names = 'kB'
    constant_expressions = '8.617343e-5'
    expression = 'M*kB*T/c'
    derivative_order = 1
  []

  [constants]
    type = GenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '1e-5'
  []
[]

[Postprocessors]
  [c_hot] # x=0, 热端
    type = PointValue
    variable = c
    point = '0 0 0'
  []
  [c_cold] # x=300, 冷端
    type = PointValue
    variable = c
    point = '300 0 0'
  []
  [c_avg]
    type = ElementAverageValue
    variable = c
  []
  [dt]
    type = TimestepSize
  []
[]

[VectorPostprocessors]
  [profile]
    type = LineValueSampler
    variable = 'c w T'
    start_point = '0 0 0'
    end_point = '300 0 0'
    num_points = 301
    sort_by = x
    execute_on = 'FINAL'
    outputs = profile_csv
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_max_its = 25
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-10
  l_max_its = 50
  l_tol = 1e-5

  start_time = 0
  end_time = 2.592e6 # 30 days

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 500
    growth_factor = 1.4
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  dtmax = 21600 # 0.25 day：保证 5/10/30 day 瞬态分辨率
[]

[Outputs]
  exodus = true
  csv = true
  [profile_csv]
    type = CSV
    execute_on = 'FINAL'
  []
  # 与 Mohanty Fig.6 对齐
  sync_times = '432000 864000 2592000'
  perf_graph = true
[]


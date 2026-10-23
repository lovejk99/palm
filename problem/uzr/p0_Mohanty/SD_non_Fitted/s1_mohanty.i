# ==============================================================================
# Mohanty et al., J. Nucl. Mater. 414 (2011) 211-216, Fig. 6
# Single-phase gamma(bcc) U-Zr thermotransport
#
# Printed equations, not the scaled/weighted variant in ../1:
#   dc/dt = div( M * grad(w) ) - div( M_T * grad(T) / T )
#   w     = df/dc                                          (kappa = 0)
#   M     = c(1-c)[ c*beta_U + (1-c)*beta_Zr ]             Eq. (2)
#   M_T   = c(1-c)[ beta_U*Q_U - beta_Zr*Q_Zr ]            Eq. (3)
#   beta  = beta0 * exp(-H/(R*T))                          Eq. (4), no scale factor
#
# Units: length um, time s, energy eV/atom, temperature K.
# SoretDiffusion hard-codes k_B in eV/K, so f and Q are in eV/atom.
# ==============================================================================

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 300
  xmin = 0
  xmax = 300 # um, Mohanty section 3
  elem_type = EDGE
[]

[Variables]
  [c]
    order = FIRST
    family = LAGRANGE
  []
  [w]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [c_IC]
    # U-39 at.% Zr, Mohanty Fig. 6
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
    # Figs. 5 and 6 draw the hot end at x = 0 and the cold end at x = 300 um.
    # The text states the end temperatures as 1405 K and 1050 K (Eq. (8) prints
    # the opposite assignment of those ends to x = 0). 1405 K is the stated
    # hot-end value; 1400 K is only the figure tick.
    type = FunctionAux
    variable = T
    function = '1405 - (1405 - 1050)/300 * x'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Kernels]
  # 论文式(1) 拆成四个内核。w 方程是质量守恒，c 方程是 w 的定义。
  #
  #   (psi, dc/dt) + (grad psi, M grad w)
  #                + (grad psi, D*Q*c/(kB*T^2) grad T) = 0
  #   (psi, w - df/dc) + (grad psi, kappa grad c) = 0
  #
  # 论文要求的是
  #   dc/dt = div(M grad w) - div(M_T gradT / T)
  #   w = df/dc ,   kappa = 0
  # 因此热迁移项必须满足
  #   D * Q * c / (kB * T^2) = - M_T / T
  # 取 D = M*kB*T/c 、Q = -M_T/M 后，上式成立。
  # 这个 Q 是 SoretDiffusion 的参数，不是论文的 Q*。
  # 论文的 Q* 只出现在下面 Materials/M_T 的式(3) 里。
  [c_dot]
    # (psi, dc/dt)，放在 w 方程上
    type = CoupledTimeDerivative
    variable = w
    v = c
  []
  [w_res]
    # (grad psi, M grad w)。M 是式(2)
    type = SplitCHWRes
    variable = w
    mob_name = M
  []
  [w_res_soret]
    # (grad psi, D_soret*Qeff*c/(kB*T^2) grad T)
    # D_soret = M*kB*T/c ，Qeff = -M_T/M，见 Materials
    type = SoretDiffusion
    variable = w
    c = c
    T = T
    diff_name = D_soret
    Q_name = Qeff
  []
  [c_res]
    # w = df/dc - kappa lap(c)。kappa_c = 0，梯度能项消失
    type = SplitCHParsed
    variable = c
    f_name = f_bulk
    kappa_name = kappa_c
    w = w
    coupled_variables = 'T'
  []
[]

[Materials]
  # --------------------------------------------------------------------------
  # gamma free energy, eV/atom.
  # Mohanty takes f from Thermo-Calc and does not print it.
  # Jung et al., Ann. Nucl. Energy 224 (2025) 111679, Eq. (24) and Appendix A.
  # T_ref = 1010 K. G0^{beta-Zr} = G0^{alpha-Zr} as Appendix (A.4).
  # Appendix A prints "ln(T)". The coefficients are Dinsdale's coefficients of
  # T*ln(T), and the appendix cites that SGTE source, so the T factor is kept.
  # Jung prints only one G0^{gamma-U} polynomial (Dinsdale's branch below
  # 1049 K); it is used over the whole 1050-1405 K bar because no other
  # branch is given.
  # J/mol -> eV/atom by /96485.3321. R/F = k_B, so the mixing term is k_B*T_ref.
  # --------------------------------------------------------------------------
  [free_energy]
    type = DerivativeParsedMaterial
    property_name = f_bulk
    coupled_variables = 'c T'
    constant_names = 'kB JeV Tref
                      aU bU cU dU eU fU
                      aZ bZ cZ dZ eZ fZ gZ
                      L0a L0b L1 L2'
    constant_expressions = '8.617333262145e-5 96485.3321 1010
                            -752.767 131.5381 -27.5152 -8.35595e-3 0.967907e-6 204611
                            4474.461 124.9457 -25.607406 -3.40084e-4 -9.792e-9 25233 -7.6143e-11
                            28369.0639 -14.741714 4992.02639 2609.10711'
    expression = '((1-c)*(aU + bU*T + cU*T*log(T) + dU*T^2 + eU*T^3 + fU/T)
                  + c*(aZ + bZ*T + cZ*T*log(T) + dZ*T^2 + eZ*T^3 + fZ/T + gZ*T^4)
                  + c*(1-c)*((L0a + L0b*Tref) + L1*(1-2*c) + L2*(1-2*c)^2))/JeV
                  + kB*Tref*(c*log(c) + (1-c)*log(1-c))'
    derivative_order = 2
  []

  # beta_i = beta0_i * exp(-H_i/(R*T)), um^2/(eV*s). No multiplicative scale.
  # H from the text (Hofman 1996): 128 and 195 kJ/mol.
  # beta0 is not tabulated. It is the Arrhenius prefactor of Figs. 2 and 4 at
  # 1400 K, using the printed Q* so that beta*Q matches Fig. 4:
  #   beta_U(1400 K)  = 1.98e8 m^2/(J*s)
  #   beta_Zr(1400 K) = 2.65e7 m^2/(J*s)
  # 1 m^2/(J*s) = 1.60217662e-7 um^2/(eV*s).
  [beta_U]
    type = DerivativeParsedMaterial
    property_name = beta_U
    coupled_variables = 'T'
    constant_names = 'b0U HU R'
    constant_expressions = '1.892434e6 128000 8.314462618'
    expression = 'b0U*exp(-HU/(R*T))'
    derivative_order = 1
  []
  [beta_Zr]
    type = DerivativeParsedMaterial
    property_name = beta_Zr
    coupled_variables = 'T'
    constant_names = 'b0Z HZ R'
    constant_expressions = '8.004853e7 195000 8.314462618'
    expression = 'b0Z*exp(-HZ/(R*T))'
    derivative_order = 1
  []

  # Eq. (2)
  [mobility]
    type = DerivativeParsedMaterial
    property_name = M
    coupled_variables = 'c T'
    material_property_names = 'beta_U(T) beta_Zr(T)'
    expression = 'c*(1-c)*(c*beta_U + (1-c)*beta_Zr)'
    derivative_order = 1
  []

  # --------------------------------------------------------------------------
  # 式(3)，按原文，括号内没有 c 或 (1-c)：
  #   M_T = c(1-c)[ beta_U*Q_U - beta_Zr*Q_Zr ]
  #
  # 论文给出的是能量，单位 J，正文写明是每个原子携带的热量：
  #   Q_U  = +2.5e-20 J
  #   Q_Zr = -17.5e-20 J
  # eV 同样是每个原子的能量。1 eV = 1.60217662e-19 J，所以只除以元电荷：
  #   Q_U  =  2.5e-20 / 1.60217662e-19 = +0.156037728 eV
  #   Q_Zr = -17.5e-20 / 1.60217662e-19 = -1.092264098 eV
  # 不能再除以 96485。96485 J/mol = 1 eV/atom，用来把 J/mol 换成 eV；
  # 这里的 10^{-20} J 已经是一个原子的能量。
  # beta 用 um^2/(eV s)，Q 用 eV，乘积才是 um^2/s，与 Fig. 4 一致。
  # Q_Zr 代入带负号的原文数值，减号只在式(3) 里出现一次。
  # --------------------------------------------------------------------------
  [M_T]
    type = DerivativeParsedMaterial
    property_name = M_T
    coupled_variables = 'c T'
    material_property_names = 'beta_U(T) beta_Zr(T)'
    constant_names = 'QU QZ'
    constant_expressions = '0.156037728 -1.092264098'
    expression = 'c*(1-c)*(beta_U*QU - beta_Zr*QZ)'
    derivative_order = 1
  []

  # 不是论文里的式子。由 D*Q*c/(kB*T^2) = -M_T/T 和 D = M*kB*T/c
  # 得到 Qeff = -M_T/M。c(1-c) 约掉，分母只剩式(2) 的迁移率组合。
  [Qeff_mat]
    type = DerivativeParsedMaterial
    property_name = Qeff
    coupled_variables = 'c T'
    material_property_names = 'M_T(c,T) M(c,T)'
    expression = '-M_T/M'
    derivative_order = 1
  []

  # D*Q*c/(kB*T^2) = M*Q/T = -M_T/T
  [soret_diffusivity]
    type = DerivativeParsedMaterial
    property_name = D_soret
    coupled_variables = 'c T'
    material_property_names = 'M(c,T)'
    constant_names = 'kB'
    constant_expressions = '8.617333262145e-5'
    expression = 'M*kB*T/c'
    derivative_order = 1
  []

  # Mohanty: the gradient-energy coefficient vanishes in a single-phase alloy.
  [constants]
    type = GenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '0'
  []
[]

[Postprocessors]
  [c_hot]
    type = PointValue
    variable = c
    point = '0 0 0'
  []
  [c_cold]
    type = PointValue
    variable = c
    point = '300 0 0'
  []
  [c_avg]
    type = ElementAverageValue
    variable = c
  []
  [c_min]
    type = ElementExtremeValue
    variable = c
    value_type = min
  []
  [c_max]
    type = ElementExtremeValue
    variable = c
    value_type = max
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
  # SoretDiffusion does not differentiate D(c,T) or Q(c,T).
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = bt

  nl_max_its = 30
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_max_its = 50
  l_tol = 1e-6

  start_time = 0
  end_time = 2592000 # 30 days

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.3
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  # Unscaled Fig. 2 mobilities: hot-end diffusion time is under a day.
  dtmax = 3600
[]

[Outputs]
  exodus = true
  csv = false
  sync_times = '432000 864000 1728000 2592000'
  print_linear_residuals = false
[]

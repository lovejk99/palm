# ==============================================================================
# Mohanty et al., J. Nucl. Mater. 414 (2011) 211-216, Fig. 6
# Single-phase gamma(bcc) U-Zr thermotransport, no fitted parameters
#
# ADThermalMigration reads M_T directly. SoretDiffusion is not used.
# Both beta_U and beta_Zr depend on temperature, as in Fig. 5c and Fig. 6.
# Figs. 5a and 5b hold beta constant and are not this case.
#
# Printed balance (Mohanty Eq. (1), with kappa = 0):
#   dc/dt = div( M * grad(w) ) - div( M_T * grad(T) / T )
#   w     = df/dc
#   M     = c(1-c)[ c*beta_U + (1-c)*beta_Zr ]          Eq. (2)
#   M_T   = c(1-c)[ beta_U*Q_U - beta_Zr*Q_Zr ]         Eq. (3)
#   beta  = beta0 * exp(-H/(R*T))                       Eq. (4)
#
# Weak form solved here, after integration by parts and with zero boundary flux:
#   (psi, dc/dt) + (grad psi, M grad w) - (grad psi, M_T/T grad T) = 0
#   (psi, df/dc - w) + (grad psi, kappa grad c) = 0
#
# Units: length um, time s, energy eV/atom, temperature K.
# ==============================================================================

[Mesh]
  # Mohanty section 3 gives the domain length, L = 300 um, and a control-volume
  # discretisation. It does not give the number of control volumes.
  # nx = 300 places one linear element on each micrometre, so dx = 1 um.
  # That choice is numerical, not a value read from the paper.
  type = GeneratedMesh
  dim = 1
  nx = 300
  xmin = 0
  xmax = 300
  elem_type = EDGE
[]

[Variables]
  [c]
    # Zr atom fraction. Mohanty uses c(x,t) for this field.
    order = FIRST
    family = LAGRANGE
  []
  [w]
    # Chemical potential df/dc, eV/atom. Split unknown; not a paper symbol.
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [c_IC]
    # Fig. 6 initial alloy: U-39 at.% Zr, uniform.
    type = ConstantIC
    variable = c
    value = 0.39
  []
[]

[AuxVariables]
  [T]
    # Prescribed temperature. Not a nonlinear unknown in this step.
  []
[]

[AuxKernels]
  [T_profile]
    # Figs. 5 and 6 draw the hot end at x = 0 and the cold end at x = 300 um.
    # The text states the two end temperatures as 1405 K and 1050 K.
    # Eq. (8) prints the opposite assignment, T(0) = T_min and T(L) = T_max.
    # The figure being reproduced fixes the direction; 1405 K is the stated
    # hot-end value. The 1400 K mark is only the figure tick.
    # Linear T follows from Laplace's equation with constant conductivity.'
    type = FunctionAux
    variable = T
    function = '1405 - (1405 - 1050)/300 * x'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

# No [BCs] block.
# Eq. (6) requires dc/dx = 0 and d^3 c/dx^3 = 0 at both ends.
# The weak form above drops the boundary integral, so the natural condition
# is zero total flux, M dw/dx - M_T (1/T) dT/dx = 0, at x = 0 and x = 300 um.
# A NeumannBC of value 0 would add that same condition a second time.
# d^3 c/dx^3 = 0 is the natural condition of the gradient-energy term. With
# kappa_c = 0 that term is absent, so the condition does not need a BC object.
# Temperature is prescribed by T_profile. It is not solved from Eq. (8).
[Kernels]
  # Four kernels, one per term of the split system. The w-equation carries
  # mass conservation. The c-equation defines w.
  [c_dot]
    # (psi, dc/dt) on the w-equation.
    #   variable = w   residual is tested with the w shape function
    #   v        = c   the time derivative is taken of c, not of w
    type = ADCoupledTimeDerivative
    variable = w
    v = c
  []
  [w_res]
    # (grad psi, M grad w). M is Eq. (2), an AD material named M.
    #   variable = w
    #   mob_name = M
    type = ADSplitCHWRes
    variable = w
    mob_name = M
  []
  [w_res_tm]
    # ADThermalMigration residual:
    #   R = - int (M_T / T) grad(T) . grad(psi) dV
    # which is the weak form of -div(M_T grad(T) / T).
    # With M_T from Eq. (3) positive, Zr moves toward the hot end.
    #   variable         = w          same equation as c_dot and w_res
    #   T                = T          auxiliary temperature
    #   thermal_mobility = M_T       AD property; the kernel does not form
    #                                 D or Q, and it does not contain k_B
    type = ADThermalMigration
    variable = w
    T = T
    thermal_mobility = M_T
  []
  [c_res]
    # (psi, df/dc - w) + (grad psi, kappa grad c) = 0.
    # kappa_c = 0, so this is w = df/dc.
    #   variable   = c
    #   f_name     = f_bulk     AD free-energy property; the kernel uses df/dc
    #   kappa_name = kappa_c
    #   w          = w
    #   coupled_variables = T   f depends on T; T is auxiliary, so this adds
    #                           no nonlinear unknown
    type = ADSplitCHParsed
    variable = c
    f_name = f_bulk
    kappa_name = kappa_c
    w = w
    coupled_variables = 'T'
  []
[]

[Materials]
  # All properties consumed by AD kernels are AD properties, so the Newton
  # Jacobian includes df/dc, M(c,T) and M_T(c,T).

  # --------------------------------------------------------------------------
  # gamma free energy, eV/atom.
  # Mohanty takes f from Thermo-Calc and does not print the function.
  # Adopted source: Jung et al., Ann. Nucl. Energy 224 (2025) 111679,
  # Eq. (24) and Appendix A. Not fitted to Fig. 6.
  #
  #   f = (1-c) G0^{gamma-U}(T) + c G0^{beta-Zr}(T)
  #       + R T_ref [c ln c + (1-c) ln(1-c)]
  #       + c(1-c)[ L0 + L1(1-2c) + L2(1-2c)^2 ]
  #   T_ref = 1010 K
  #   L0 = 28369.0639 - 14.741714 T_ref
  #   L1 = 4992.02639
  #   L2 = 2609.10711
  #   G0^{beta-Zr} = G0^{alpha-Zr}          Appendix (A.4)
  #
  # Appendix A prints "ln(T)". Those coefficients are Dinsdale's coefficients'
  # of T ln(T), and the appendix cites that SGTE source, so T*ln(T) is used.
  # Jung prints only one G0^{gamma-U} polynomial (the branch valid below
  # 1049 K). It is used on the whole 1050-1405 K bar; no second branch is given.
  #
  # J/mol is divided by 96485.3321 to reach eV/atom. R/F = k_B, so the mixing
  # term is k_B*T_ref*(...).
  # --------------------------------------------------------------------------
  [free_energy]
    type = ADDerivativeParsedMaterial
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

  # Eq. (4), both species, temperature dependent. This is Fig. 5c / Fig. 6,
  # not the constant-mobility curves in Figs. 5a and 5b. Q* stays constant.
  #
  #   beta_i = beta0_i * exp(-H_i/(R*T))
  #   R = 8.314462618 J/(mol*K)
  #   H_U = 128000 J/mol, H_Zr = 195000 J/mol
  #
  # Prefactors fitted to Fig. 2, used here without a further scale factor:
  #   beta0_U  = 119945.49892e8  = 1.1994549892e13
  #   beta0_Zr = 5194864.7e8     = 5.1948647e14
  #
  # Fig. 2 prints the unit as "m^2·sec/J", which as written is m^2*s/J.
  # Fig. 4 prints beta*Q* in m^2/s. With Q* in J, that product has units
  # m^2/s only when beta is m^2/(J*s). Eq. (1) needs the same dimension:
  # beta * grad(w) must be a velocity when w is an energy. The reduced time
  # in the paper, proportional to beta * t / L^2, is dimensionless only for
  # m^2/(J*s). The input therefore takes the Fig. 2 magnitudes as m^2/(J*s).
  #
  # 1 m^2/(J*s) = 1.60217662e-7 um^2/(eV*s):
  #   1 m^2 = 1e12 um^2, and 1 J = 1/1.60217662e-19 eV,
  #   so multiply by 1e12 * 1.60217662e-19.
  #   beta0_U  = 1.92173874e6 um^2/(eV*s)
  #   beta0_Zr = 8.32309077e7 um^2/(eV*s)
  # At 1400 K this gives beta_U*Q_U = 5.03e-12 m^2/s and
  # beta_Zr*Q_Zr = -4.82e-12 m^2/s, on the scale of Fig. 4.
  [beta_U]
    type = ADDerivativeParsedMaterial
    property_name = beta_U
    coupled_variables = 'T'
    constant_names = 'b0U HU R'
    constant_expressions = '1.92173874e6 128000 8.314462618'
    expression = 'b0U*exp(-HU/(R*T))'
    derivative_order = 2
  []
  [beta_Zr]
    type = ADDerivativeParsedMaterial
    property_name = beta_Zr
    coupled_variables = 'T'
    constant_names = 'b0Z HZ R'
    constant_expressions = '8.32309077e7 195000 8.314462618'
    expression = 'b0Z*exp(-HZ/(R*T))'
    derivative_order = 2
  []

  # Eq. (2). Chemical mobility passed to ADSplitCHWRes as mob_name = M.
  [mobility]
    type = ADDerivativeParsedMaterial
    property_name = M
    coupled_variables = 'c T'
    material_property_names = 'beta_U(T) beta_Zr(T)'
    expression = 'c*(1-c)*(c*beta_U + (1-c)*beta_Zr)'
    derivative_order = 2
  []

  # --------------------------------------------------------------------------
  # Eq. (3), as printed. No c or (1-c) inside the bracket, and no scale.
  #   M_T = c(1-c)[ beta_U*Q_U - beta_Zr*Q_Zr ]
  #
  # The text gives the heats per atom, in joule, not in kJ/mol:
  #   Q_U  = +2.5e-20 J
  #   Q_Zr = -17.5e-20 J
  # 1 eV = 1.60217662e-19 J, so
  #   Q_U  = +0.156037728 eV
  #   Q_Zr = -1.092264098 eV
  # 96485 J/mol = 1 eV/atom converts molar energies. These 10^{-20} J values
  # are already per atom, so they are not divided by 96485.
  # Q_Zr is stored with its printed negative sign. The minus in Eq. (3)
  # appears once.
  # Passed to ADThermalMigration as thermal_mobility = M_T.
  # --------------------------------------------------------------------------
  [M_T]
    type = ADDerivativeParsedMaterial
    property_name = M_T
    coupled_variables = 'c T'
    material_property_names = 'beta_U(T) beta_Zr(T)'
    constant_names = 'QU QZ'
    constant_expressions = '0.156037728 -1.092264098'
    expression = 'c*(1-c)*(beta_U*QU - beta_Zr*QZ)'
    derivative_order = 2
  []

  # Mohanty: the gradient-energy coefficient vanishes in a single-phase alloy.
  [constants]
    type = ADGenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '0'
  []
[]

[Postprocessors]
  [c_hot]
    # x = 0, hot end
    type = PointValue
    variable = c
    point = '0 0 0'
  []
  [c_cold]
    # x = 300 um, cold end
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
  # AD kernels supply the Jacobian of M(c,T), M_T(c,T) and df/dc.
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = bt

  nl_max_its = 30
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_max_its = 50
  l_tol = 1e-6

  start_time = 0
  end_time = 2592000 # 30 days = 30*86400 s. Figs. 6 and 7.

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.3
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  # Fig. 2 mobilities, unscaled. The hot-end diffusion time is under one day,
  # so the step is capped at 1 hour to resolve the early transient.
  dtmax = 3600
[]

[Outputs]
  exodus = true
  csv = false
  # 5, 10, 20 and 30 days, for comparison with Fig. 6.
  sync_times = '432000 864000 1728000 2592000'
  print_linear_residuals = false
[]

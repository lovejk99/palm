# Equivalence of ADThermalMigration with the built-in SoretDiffusion under the
# constant-coefficient mapping used in palm/problem/uzr/s1_mohanty/calibrated.i:
#
#   SoretDiffusion residual :  + D_soret * Qeff * c / (kB T^2) grad T . grad psi
#   ADThermalMigration      :  - (M_T / T) grad T . grad psi
#
# They coincide for  D_soret = M kB T / c  and  Qeff = - M_T / M.
# Default `active` list uses ADThermalMigration; the test spec re-runs the file with
#   Kernels/active='c_dot w_res w_res_soret c_res'
# and compares both CSV outputs against the same gold file.
#
# Units: um, s, eV, K (kB = 8.617343e-5 eV/K is hard-coded inside SoretDiffusion).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 40
  xmin = 0
  xmax = 100
[]

[Variables]
  [c]
    initial_condition = 0.39
  []
  [w]
  []
[]

[AuxVariables]
  [T]
  []
[]

[AuxKernels]
  [T_aux]
    type = FunctionAux
    variable = T
    function = '1400 - 3.5*x'
    execute_on = 'INITIAL'
  []
[]

[Kernels]
  active = 'c_dot w_res w_res_tm c_res'
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
  [w_res_tm]
    type = ADThermalMigration
    variable = w
    T = T
    thermal_mobility = MT
  []
  [w_res_soret]
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
    f_name = f
    kappa_name = kappa
    w = w
  []
[]

[Materials]
  # regular solution, eV/atom
  [f]
    type = DerivativeParsedMaterial
    property_name = f
    coupled_variables = 'c T'
    constant_names = 'kB L0'
    constant_expressions = '8.617343e-5 0.05'
    expression = 'kB*T*(c*log(c) + (1-c)*log(1-c)) + L0*c*(1-c)'
    derivative_order = 2
  []
  # constant atomic mobility beta and constant beta*Q product B, so that
  #   M   = c(1-c) beta,   M_T = c(1-c) B,   Qeff = -B/beta,   D_soret = (1-c) beta kB T
  [M]
    type = DerivativeParsedMaterial
    property_name = M
    coupled_variables = 'c'
    constant_names = 'beta'
    constant_expressions = '2.0'
    expression = 'c*(1-c)*beta'
    derivative_order = 2
  []
  [MT_ad]
    type = ADParsedMaterial
    property_name = MT
    coupled_variables = 'c'
    constant_names = 'B'
    constant_expressions = '0.8'
    expression = 'c*(1-c)*B'
  []
  [Qeff]
    type = GenericConstantMaterial
    prop_names = 'Qeff kappa'
    prop_values = '-0.4 0'
  []
  [D_soret]
    type = DerivativeParsedMaterial
    property_name = D_soret
    coupled_variables = 'c T'
    constant_names = 'beta kB'
    constant_expressions = '2.0 8.617343e-5'
    expression = '(1-c)*beta*kB*T'
    derivative_order = 1
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
    point = '100 0 0'
  []
  [c_int]
    type = ElementIntegralVariablePostprocessor
    variable = c
  []
[]

[Preconditioning]
  [smp]
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
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  dt = 50
  num_steps = 10
[]

[Outputs]
  csv = true
[]

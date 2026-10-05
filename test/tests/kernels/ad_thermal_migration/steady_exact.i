# Exact-solution test for ADThermalMigration with constant M_T.
#
#   0 = div( D grad c ) - div( M_T grad T / T ),   T = T0 + g x  (prescribed, linear)
#
# Left: Dirichlet c(0) = c0.  Right: natural BC of the weak form = zero *total* flux,
# which forces D c' - M_T T'/T = 0 everywhere, so
#
#   c(x) = c0 + (M_T / D) ln( T(x) / T0 ).
#
# With D = 1, M_T = 0.5, T0 = 1000, g = 500:  c(1) = 1 + 0.5 ln 1.5 = 1.202733.
# The test checks the L2 error against the exact solution and the right-end value.
# Running with  AuxKernels/T_aux/function='1500-500*x'  (reversed gradient) gives
# c(1) = 1 + 0.5 ln(2/3) = 0.797267, i.e. the migration direction reverses.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 50
  xmin = 0
  xmax = 1
[]

[Variables]
  [c]
    initial_condition = 1
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
    function = '1000 + 500*x'
    execute_on = 'INITIAL'
  []
[]

[Functions]
  [c_exact]
    type = ParsedFunction
    # c0 + (MT/D) * ln(T/T0); T itself is re-evaluated here so the reversed case
    # must also override this function (see tests spec).
    expression = '1 + 0.5*log((1000 + 500*x)/1000)'
  []
[]

[Kernels]
  [diff]
    type = ADMatDiffusion
    variable = c
    diffusivity = D
  []
  [thermal_migration]
    type = ADThermalMigration
    variable = c
    T = T
    thermal_mobility = MT
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = c
    boundary = left
    value = 1
  []
[]

[Materials]
  [props]
    type = ADGenericConstantMaterial
    prop_names = 'D MT'
    prop_values = '1 0.5'
  []
[]

[Postprocessors]
  [l2_error]
    type = ElementL2Error
    variable = c
    function = c_exact
  []
  [c_right]
    type = PointValue
    variable = c
    point = '1 0 0'
  []
  [c_left]
    type = PointValue
    variable = c
    point = '0 0 0'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-12
[]

[Outputs]
  csv = true
[]

# Jacobian test for ADThermalMigration with M_T = M_T(c, T) and T a nonlinear variable.
# Both the c-dependence (through the parsed material) and the T-dependence
# (through 1/T, grad T and the material) must be captured by automatic differentiation.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 3
[]

[Variables]
  [c]
    [InitialCondition]
      type = RandomIC
      min = 0.2
      max = 0.6
    []
  []
  [T]
    [InitialCondition]
      type = RandomIC
      min = 1000
      max = 1400
    []
  []
[]

[Kernels]
  [c_time]
    type = ADTimeDerivative
    variable = c
  []
  [c_diff]
    type = ADMatDiffusion
    variable = c
    diffusivity = D
  []
  [c_tm]
    type = ADThermalMigration
    variable = c
    T = T
    thermal_mobility = MT
  []
  [T_diff]
    type = ADDiffusion
    variable = T
  []
[]

[Materials]
  [D]
    type = ADGenericConstantMaterial
    prop_names = 'D'
    prop_values = '1'
  []
  [MT]
    type = ADParsedMaterial
    property_name = MT
    coupled_variables = 'c T'
    constant_names = 'b0U b0Zr HU HZr R QU QZr'
    constant_expressions = '1.7204e6 6.9045e7 128000 195000 8.314462 0.156037728 -1.092264098'
    expression = 'c*(1-c)*(b0U*exp(-HU/(R*T))*QU - b0Zr*exp(-HZr/(R*T))*QZr)'
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
  solve_type = NEWTON
  num_steps = 1
  dt = 1
[]

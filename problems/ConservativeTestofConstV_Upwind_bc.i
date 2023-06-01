[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
    xmin = 0.0
    ymin = 0.0
    xmax = 1.0
    ymax = 1.0
  []
[]

[Variables]
  [rho]
  []
[]

[ICs]
  [rho_ic]
    type = FunctionIC
    variable = rho
    function = 'if(x<0.5,1,0)'
  []
[]

[Kernels]
  [Time_Deri]
    type = MassLumpedTimeDerivative
    variable = rho
  []
  [Flux]
    type = ConservativeAdvection
    variable = rho
    velocity = '1.0 0.0 0.0'
    upwinding_type = full
  []
[]

[Materials]
[]

[BCs]
  [right]
    type = DirichletBC
    variable = rho
    boundary = right
    value = 1.0
  []
[]

# Transient (time-dependent) details for simulations go here:
[Executioner]
  type = Transient   # Here we use the Transient Executioner (instead of steady)
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg          31'
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_tol = 1e-8

  start_time = 0.0
  end_time = 1.0
  dt = 0.01
  dtmin = 0.001
[]

[VectorPostprocessors]
  [rho]
    type = LineValueSampler
    variable = rho
    start_point = '0 0.5 0'
    end_point = '1 0.5 0'
    num_points = 41
    sort_by = x
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    file_base = rho_x_out_upwind
    execute_on = final
  []
[]

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
    initial_condition = 8.e3
  []
[]

[Kernels]
  [Time_Deri]
    type = TestTimeDerivative
    variable = rho
  []
  [Flux]
    type = ConservativeAdvection
    variable = rho
    velocity = '0.75 0.0 0.0'
    upwinding_type = none
  []
  [Diffusion]
    type = CoefDiffusion
    variable = rho
  use_displaced_mesh = false
    coef = 0.01
  []
[]

[Materials]
[]

[BCs]
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
    file_base = rho_x_out_adm
    execute_on = final
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
  []
  [block1]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '0.01 1 0'
    input = gen
  []
  [block2]
    type = SubdomainBoundingBoxGenerator
    block_id = 2
    bottom_left = '0.01 0 0'
    top_right = '0.99 1 0'
    input = block1
  []
  [block3]
    type = SubdomainBoundingBoxGenerator
    block_id = 3
    bottom_left = '0.99 0 0'
    top_right = '1 1 0'
    input = block2
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
    type = ConservativeAdvectionSchmid
    variable = rho
    upwinding_type = none
      dislo_sign = positive
      slip_sys_index = 0
  []
  [Diffusion]
    type = CoefDiffusion
    variable = rho
  use_displaced_mesh = false
    coef = 0.01
  []
[]

[Materials]
  [leftBC]
    type = BCMaterial
    block = '1'
    nss = 1
  []
  [vel]
    type = VelocityMaterial
    block = '2'
    nss = 1
  []
  [rightBC]
    type = BCMaterial
    block = '3'
    nss = 1
  []
[]

[BCs]
  [right]
    type = DirichletBC
    variable = rho
    boundary = right
    value = 8.e3
  []
  [left]
    type = DirichletBC
    variable = rho
    boundary = left
    value = 8.e3
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
  end_time = 10.0
  dt = 0.05
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
    file_base = rho_x_out_adm_bc_none_coef_1e-3
    execute_on = final
  []
[]

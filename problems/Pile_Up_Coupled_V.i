[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 200
    xmin = 0.0
    xmax = 1
    ymin = 0.0
    ymax = 1
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
  [rhoep]
    initial_condition = 8.e3
  []
  [rhoen]
    initial_condition = 8.e3
  []
[]

[Kernels]
  [Edeg_Pos_Time_Deri]
    type = MassLumpedTimeDerivative
    variable = rhoep
  []
  [Edge_Pos_Flux]
    type = ConservativeAdvectionSchmid
    variable = rhoep
    upwinding_type = full
      dislo_sign = positive
      slip_sys_index = 0
  []
  [Edeg_Neg_Time_Deri]
    type = MassLumpedTimeDerivative
    variable = rhoen
  []
  [Edge_Neg_Flux]
    type = ConservativeAdvectionSchmid
    variable = rhoen
    upwinding_type = full
      dislo_sign = negative
      slip_sys_index = 0
  []
[]

[Materials]
  [leftBC]
    type = DisloVelocityBC
    block = '1'
    nss = 1
    rhoen = rhoen
    rhoep = rhoep
  []
  [vel]
    type = DisloVelocityCoupled
    block = '2'
    nss = 1
    rhoen = rhoen
    rhoep = rhoep
  []
  [rightBC]
    type = DisloVelocityBC
    block = '3'
    nss = 1
    rhoen = rhoen
    rhoep = rhoep
  []
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
  end_time = 0.001
  dt = 5.e-7
  dtmin = 1.e-9
[]

[VectorPostprocessors]
  [rhoep]
    type = LineValueSampler
    variable = rhoep
    start_point = '0 0.05 0'
    end_point = '0.1 0.05 0'
    num_points = 101
    sort_by = x
  []
  [rhoen]
    type = LineValueSampler
    variable = rhoen
    start_point = '0 0.05 0'
    end_point = '0.1 0.05 0'
    num_points = 101
    sort_by = x
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    file_base = rhoe_x_out_l1_impbc
    execute_on = final
  []
[]

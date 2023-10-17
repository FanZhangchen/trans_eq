[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 50
    xmin = 0.0
    ymin = 0.0
    xmax = 0.1
    ymax = 0.1    
  []
[]

[Variables]
  [rhoe1]
    initial_condition = 8.e3
  []
  [rhoe2]
    initial_condition = 8.e3
  []
  [rhoe3]
    initial_condition = 8.e3
  []
  [rhoe4]
    initial_condition = 8.e3
  []
  [rhos1]
    initial_condition = 8.e3
  []
  [rhos2]
    initial_condition = 8.e3
  []
  [rhos3]
    initial_condition = 8.e3
  []
  [rhos4]
    initial_condition = 8.e3
  []
[]

[Kernels]

  [Edeg_Time_Deri_1]
    type = MassLumpedTimeDerivative
    variable = rhoe1
  []
  [Edge_Flux_1]
    type = ConservativeAdvectionSchmid2
    variable = rhoe1
    upwinding_type = full
      dislo_character = edge
      dislo_sign = positive
      slip_sys_index = 0
  []
  [Screw_Time_Deri_1]
    type = MassLumpedTimeDerivative
    variable = rhos1
  []
  [Screw_Flux_1]
    type = ConservativeAdvectionSchmid2
    variable = rhos1
    upwinding_type = full
      dislo_character = screw
      dislo_sign = positive
      slip_sys_index = 0
  []

  [Edeg_Time_Deri_2]
    type = MassLumpedTimeDerivative
    variable = rhoe2
  []
  [Edge_Flux_2]
    type = ConservativeAdvectionSchmid2
    variable = rhoe2
    upwinding_type = full
      dislo_character = edge
      dislo_sign = negative
      slip_sys_index = 0
  []
  [Screw_Time_Deri_2]
    type = MassLumpedTimeDerivative
    variable = rhos2
  []
  [Screw_Flux_2]
    type = ConservativeAdvectionSchmid2
    variable = rhos2
    upwinding_type = full
      dislo_character = screw
      dislo_sign = positive
      slip_sys_index = 0
  []

  [Edeg_Time_Deri_3]
    type = MassLumpedTimeDerivative
    variable = rhoe3
  []
  [Edge_Flux_3]
    type = ConservativeAdvectionSchmid2
    variable = rhoe3
    upwinding_type = full
      dislo_character = edge
      dislo_sign = negative
      slip_sys_index = 0
  []
  [Screw_Time_Deri_3]
    type = MassLumpedTimeDerivative
    variable = rhos3
  []
  [Screw_Flux_3]
    type = ConservativeAdvectionSchmid2
    variable = rhos3
    upwinding_type = full
      dislo_character = screw
      dislo_sign = negative
      slip_sys_index = 0
  []

  [Edeg_Time_Deri_4]
    type = MassLumpedTimeDerivative
    variable = rhoe4
  []
  [Edge_Flux_4]
    type = ConservativeAdvectionSchmid2
    variable = rhoe4
    upwinding_type = full
      dislo_character = edge
      dislo_sign = positive
      slip_sys_index = 0
  []
  [Screw_Time_Deri_4]
    type = MassLumpedTimeDerivative
    variable = rhos4
  []
  [Screw_Flux_4]
    type = ConservativeAdvectionSchmid2
    variable = rhos4
    upwinding_type = full
      dislo_character = screw
      dislo_sign = negative
      slip_sys_index = 0
  []
[]


[Materials]
  [vel]
    type = DisloVelocityCompleted
    nss = 1
    rhoe1 = rhoe1
    rhoe2 = rhoe2
    rhoe3 = rhoe3
    rhoe4 = rhoe4
    rhos1 = rhos1
    rhos2 = rhos2
    rhos3 = rhos3
    rhos4 = rhos4
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
  end_time = 0.2
  dt = 2.e-6
  dtmin = 1.e-9
[]

[VectorPostprocessors]
  [rhoe1]
    type = LineValueSampler
    variable = rhoe1
    start_point = '0 0.05 0'
    end_point = '0.1 0.05 0'
    num_points = 41
    sort_by = x
  []
  [rhoe2]
    type = LineValueSampler
    variable = rhoe2
    start_point = '0 0.05 0'
    end_point = '0.1 0.05 0'
    num_points = 41
    sort_by = x
  []
  [rhoe3]
    type = LineValueSampler
    variable = rhoe3
    start_point = '0 0.05 0'
    end_point = '0.1 0.05 0'
    num_points = 41
    sort_by = x
  []
  [rhoe4]
    type = LineValueSampler
    variable = rhoe4
    start_point = '0 0.05 0'
    end_point = '0.1 0.05 0'
    num_points = 41
    sort_by = x
  []
  [rhos1]
    type = LineValueSampler
    variable = rhos1
    start_point = '0.05 0 0'
    end_point = '0.05 0.1 0'
    num_points = 41
    sort_by = y
  []
  [rhos2]
    type = LineValueSampler
    variable = rhos2
    start_point = '0.05 0 0'
    end_point = '0.05 0.1 0'
    num_points = 41
    sort_by = y
  []
  [rhos3]
    type = LineValueSampler
    variable = rhos3
    start_point = '0.05 0 0'
    end_point = '0.05 0.1 0'
    num_points = 41
    sort_by = y
  []
  [rhos4]
    type = LineValueSampler
    variable = rhos4
    start_point = '0.05 0 0'
    end_point = '0.05 0.1 0'
    num_points = 41
    sort_by = y
  []
[]

[Outputs]
  exodus = true
  interval = 100
  [csv]
    type = CSV
    file_base = rho_x_out_l1e-1_2d_8dofs
    execute_on = final
  []
[]

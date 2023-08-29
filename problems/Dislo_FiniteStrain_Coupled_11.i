[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 10
  xmin = 0.0
  ymin = 0.0
  xmax = 0.1
  ymax = 0.1  
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
      family = LAGRANGE
  []
  [disp_y]
    order = FIRST
      family = LAGRANGE
  []
  [rho_edge_pos_1]
    initial_condition = 8.e3
  []
  [rho_edge_neg_1]
    initial_condition = 8.e3
  []
[]

[AuxVariables]
  [rhot]
    order = CONSTANT
      family = MONOMIAL_VEC
  []
  [rhognd]
    order = CONSTANT
      family = MONOMIAL_VEC
  []

  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./e_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [slip_rate]
    order = CONSTANT
    family = MONOMIAL
  []

  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []

  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []

  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []

  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []

[]

[Functions]
  [disp_load]
    type = ParsedFunction
    value = '0.05*t'
  []
[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFile
    prop_file_name = 'euler_ang_test.inp'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    ngrain = 1
    read_type = indexgrain
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    use_displaced_mesh = true
    add_variables = true
    generate_output = 'deformation_gradient_xx deformation_gradient_xy deformation_gradient_yy'
  [../]
  [Edeg_Pos_Time_Deri]
    type = MassLumpedTimeDerivative
    variable = rho_edge_pos_1
  []
  [Edge_Pos_Flux]
    type = ConservativeAdvectionSchmid
    variable = rho_edge_pos_1
    upwinding_type = full
      dislo_sign = positive
      slip_sys_index = 0
  []
  [Edeg_Neg_Time_Deri]
    type = MassLumpedTimeDerivative
    variable = rho_edge_neg_1
  []
  [Edge_Neg_Flux]
    type = ConservativeAdvectionSchmid
    variable = rho_edge_neg_1
    upwinding_type = full
      dislo_sign = negative
      slip_sys_index = 0
  []
[]

[AuxKernels]
  [rhot]
    type = TotalDislocationDensity
    variable = rhot
    execute_on = timestep_end
    rhoep = rho_edge_pos_1
    rhoen = rho_edge_neg_1
  []
  [rhognd]
    type = GNDDislocationDensity
    variable = rhognd
    execute_on = timestep_end
    rhoep = rho_edge_pos_1
    rhoen = rho_edge_neg_1
  []

  [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]

 [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]

  [./e_xy]
    type = RankTwoAux
    variable = e_xy
    rank_two_tensor = lage
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]

  [slip_rate]
    type = MaterialStdVectorAux
    variable = slip_rate
    property = slip_rate
    index = 0
    execute_on = timestep_end
  []

  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
  [../]

[]

[Materials]
  [./cp]
    type = FiniteStrainCrystalPlasticityDislo
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt # no need to normalize vectors
    nss = 1 # Number of slip systems
    num_slip_sys_flowrate_props = 2 #Number of flow rate properties in a slip system
    flowprops = '1 1 0.001 0.01' # slip rate equations parameters
  hprops = '1.0 3629.0 216.0 300.5 2.5' # hardening properties
    gprops = '1 1 0.44' # initial values of slip system resistances (start_slip_sys, end_slip_sys, value)
    tan_mod_type = exact
  rho_edge_pos_1 = rho_edge_pos_1
  rho_edge_neg_1 = rho_edge_neg_1
    maxiter = 150
    maximum_substep_iteration = 10
  [../]
  [elastic_tensor]
    type = ComputeElasticityTensorCP
     C_ijkl = '2.50e5 2.00e5 2.00e5 2.50e5 2.00e5 2.50e5 1.00e5 1.00e5 1.00e5'
    fill_method = symmetric9
  []
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
[]

[BCs]
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0.0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  []
  [top_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'top'
    function = disp_load
  []
  [top_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 0.0 
  []

  [./Periodic]

    [./auto_rho_edge_pos_boundary_x]
      variable = rho_edge_pos_1
      primary = 'left'
    secondary = 'right'
    translation = '0.1 0.0 0.0'
    [../]

    [./auto_rho_edge_neg_boundary_x]
      variable = rho_edge_neg_1
      primary = 'left'
    secondary = 'right'
    translation = '0.1 0.0 0.0'
    [../]

  [../]

[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

# Transient (time-dependent) details for simulations go here:
[Executioner]

  type = Transient
  solve_type = 'PJFNK'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg          31'
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_tol = 1e-8

  start_time = 0.0
  end_time = 0.1 #0.01
  dt = 5.e-6
  dtmin = 1.e-9
[]

[VectorPostprocessors]
  [rhoep]
    type = LineValueSampler
    variable = rho_edge_pos_1
    start_point = '0 0.005 0'
    end_point = '0.1 0.005 0'
    num_points = 6
    sort_by = x
  []
  [rhoen]
    type = LineValueSampler
    variable = rho_edge_neg_1
    start_point = '0 0.005 0'
    end_point = '0.1 0.005 0'
    num_points = 6
    sort_by = x
  []
[]

[Outputs]
  exodus = true
  interval = 10
  [csv]
    type = CSV
    file_base = rhoe_x_out_l1
    execute_on = final
  []
[]

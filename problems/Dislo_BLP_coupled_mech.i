[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
    xmin = 0.0
    ymin = 0.0
    xmax = 0.01
    ymax = 0.01  
  []
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
  [rhoep]
    initial_condition = 8.e3
  []
  [rhoen]
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
[]

[Modules/TensorMechanics/Master]
  [all]
    # This block adds all of the proper Kernels, strain calculators, and Variables
    # for TensorMechanics in the correct coordinate system (autodetected)
    add_variables = true
    strain = FINITE
    generate_output = 'vonmises_stress stress_xx stress_yy strain_xx strain_yy'
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

[AuxKernels]
  [rhot]
    type = TotalDislocationDensity
    variable = rhot
    execute_on = timestep_end
    rhoep = rhoep
    rhoen = rhoen
  []
  [rhognd]
    type = GNDDislocationDensity
    variable = rhognd
    execute_on = timestep_end
    rhoep = rhoep
    rhoen = rhoen
  []
[]

[Materials]
  [vel]
    type = DisloVelocityCoupled
    nss = 1
    rhoen = rhoen
    rhoep = rhoep
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [xdisp]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0.00000001
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
  end_time = 0.002
  dt = 5.e-7
  dtmin = 1.e-9
[]

[VectorPostprocessors]
  [rhoep]
    type = LineValueSampler
    variable = rhoep
    start_point = '0 0.005 0'
    end_point = '0.01 0.005 0'
    num_points = 11
    sort_by = x
  []
  [rhoen]
    type = LineValueSampler
    variable = rhoen
    start_point = '0 0.005 0'
    end_point = '0.01 0.005 0'
    num_points = 11
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

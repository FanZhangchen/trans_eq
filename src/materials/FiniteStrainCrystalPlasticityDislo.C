#include "FiniteStrainCrystalPlasticityDislo.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TensorMechanicsApp", FiniteStrainCrystalPlasticityDislo);

InputParameters
FiniteStrainCrystalPlasticityDislo::validParams()
{
  InputParameters params = FiniteStrainCrystalPlasticity::validParams();
  params.addClassDescription("Crystal Plasticity with thermal eigenstrain. "
                             "Temperature dependence of the CRSS. "
                             "Dislocation based model. "
                             "Stress dependent dislocation velocity. "
                             "CRSS with Taylor hardening law and bow-out line tension. ");
  // params.addCoupledVar("temp",303.0,"Temperature");
  // params.addCoupledVar("q_t",0.0,"Curvature density (only one slip system)");
  params.addCoupledVar("rho_edge_pos_1", 0.0, "Positive edge dislocation density: slip system 1");
  params.addCoupledVar("rho_edge_neg_1", 0.0, "Negative edge dislocation density: slip system 1");
  params.addCoupledVar("rho_edge_pos_2", 0.0, "Positive edge dislocation density: slip system 2");
  params.addCoupledVar("rho_edge_neg_2", 0.0, "Negative edge dislocation density: slip system 2");
  params.addCoupledVar("rho_forest", 0.0, "Forest dislocation density");
  params.addParam<Real>("thermal_expansion", 0.0, "Thermal expansion coefficient");
  params.addParam<Real>(
      "reference_temperature", 303.0, "reference temperature for thermal expansion");
  params.addParam<Real>("dCRSS_dT_A",
                        1.0,
                        "A coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dCRSS_dT_B",
                        0.0,
                        "B coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dCRSS_dT_C",
                        0.0,
                        "C coefficient for the exponential decrease of the critical "
                        "resolved shear stress with temperature: A + B exp(- C * (T - 293.0))");
  params.addParam<Real>("dislo_mobility", 0.0, "Dislocation mobility");
  params.addParam<Real>("reduced_mobility", 0.0, "Ratio between mobility above vmax and mobility");
  params.addParam<Real>("burgers_vector_mag", 0.257e-6, "Magnitude of the Burgers vector");
  params.addParam<Real>(
      "shear_modulus_hardening", 86000.0, "Shear modulus in Taylor hardening law");
  params.addParam<Real>("dislo_max_velocity", 1000.0, "Maximum dislocation velocity (phonon drag)");
  params.addParam<Real>(
      "bowout_coef", 0.0, "bow-out coefficient: alpha in 4.30 of Hull-Bacon book");
  params.addParam<Real>("bowout_rho_threshold", 0.2, "dislo density threshold to apply bow-out");
  params.addParam<Real>(
      "rho_v_thres", 0.001, "Dislo density threshold below which velocity goes to zero");
  params.addParam<bool>(
      "rho_v_thres_flag", false, "Flag to determine whether to apply the previous threshold");

  params.addParam<Real>("boltzmann", 1.38065e-23, "The Boltzmann Constant");

  params.addParam<Real>("abstemp", 298, "The absolute temperature");

  params.addParam<Real>("p", 0.2, "The flow rule parameter p");

  params.addParam<Real>("q", 1.2, "The flow rule parameter q");

  params.addParam<Real>(
      "tau0hat", 20, "Obtained by extrapolating the lattice friction stress at 0K");

  params.addParam<Real>("gamma0dot", 1.e6, "The flow rule parameter gamma0");

  params.addParam<Real>("F0", 2.77e-19, "Helmholtz free energy of activation");

  params.addParam<Real>("lambda",
                        0.3,
                        "A statistical coefficient which accounts for the deviation from regular "
                        "spatial arrangements of the dislocation");

  params.addParam<Real>("mu", 45.e3, "Shear moduli");

  params.addParam<Real>("burgersvector", 0.257e-6, "The Burgers Vector");

  params.addParam<Real>("taualpha", 2.63, "The resolved shear stress");
  return params;
}

FiniteStrainCrystalPlasticityDislo::FiniteStrainCrystalPlasticityDislo(
    const InputParameters & parameters)
  : FiniteStrainCrystalPlasticity(parameters),

    _rho_edge_pos_1(coupledValue("rho_edge_pos_1")),

    _grad_rhoep1(coupledGradient("rho_edge_pos_1")), // Coupled rhoep gradient

    _rho_edge_neg_1(coupledValue("rho_edge_neg_1")),

    _grad_rhoen1(coupledGradient("rho_edge_neg_1")), // Coupled rhoen gradient

    _rho_edge_pos_2(coupledValue("rho_edge_pos_2")),

    _grad_rhoep2(coupledGradient("rho_edge_pos_2")), // Coupled rhoep gradient

    _rho_edge_neg_2(coupledValue("rho_edge_neg_2")),

    _grad_rhoen2(coupledGradient("rho_edge_neg_2")), // Coupled rhoen gradient

    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _dCRSS_dT_A(getParam<Real>("dCRSS_dT_A")),
    _dCRSS_dT_B(getParam<Real>("dCRSS_dT_B")),
    _dCRSS_dT_C(getParam<Real>("dCRSS_dT_C")),
    _dislo_mobility(getParam<Real>("dislo_mobility")),
    _reduced_mobility(getParam<Real>("reduced_mobility")),
    _burgers_vector_mag(getParam<Real>("burgers_vector_mag")), // Magnitude of the Burgers vector
    _shear_modulus_hardening(
        getParam<Real>("shear_modulus_hardening")), // Shear modulus in Taylor hardening law
    _dislo_max_velocity(
        getParam<Real>("dislo_max_velocity")), // Maximum dislocation velocity (phonon drag)
    _bowout_coef(getParam<Real>("bowout_coef")),
    _bowout_rho_threshold(getParam<Real>("bowout_rho_threshold")),
    _rho_v_thres(
        getParam<Real>("rho_v_thres")), // Dislo density threshold below which velocity goes to zero
    _rho_v_thres_flag(getParam<bool>(
        "rho_v_thres_flag")), // Flag to determine whether to apply the previous threshold

    _gssT(_nss),

    _edge_slip_direction(
        declareProperty<std::vector<Real>>("edge_slip_direction")), // Edge slip directions
    _screw_slip_direction(
        declareProperty<std::vector<Real>>("screw_slip_direction")), // Screw slip direction
    _slip_incr_out(
        declareProperty<std::vector<Real>>("slip_incr_out")), // Slip system increment for output
    _dislo_velocity(declareProperty<std::vector<Real>>("dislo_velocity")), // Dislocation velocity
    _ddislo_velocity_dtau(
        declareProperty<std::vector<Real>>("ddislo_velocity_dtau")), // Derivative of dislo velocity
    _tau_out(declareProperty<std::vector<Real>>("tau_out")), // resolved shear stress for output

    _boltzmann(getParam<Real>("boltzmann")),

    _abstemp(getParam<Real>("abstemp")),

    _p(getParam<Real>("p")),

    _q(getParam<Real>("q")),

    _tau0hat(getParam<Real>("tau0hat")),

    _gamma0dot(getParam<Real>("gamma0dot")),

    _F0(getParam<Real>("F0")),

    _lambda(getParam<Real>("lambda")),

    _mu(getParam<Real>("mu")),

    _taualpha(getParam<Real>("taualpha")),

    _slip_rate_out(declareProperty<std::vector<Real>>("slip_rate")),

    _tau_backstress(_nss),

    _slip_rate(_nss),

    _slip_accum_out(declareProperty<Real>("slip_accum")),

    _slip_accum_out_old(getMaterialPropertyOld<Real>("slip_accum"))
{
}

void
FiniteStrainCrystalPlasticityDislo::calcResidual(RankTwoTensor & resid)
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  // Real temp = _temp[_qp];
  // Real thermal_expansion = _thermal_expansion;
  // Real reference_temperature = _reference_temperature;

  iden.zero();
  iden.addIa(1.0);

  // mooseWarning("_pk2_old", _pk2_old[_qp]);
  // mooseWarning("_fp[_qp]", _fp[_qp]);
  // mooseWarning("_fp_old[_qp]", _fp_old[_qp]);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv
  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // if(_qp == 0)
  // {
  //   mooseWarning("_qp", _qp);
  //   mooseWarning("_dfgrd_tmp", _dfgrd_tmp);
  //   mooseWarning("_fp_prev_inv", _fp_prev_inv);
  //   mooseWarning("_fe", _fe);
  //   mooseWarning("ce", ce);
  //   mooseWarning("_pk2_tmp", _pk2_tmp);
  //   mooseWarning("_fe.det()", _fe.det());
  //   mooseWarning("ce_pk2", ce_pk2);
  // }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);

  // store critical resolved shear stress for output
  _tau_out[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _tau_out[_qp][i] = _tau(i);
  }

  // Introduce the CRSS that will be used in the getSlipIncrements()
  TmpDependCRSS();

  // necessary to call it here because getDisloVelocity
  // depends on slip rate, changing at each iteration
  // of the CP algorithm
  getSlipIncrements(); // Calculate slip rate,dslip,dslipdtau
  // mooseWarning("_slip_incr(0) ", _slip_incr(0));

  // calculate dislocation velocity
  // and store it for advection kernel
  getDisloVelocity();

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
  {
    eqv_slip_incr += _s0[i] * _slip_incr(i);
    // if(_qp == 0)
    //   {
    //     mooseWarning("_s0[i] ", _s0[i]);
    //     // mooseWarning("_slip_incr(i) ", _slip_incr(i));
    //   }
  }

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;

  // if(_qp == 0)
  // {
  //   mooseWarning("_fp_inv", _fp_inv);
  //   mooseWarning("_fp_old_inv", _fp_old_inv);
  //   mooseWarning("eqv_slip_incr", eqv_slip_incr);
  // }

  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;
  // RankTwoTensor thermal_eigenstrain;
  // thermal_eigenstrain =
  //     (1.0 / 2.0) *
  //     (std::exp((2.0 / 3.0) * thermal_expansion * (temp - reference_temperature)) - 1.0) * iden;
  pk2_new = _elasticity_tensor[_qp] * ee;

  // mooseWarning("ee", ee);
  // mooseWarning("pk2_new", pk2_new);

  resid = _pk2_tmp - pk2_new;

  // It would be better to call this function in postSolveQp()
  // so it is not called more times than necessary
  OutputSlipDirection();
}

// Critical resolved shear stress decreases exponentially with temperature
// A + B exp(- C * (T - 303.0))
void
FiniteStrainCrystalPlasticityDislo::TmpDependCRSS()
{
  // Real temp = _temp[_qp];
  // updateGss();

  // Critical resolved shear stress in the input file
  // refers always to room temperature
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gssT[i] = _gss_tmp[i];
  }
}

// Calculate slip increment,dslipdtau
void
FiniteStrainCrystalPlasticityDislo::getSlipIncrements()
{
  std::vector<Real> rho_edge_pos(_nss);
  std::vector<Real> rho_edge_neg(_nss);
  // std::vector<Real> rho_screw_pos(_nss);
  // std::vector<Real> rho_screw_neg(_nss);

  std::vector<Real> rho_edge_pos_grad(_nss);
  std::vector<Real> rho_edge_neg_grad(_nss);

  Real RhoTotSlip; // total dislocation density in the current slip system

  Real _tau_judge;

  // Assign dislocation density vectors
  rho_edge_pos[0] = _rho_edge_pos_1[_qp];
  rho_edge_pos[1] = _rho_edge_pos_2[_qp];
  // rho_edge_pos[2] = _rho_edge_pos_3[_qp];
  // rho_edge_pos[3] = _rho_edge_pos_4[_qp];
  // rho_edge_pos[4] = _rho_edge_pos_5[_qp];
  // rho_edge_pos[5] = _rho_edge_pos_6[_qp];
  // rho_edge_pos[6] = _rho_edge_pos_7[_qp];
  // rho_edge_pos[7] = _rho_edge_pos_8[_qp];
  // rho_edge_pos[8] = _rho_edge_pos_9[_qp];
  // rho_edge_pos[9] = _rho_edge_pos_10[_qp];
  // rho_edge_pos[10] = _rho_edge_pos_11[_qp];
  // rho_edge_pos[11] = _rho_edge_pos_12[_qp];

  rho_edge_neg[0] = _rho_edge_neg_1[_qp];
  rho_edge_neg[1] = _rho_edge_neg_2[_qp];
  // rho_edge_neg[2] = _rho_edge_neg_3[_qp];
  // rho_edge_neg[3] = _rho_edge_neg_4[_qp];
  // rho_edge_neg[4] = _rho_edge_neg_5[_qp];
  // rho_edge_neg[5] = _rho_edge_neg_6[_qp];
  // rho_edge_neg[6] = _rho_edge_neg_7[_qp];
  // rho_edge_neg[7] = _rho_edge_neg_8[_qp];
  // rho_edge_neg[8] = _rho_edge_neg_9[_qp];
  // rho_edge_neg[9] = _rho_edge_neg_10[_qp];
  // rho_edge_neg[10] = _rho_edge_neg_11[_qp];
  // rho_edge_neg[11] = _rho_edge_neg_12[_qp];

  // Assigin dislocation density gradient vectors
  rho_edge_pos_grad[0] = _grad_rhoep1[_qp](1);
  rho_edge_pos_grad[1] = _grad_rhoep2[_qp](1);

  rho_edge_neg_grad[0] = _grad_rhoen1[_qp](1);
  rho_edge_neg_grad[1] = _grad_rhoen2[_qp](1);

  // Positive and negative dislocation give the same
  // contribution to Lp even if their velocity is opposite
  // same for edge and screw
  for (unsigned int i = 0; i < _nss; ++i)
  {
    // temporary variable
    RhoTotSlip = rho_edge_pos[i] + rho_edge_neg[i]; // + rho_screw_pos[i] + rho_screw_neg[i];

    // calculate the backstress term
    _tau_backstress(i) =
        _burgers_vector_mag * _mu * (rho_edge_pos_grad[i] - rho_edge_neg_grad[i]) / RhoTotSlip;
    // mooseWarning("_tau_backstress(i) ", _tau_backstress(i));

    if (RhoTotSlip > 0.0)
    {

      _tau_judge = std::abs(_tau(i) - _tau_backstress(i));

      if (_tau_judge > _gssT[i])
      {
        _slip_rate(i) =
            _gamma0dot *
            std::exp(-_F0 / _boltzmann / _abstemp *
                     std::pow((1.0 - std::pow(((std::abs(_tau(i) - _tau_backstress(i)) - _gssT[i]) /
                                               _tau0hat),
                                              _p)),
                              _q)) *
            std::copysign(1.0, (_tau(i) - _tau_backstress(i)));
        // mooseWarning("_gssT[i] ", _gssT[i]);
        // mooseWarning("_slip_rate(i)11111 ", _slip_rate(i));
      }
      else
      {
        _slip_rate(i) = 0.0;
        // mooseWarning("_slip_rate(i)22222 ", _slip_rate(i));
      }

      // mooseWarning("_dt11111 ", _dt);
      _slip_incr(i) = _slip_rate(i) * _dt;
      // mooseWarning("_slip_incr(i)11111 ", _slip_incr(i));

      // mooseWarning("slip resistance ", _gssT[i]);
      // mooseWarning("tau ", _tau(i));
      // mooseWarning("The slip rate value is ", _slip_rate(i));

      // Derivative is always positive
      if (_tau_judge > _gssT[i])
      {
        _dslipdtau(i) =
            _gamma0dot * _p * _q * _F0 / _boltzmann / _abstemp *
            std::exp(-_F0 / _boltzmann / _abstemp *
                     std::pow((1.0 - std::pow(((std::abs(_tau(i) - _tau_backstress(i)) - _gssT[i]) /
                                               _tau0hat),
                                              _p)),
                              _q)) *
            std::pow(
                (1.0 -
                 std::pow(((std::abs(_tau(i) - _tau_backstress(i)) - _gssT[i]) / _tau0hat), _p)),
                _q - 1.0) *
            std::pow(((std::abs(_tau(i) - _tau_backstress(i)) - _gssT[i]) / _tau0hat), _p - 1.0) *
            std::copysign(1.0, (_tau(i) - _tau_backstress(i))) * _dt;
      }
      else
      {
        _dslipdtau(i) = 0.0;
      }
    }
    else
    {

      _slip_incr(i) = 0.0;

      _dslipdtau(i) = 0.0;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (std::abs(_slip_incr(i)) > _slip_incr_tol)
    {
      //_err_tol = true;
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
      _slip_incr(i) = _slip_incr_tol * std::copysign(1.0, _tau(i));
      _dslipdtau(i) = 0.0;
      // mooseWarning("_slip_incr(0)22222 ", _slip_incr(i));
    }
  }

  // store slip increment and slip rate for output
  _slip_incr_out[_qp].resize(_nss);

  _slip_rate_out[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_rate_out[_qp][i] = _slip_rate(i);

    _slip_incr_out[_qp][i] = _slip_incr(i);

    _acc_slip[_qp] = std::abs(_slip_incr(i)) + _acc_slip_old[_qp];
  }
}

// Calculate dislocation velocity (edge and screw) as a function
// of the resolved shear stress and its derivative
void
FiniteStrainCrystalPlasticityDislo::getDisloVelocity()
{
  Real tau0; // resolved shear stress at max velocity _dislo_max_velocity

  std::vector<Real> rho_edge_pos(_nss);
  std::vector<Real> rho_edge_neg(_nss);
  // std::vector<Real> rho_screw_pos(_nss);
  // std::vector<Real> rho_screw_neg(_nss);

  Real RhoTotSlip = 0.0; // total dislocation density in the current slip system

  // Real rho_v_thres = _rho_v_thres; // below this threshold rho_tot the velocity decreases to zero

  _dislo_velocity[_qp].resize(_nss);
  _ddislo_velocity_dtau[_qp].resize(_nss);

  // initialize the dislocation velocity
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] = 0.0;
    _ddislo_velocity_dtau[_qp][i] = 0.0;
  }

  // Assign dislocation density vectors
  rho_edge_pos[0] = _rho_edge_pos_1[_qp];
  rho_edge_pos[1] = _rho_edge_pos_2[_qp];
  // rho_edge_pos[2] = _rho_edge_pos_3[_qp];
  // rho_edge_pos[3] = _rho_edge_pos_4[_qp];
  // rho_edge_pos[4] = _rho_edge_pos_5[_qp];
  // rho_edge_pos[5] = _rho_edge_pos_6[_qp];
  // rho_edge_pos[6] = _rho_edge_pos_7[_qp];
  // rho_edge_pos[7] = _rho_edge_pos_8[_qp];
  // rho_edge_pos[8] = _rho_edge_pos_9[_qp];
  // rho_edge_pos[9] = _rho_edge_pos_10[_qp];
  // rho_edge_pos[10] = _rho_edge_pos_11[_qp];
  // rho_edge_pos[11] = _rho_edge_pos_12[_qp];

  rho_edge_neg[0] = _rho_edge_neg_1[_qp];
  rho_edge_neg[1] = _rho_edge_neg_2[_qp];
  // rho_edge_neg[2] = _rho_edge_neg_3[_qp];
  // rho_edge_neg[3] = _rho_edge_neg_4[_qp];
  // rho_edge_neg[4] = _rho_edge_neg_5[_qp];
  // rho_edge_neg[5] = _rho_edge_neg_6[_qp];
  // rho_edge_neg[6] = _rho_edge_neg_7[_qp];
  // rho_edge_neg[7] = _rho_edge_neg_8[_qp];
  // rho_edge_neg[8] = _rho_edge_neg_9[_qp];
  // rho_edge_neg[9] = _rho_edge_neg_10[_qp];
  // rho_edge_neg[10] = _rho_edge_neg_11[_qp];
  // rho_edge_neg[11] = _rho_edge_neg_12[_qp];

  for (unsigned int i = 0; i < _nss; ++i)
  {

    tau0 = 0.0;

    RhoTotSlip = rho_edge_pos[i] + rho_edge_neg[i]; // + rho_screw_pos[i] + rho_screw_neg[i];

    // if (_dislo_mobility > 0.0) {
    //   tau0 = _dislo_max_velocity / _dislo_mobility; // temporary variable for this slip system
    tau0 += _gssT[i];
    // }

    if (std::abs(_tau(i)) > tau0)
    { // Case above _dislo_max_velocity: use reduced mobility

      _dislo_velocity[_qp][i] = _slip_rate(i) / _burgers_vector_mag / RhoTotSlip;

      // Derivative is always positive
      // _ddislo_velocity_dtau[_qp][i] = _reduced_mobility;

      // } else if (std::abs(_tau(i)) > _gssT[i]) { // Case below _dislo_max_velocity

      //     _dislo_velocity[_qp][i] = _dislo_mobility * (std::abs(_tau(i)) - _gssT[i])
      //                           * std::copysign(1.0, _tau(i));

      //   // Derivative is always positive
      //   _ddislo_velocity_dtau[_qp][i] = _dislo_mobility;

      //   if (_rho_v_thres_flag) { // Case with density below threshold

      //       if (RhoTotSlip < rho_v_thres) { // rescale dislocation velocity and derivative by a
      //       factor

      //         _dislo_velocity[_qp][i] *= (RhoTotSlip / rho_v_thres);
      // 	  _ddislo_velocity_dtau[_qp][i] *= (RhoTotSlip / rho_v_thres);

      //       }

      //     }
    }
    else
    { // Case below critical resolved shear stress

      _dislo_velocity[_qp][i] = 0.0;
      // _ddislo_velocity_dtau[_qp][i] = 0.0;
    }

  } // end cycle over slip systems
}

// Store slip direction
// to couple with dislocation transport
void
FiniteStrainCrystalPlasticityDislo::OutputSlipDirection()
{
  DenseVector<Real> mo(LIBMESH_DIM * _nss);
  DenseVector<Real> no(LIBMESH_DIM * _nss);

  // Temporary directions and normals to calculate
  // screw dislocation slip direction
  RealVectorValue temp_mo;
  RealVectorValue temp_no;
  RealVectorValue temp_screw_mo;

  // Update slip direction with crystal orientation
  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(i * LIBMESH_DIM + j) =
            mo(i * LIBMESH_DIM + j) + _crysrot[_qp](j, k) * _mo(i * LIBMESH_DIM + k);
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(i * LIBMESH_DIM + j) =
            no(i * LIBMESH_DIM + j) + _crysrot[_qp](j, k) * _no(i * LIBMESH_DIM + k);
    }
  }

  _edge_slip_direction[_qp].resize(LIBMESH_DIM * _nss);
  _screw_slip_direction[_qp].resize(LIBMESH_DIM * _nss);

  // Store slip direction (already normalized)
  // for edge and screw dislocations
  // to couple with dislocation transport
  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      temp_mo(j) = mo(i * LIBMESH_DIM + j);
      temp_no(j) = no(i * LIBMESH_DIM + j);
    }

    temp_screw_mo = temp_mo.cross(temp_no);

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      _edge_slip_direction[_qp][i * LIBMESH_DIM + j] = mo(i * LIBMESH_DIM + j);
      _screw_slip_direction[_qp][i * LIBMESH_DIM + j] = temp_screw_mo(j);
    }
  }
}

/**
 * Calculate slip system resistance (CRSS)
 * based on Taylor hardening model
 */
void
FiniteStrainCrystalPlasticityDislo::updateGss()
{
  std::vector<Real> sres(_nss);           // Taylor hardening + bow-out line tension
  std::vector<Real> TotalRho(_nss); // total dislocation density
  // Real rho_forest; // forest dislocation density

  // Real q_t = _q_t[_qp]; // curvature density of the active slip system (only 1)

  std::vector<Real> rho_edge_pos(_nss);
  std::vector<Real> rho_edge_neg(_nss);
  // std::vector<Real> rho_screw_pos(_nss);
  // std::vector<Real> rho_screw_neg(_nss);

  // Assign dislocation density vectors
  rho_edge_pos[0] = _rho_edge_pos_1[_qp];
  rho_edge_pos[1] = _rho_edge_pos_2[_qp];
  // rho_edge_pos[2] = _rho_edge_pos_3[_qp];
  // rho_edge_pos[3] = _rho_edge_pos_4[_qp];
  // rho_edge_pos[4] = _rho_edge_pos_5[_qp];
  // rho_edge_pos[5] = _rho_edge_pos_6[_qp];
  // rho_edge_pos[6] = _rho_edge_pos_7[_qp];
  // rho_edge_pos[7] = _rho_edge_pos_8[_qp];
  // rho_edge_pos[8] = _rho_edge_pos_9[_qp];
  // rho_edge_pos[9] = _rho_edge_pos_10[_qp];
  // rho_edge_pos[10] = _rho_edge_pos_11[_qp];
  // rho_edge_pos[11] = _rho_edge_pos_12[_qp];

  rho_edge_neg[0] = _rho_edge_neg_1[_qp];
  rho_edge_neg[1] = _rho_edge_neg_2[_qp];
  // rho_edge_neg[2] = _rho_edge_neg_3[_qp];
  // rho_edge_neg[3] = _rho_edge_neg_4[_qp];
  // rho_edge_neg[4] = _rho_edge_neg_5[_qp];
  // rho_edge_neg[5] = _rho_edge_neg_6[_qp];
  // rho_edge_neg[6] = _rho_edge_neg_7[_qp];
  // rho_edge_neg[7] = _rho_edge_neg_8[_qp];
  // rho_edge_neg[8] = _rho_edge_neg_9[_qp];
  // rho_edge_neg[9] = _rho_edge_neg_10[_qp];
  // rho_edge_neg[10] = _rho_edge_neg_11[_qp];
  // rho_edge_neg[11] = _rho_edge_neg_12[_qp];

  // rho_forest = _rho_forest[_qp];

  // Is this update necessary
  // for the constitutive model
  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i = 0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr(i));

  for (unsigned int i = 0; i < _nss; ++i)
    {
      for (unsigned int j = 0; j < _nss; ++j)
      {
        if (j == i)
        TotalRho[j] = (rho_edge_pos[i] + rho_edge_neg[i]);
      }
    }

  // TotalRho += rho_forest;
  for (unsigned int i = 0; i < _nss; ++i)
    {
      if (TotalRho[i] >= 0.0)
        {
          sres[i] = _lambda * _mu * _burgers_vector_mag * std::sqrt(TotalRho[i]);
        }
        else
        {
          sres[i] = 0.0;
        }
    }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _gss_tmp[i] = sres[i];
  }
}

void
FiniteStrainCrystalPlasticityDislo::postSolveQp()
{
  if (_err_tol)
  {
    _err_tol = false;
    if (_gen_rndm_stress_flag)
    {
      if (!_input_rndm_scale_var)
        _rndm_scale_var = _elasticity_tensor[_qp](0, 0, 0, 0);

      _stress[_qp] = RankTwoTensor::genRandomSymmTensor(_rndm_scale_var, 1.0);
    }
    else
      mooseError("FiniteStrainCrystalPlasticity: Constitutive failure");
  }
  else
  {
    _stress[_qp] = _fe * _pk2[_qp] * _fe.transpose() / _fe.det();

    _Jacobian_mult[_qp] += calcTangentModuli(); // Calculate jacobian for preconditioner

    RankTwoTensor iden(RankTwoTensor::initIdentity);

    _lag_e[_qp] = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] - iden;
    _lag_e[_qp] = _lag_e[_qp] * 0.5;

    RankTwoTensor rot;
    rot = get_current_rotation(_deformation_gradient[_qp]); // Calculate material rotation
    _update_rot[_qp] = rot * _crysrot[_qp];
  }
}

// Calls getMatRot to perform RU factorization of a tensor.
RankTwoTensor
FiniteStrainCrystalPlasticityDislo::get_current_rotation(const RankTwoTensor & a)
{
  return getMatRot(a);
}

// Performs RU factorization of a tensor
// Added debug information when RU decomposition fails
RankTwoTensor
FiniteStrainCrystalPlasticityDislo::getMatRot(const RankTwoTensor & a)
{
  RankTwoTensor rot;
  RankTwoTensor c, diag, evec;
  PetscScalar cmat[LIBMESH_DIM][LIBMESH_DIM], work[10];
  PetscReal w[LIBMESH_DIM];
  PetscBLASInt nd = LIBMESH_DIM, lwork = 10, info;

  c = a.transpose() * a;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      cmat[i][j] = c(i, j);

  LAPACKsyev_("V", "U", &nd, &cmat[0][0], &nd, w, work, &lwork, &info);

  if (info != 0)
  {
    mooseWarning("Deformation gradient components (0,0)", _deformation_gradient[_qp](0, 0));
    mooseWarning("Deformation gradient components (0,1)", _deformation_gradient[_qp](0, 1));
    mooseWarning("Deformation gradient components (0,2)", _deformation_gradient[_qp](0, 2));
    mooseWarning("Deformation gradient components (1,0)", _deformation_gradient[_qp](1, 0));
    mooseWarning("Deformation gradient components (1,1)", _deformation_gradient[_qp](1, 1));
    mooseWarning("Deformation gradient components (1,2)", _deformation_gradient[_qp](1, 2));
    mooseWarning("Deformation gradient components (2,0)", _deformation_gradient[_qp](2, 0));
    mooseWarning("Deformation gradient components (2,1)", _deformation_gradient[_qp](2, 1));
    mooseWarning("Deformation gradient components (2,2)", _deformation_gradient[_qp](2, 2));
    mooseWarning("Rho edge pos 1 ", _rho_edge_pos_1[_qp]);
    mooseWarning("Rho edge pos 2 ", _rho_edge_pos_2[_qp]);
    // mooseWarning("Rho edge pos 3 ", _rho_edge_pos_3[_qp]);
    // mooseWarning("Rho edge pos 4 ", _rho_edge_pos_4[_qp]);
    // mooseWarning("Rho edge pos 5 ", _rho_edge_pos_5[_qp]);
    // mooseWarning("Rho edge pos 6 ", _rho_edge_pos_6[_qp]);
    // mooseWarning("Rho edge pos 7 ", _rho_edge_pos_7[_qp]);
    // mooseWarning("Rho edge pos 8 ", _rho_edge_pos_8[_qp]);
    // mooseWarning("Rho edge pos 9 ", _rho_edge_pos_9[_qp]);
    // mooseWarning("Rho edge pos 10 ", _rho_edge_pos_10[_qp]);
    // mooseWarning("Rho edge pos 11 ", _rho_edge_pos_11[_qp]);
    // mooseWarning("Rho edge pos 12 ", _rho_edge_pos_12[_qp]);
    mooseWarning("Rho edge neg 1 ", _rho_edge_neg_1[_qp]);
    mooseWarning("Rho edge neg 2 ", _rho_edge_neg_2[_qp]);
    // mooseWarning("Rho edge neg 3 ", _rho_edge_neg_3[_qp]);
    // mooseWarning("Rho edge neg 4 ", _rho_edge_neg_4[_qp]);
    // mooseWarning("Rho edge neg 5 ", _rho_edge_neg_5[_qp]);
    // mooseWarning("Rho edge neg 6 ", _rho_edge_neg_6[_qp]);
    // mooseWarning("Rho edge neg 7 ", _rho_edge_neg_7[_qp]);
    // mooseWarning("Rho edge neg 8 ", _rho_edge_neg_8[_qp]);
    // mooseWarning("Rho edge neg 9 ", _rho_edge_neg_9[_qp]);
    // mooseWarning("Rho edge neg 10 ", _rho_edge_neg_10[_qp]);
    // mooseWarning("Rho edge neg 11 ", _rho_edge_neg_11[_qp]);
    // mooseWarning("Rho edge neg 12 ", _rho_edge_neg_12[_qp]);
    // mooseWarning("Rho screw pos 1 ", _rho_screw_pos_1[_qp]);
    // mooseWarning("Rho screw pos 2 ", _rho_screw_pos_2[_qp]);
    // mooseWarning("Rho screw pos 3 ", _rho_screw_pos_3[_qp]);
    // mooseWarning("Rho screw pos 4 ", _rho_screw_pos_4[_qp]);
    // mooseWarning("Rho screw pos 5 ", _rho_screw_pos_5[_qp]);
    // mooseWarning("Rho screw pos 6 ", _rho_screw_pos_6[_qp]);
    // mooseWarning("Rho screw pos 7 ", _rho_screw_pos_7[_qp]);
    // mooseWarning("Rho screw pos 8 ", _rho_screw_pos_8[_qp]);
    // mooseWarning("Rho screw pos 9 ", _rho_screw_pos_9[_qp]);
    // mooseWarning("Rho screw pos 10 ", _rho_screw_pos_10[_qp]);
    // mooseWarning("Rho screw pos 11 ", _rho_screw_pos_11[_qp]);
    // mooseWarning("Rho screw pos 12 ", _rho_screw_pos_12[_qp]);
    // mooseWarning("Rho screw neg 1 ", _rho_screw_neg_1[_qp]);
    // mooseWarning("Rho screw neg 2 ", _rho_screw_neg_2[_qp]);
    // mooseWarning("Rho screw neg 3 ", _rho_screw_neg_3[_qp]);
    // mooseWarning("Rho screw neg 4 ", _rho_screw_neg_4[_qp]);
    // mooseWarning("Rho screw neg 5 ", _rho_screw_neg_5[_qp]);
    // mooseWarning("Rho screw neg 6 ", _rho_screw_neg_6[_qp]);
    // mooseWarning("Rho screw neg 7 ", _rho_screw_neg_7[_qp]);
    // mooseWarning("Rho screw neg 8 ", _rho_screw_neg_8[_qp]);
    // mooseWarning("Rho screw neg 9 ", _rho_screw_neg_9[_qp]);
    // mooseWarning("Rho screw neg 10 ", _rho_screw_neg_10[_qp]);
    // mooseWarning("Rho screw neg 11 ", _rho_screw_neg_11[_qp]);
    // mooseWarning("Rho screw neg 12 ", _rho_screw_neg_12[_qp]);
    mooseWarning("Slip increment ", _slip_incr(0));
    // mooseWarning("Slip increment ", _slip_incr(1));
    // mooseWarning("Slip increment ", _slip_incr(2));
    // mooseWarning("Slip increment ", _slip_incr(3));
    // mooseWarning("Slip increment ", _slip_incr(4));
    // mooseWarning("Slip increment ", _slip_incr(5));
    // mooseWarning("Slip increment ", _slip_incr(6));
    // mooseWarning("Slip increment ", _slip_incr(7));
    // mooseWarning("Slip increment ", _slip_incr(8));
    // mooseWarning("Slip increment ", _slip_incr(9));
    // mooseWarning("Slip increment ", _slip_incr(10));
    // mooseWarning("Slip increment ", _slip_incr(11));
    mooseError(
        "FiniteStrainCrystalPlasticityDislo: DSYEV function call in getMatRot function failed");
  }

  diag.zero();

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    diag(i, i) = std::sqrt(w[i]);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      evec(i, j) = cmat[i][j];

  rot = a * ((evec.transpose() * diag * evec).inverse());

  return rot;
}

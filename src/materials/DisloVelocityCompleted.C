// Zhangchen Fan
// HITSZ - CMMC
// May 2023

#include "DisloVelocityCompleted.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TransEQApp", DisloVelocityCompleted);

InputParameters
DisloVelocityCompleted::validParams()
{
  InputParameters params = Material::validParams();

  params.addClassDescription("Dislocation velocity for Dislocation transport equation. ");

  params.addRequiredParam<int>("nss", "Number of slip systems");

  params.addRequiredCoupledVar("rhoe1", "positive edge dislocation density");

  params.addRequiredCoupledVar("rhoe2", "negative edge dislocation density");

  params.addRequiredCoupledVar("rhoe3", "positive edge dislocation density");

  params.addRequiredCoupledVar("rhoe4", "negative edge dislocation density");

  params.addRequiredCoupledVar("rhos1", "positive screw dislocation density");

  params.addRequiredCoupledVar("rhos2", "negative screw dislocation density");

  params.addRequiredCoupledVar("rhos3", "positive screw dislocation density");

  params.addRequiredCoupledVar("rhos4", "negative screw dislocation density");

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

  params.addParam<Real>("taualpha", 2.36, "The resolved shear stress");

  // params.addParam<std::vector<Real>>("rho_edge", 16000, "The total edge dislocation density");

  return params;
}

DisloVelocityCompleted::DisloVelocityCompleted(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),

    _nss(getParam<int>("nss")),

    _gssT(_nss),

    // _slip_rate(getParam<Real>("slip_rate")),

    _dislo_velocity(declareProperty<std::vector<Real>>(
        "dislo_velocity")), // Dislocation velocity at current time step t

    _velocity_old(
        getMaterialPropertyOld<std::vector<Real>>("dislo_velocity")), // Dislocation velocity at t-1

    _rhoe1(coupledValue("rhoe1")), // Coupled rhoep

    _grad_rhoe1(coupledGradient("rhoe1")), // Coupled rhoep gradient

    _rhoe2(coupledValue("rhoe2")), // Coupled rhoen

    _grad_rhoe2(coupledGradient("rhoe3")), // Coupled rhoen gradient

    _rhoe3(coupledValue("rhoe3")), // Coupled rhoep

    _grad_rhoe3(coupledGradient("rhoe3")), // Coupled rhoep gradient

    _rhoe4(coupledValue("rhoe4")), // Coupled rhoen

    _grad_rhoe4(coupledGradient("rhoe4")), // Coupled rhoen gradient

    _rhos1(coupledValue("rhos1")), // Coupled rhosp

    _grad_rhos1(coupledGradient("rhos1")), // Coupled rhosp gradient

    _rhos2(coupledValue("rhos2")), // Coupled rhosn

    _grad_rhos2(coupledGradient("rhos2")), // Coupled rhosn gradient

    _rhos3(coupledValue("rhos3")), // Coupled rhosp

    _grad_rhos3(coupledGradient("rhos3")), // Coupled rhosp gradient

    _rhos4(coupledValue("rhos4")), // Coupled rhosn

    _grad_rhos4(coupledGradient("rhos4")), // Coupled rhosn gradient

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

    _burgersvector(getParam<Real>("burgersvector")),

    _rho_edge(declareProperty<Real>("rho_edge")),

    _rho_screw(declareProperty<Real>("rho_screw")),

    _rhot(declareProperty<Real>("rhot")),

    _tau_backstress(declareProperty<Real>("tau_backstress")),

    _slip_rate(declareProperty<Real>("slip_rate"))

{
}

void
DisloVelocityCompleted::computeQpProperties()
{

  _dislo_velocity[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] = 0.0;
  }

  // initialize the edge dislocation density

  _rho_edge[_qp] = _rhoe1[_qp] + _rhoe1[_qp] + _rhoe3[_qp] + _rhoe4[_qp];

  _rho_screw[_qp] = _rhos1[_qp] + _rhos2[_qp] + _rhos3[_qp] + _rhos4[_qp];

  _rhot[_qp] = _rho_edge[_qp] + _rho_screw[_qp];

  _tau_backstress[_qp] =
      _burgersvector * _mu * (_grad_rhoe1[_qp](0) - _grad_rhoe2[_qp](0) - _grad_rhoe3[_qp](0) + _grad_rhoe4[_qp](0) + _grad_rhos1[_qp](1) + _grad_rhos2[_qp](1) - _grad_rhos3[_qp](1) - _grad_rhos4[_qp](1)) / _rhot[_qp];

  _slip_rate[_qp] =
      _gamma0dot *
      std::exp(-_F0 / _boltzmann / _abstemp *
               std::pow((1 - std::pow(((std::abs(_taualpha - _tau_backstress[_qp]) -
                                        _lambda * _mu * _burgersvector * std::sqrt(_rhot[_qp])) /
                                       _tau0hat),
                                      _p)),
                        _q)) *
      std::copysign(1.0, (_taualpha - _tau_backstress[_qp]));

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] =
        _slip_rate[_qp] / _burgersvector / (_rho_edge[_qp] + 0.5 * _rho_screw[_qp]);
  }
}

void
DisloVelocityCompleted::initQpStatefulProperties()
{

  _dislo_velocity[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] = 0.0;
  }

  // initialize the edge dislocation density

  _rho_edge[_qp] = _rhoe1[_qp] + _rhoe1[_qp] + _rhoe3[_qp] + _rhoe4[_qp];

  _rho_screw[_qp] = _rhos1[_qp] + _rhos2[_qp] + _rhos3[_qp] + _rhos4[_qp];

  _rhot[_qp] = _rho_edge[_qp] + _rho_screw[_qp];

  _tau_backstress[_qp] =
      _burgersvector * _mu * (_grad_rhoe1[_qp](0) - _grad_rhoe2[_qp](0) - _grad_rhoe3[_qp](0) + _grad_rhoe4[_qp](0) + _grad_rhos1[_qp](1) + _grad_rhos2[_qp](1) - _grad_rhos3[_qp](1) - _grad_rhos4[_qp](1)) / _rhot[_qp];

  _slip_rate[_qp] =
      _gamma0dot *
      std::exp(-_F0 / _boltzmann / _abstemp *
               std::pow((1 - std::pow(((std::abs(_taualpha - _tau_backstress[_qp]) -
                                        _lambda * _mu * _burgersvector * std::sqrt(_rhot[_qp])) /
                                       _tau0hat),
                                      _p)),
                        _q)) *
      std::copysign(1.0, (_taualpha - _tau_backstress[_qp]));

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] =
        _slip_rate[_qp] / _burgersvector / (_rho_edge[_qp] + 0.5 * _rho_screw[_qp]);
  }
}
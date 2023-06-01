// Zhangchen Fan
// HITSZ - CMMC
// May 2023

#include "DisloVelocityCoupled_ConstSlipRate.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TransEQApp", DisloVelocityCoupled_ConstSlipRate);

InputParameters
DisloVelocityCoupled_ConstSlipRate::validParams()
{
  InputParameters params = Material::validParams();

  params.addClassDescription("Dislocation velocity for Dislocation transport equation. ");

  params.addRequiredParam<int>("nss", "Number of slip systems");

  params.addParam<Real>("slip_rate", 10e6, "The velocity of the flux");

  params.addRequiredCoupledVar("rhoep", "positive edge dislocation density");

  params.addRequiredCoupledVar("rhoen", "negative edge dislocation density");

  return params;
}

DisloVelocityCoupled_ConstSlipRate::DisloVelocityCoupled_ConstSlipRate(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),

    _nss(getParam<int>("nss")),

      _gssT(_nss),

    _slip_rate(getParam<Real>("slip_rate")),

    _dislo_velocity(declareProperty<std::vector<Real>>("dislo_velocity")), // Dislocation velocity at current time step t

    _velocity_old(getMaterialPropertyOld<std::vector<Real>>("dislo_velocity")), // Dislocation velocity at t-1

    _rhoep(coupledValue("rhoep")), // Coupled rhoep

    _rhoen(coupledValue("rhoen")) // Coupled rhoen
{
}

void
DisloVelocityCoupled_ConstSlipRate::initQpStatefulProperties()
{
  
  _dislo_velocity[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] = _slip_rate/(_rhoep[_qp] + _rhoen[_qp]);
  }
  
}

void
DisloVelocityCoupled_ConstSlipRate::computeQpProperties()
{
  // Real tau0; // resolved shear stress at max velocity _dislo_max_velocity 
  
  // std::vector<Real> rho_edge_pos(_nss);
  // std::vector<Real> rho_edge_neg(_nss);
  // std::vector<Real> rho_screw_pos(_nss);
  // std::vector<Real> rho_screw_neg(_nss);

  // Real RhoTotSlip = 0.0; // total dislocation density in the current slip system

  // Real rho_v_thres = _rho_v_thres; // below this threshold rho_tot the velocity decreases to zero  
  
  _dislo_velocity[_qp].resize(_nss);
  // _ddislo_velocity_dtau[_qp].resize(_nss);
  
  // for (unsigned int i = 0; i < _nss; ++i) // initial the _dislo_velocity
  // {
  //   _dislo_velocity[_qp][i] = 0.0;
  //   // _ddislo_velocity_dtau[_qp][i] = 0.0;
  // }
  
  for (unsigned int i = 0; i < _nss; ++i)
  {

    _dislo_velocity[_qp][i] = _velocity_old[_qp][i];

  } // end cycle over slip systems
}
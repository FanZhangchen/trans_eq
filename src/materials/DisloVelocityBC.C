// Zhangchen Fan
// HITSZ - CMMC
// May 2023

#include "DisloVelocityBC.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

#include <fstream>
#include <cmath>

registerMooseObject("TransEQApp", DisloVelocityBC);

InputParameters
DisloVelocityBC::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<int>("nss", "Number of slip systems");

  params.addParam<Real>("initial_velocity", 0.0, "The velocity of the flux");

  params.addRequiredCoupledVar("rhoep", "positive edge dislocation density");

  params.addRequiredCoupledVar("rhoen", "negative edge dislocation density");

  return params;
}

DisloVelocityBC::DisloVelocityBC(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // _initial_velocity(getParam<Real>("initial_velocity")),

    _nss(getParam<int>("nss")),

    _gssT(_nss),

    _initial_velocity(getParam<Real>("initial_velocity")),

    _dislo_velocity(declareProperty<std::vector<Real>>("dislo_velocity")), // Dislocation velocity

    _velocity_old(getMaterialPropertyOld<std::vector<Real>>("dislo_velocity")),

    _rhoep(coupledValue("rhoep")), // Coupled rhoep

    _rhoen(coupledValue("rhoen")) // Coupled rhoen
{
}

void
DisloVelocityBC::initQpStatefulProperties()
{

  _dislo_velocity[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dislo_velocity[_qp][i] = _initial_velocity;
  }
}

void
DisloVelocityBC::computeQpProperties()
{

  _dislo_velocity[_qp].resize(_nss);
  // _ddislo_velocity_dtau[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i) // initial the _dislo_velocity
  {
    _dislo_velocity[_qp][i] = _velocity_old[_qp][i];
    // _ddislo_velocity_dtau[_qp][i] = 0.0;
  }

  // for (unsigned int i = 0; i < _nss; ++i)
  // {

  //   _dislo_velocity[_qp][i] = 0.0;
  // }
}
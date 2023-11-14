//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CDDSource.h"

registerMooseObject("TransEQApp", CDDSource);

InputParameters
CDDSource::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addCoupledVar("rho_screw_1",
                       "The variable representing the quandrant 1 edge dislocation density.");
  return params;
}

CDDSource::CDDSource(const InputParameters & parameters)
  : Kernel(parameters),
  _velocity(getMaterialProperty<Real>("velocity"))
{
}

Real
CDDSource::computeQpResidual()
{
  return -_grad_test[_i][_qp](0)*_u[_qp]*_velocity[_qp];
}

Real
CDDSource::computeQpJacobian()
{
  return _grad_test[_i][_qp](0)*_grad_phi[_j][_qp](0)*_velocity[_qp];
}

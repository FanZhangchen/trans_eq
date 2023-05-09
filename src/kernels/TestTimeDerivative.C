//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TestTimeDerivative.h"

#include "Material.h"

registerMooseObject("TransEQApp", TestTimeDerivative);

InputParameters
TestTimeDerivative::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addParam<Real>("time_coefficient", 1.0, "Time Coefficient");
  return params;
}

TestTimeDerivative::TestTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters), _time_coefficient(getParam<Real>("time_coefficient"))
{
}

Real
TestTimeDerivative::computeQpResidual()
{
  return _time_coefficient * TimeDerivative::computeQpResidual();
}

Real
TestTimeDerivative::computeQpJacobian()
{
  return _time_coefficient * TimeDerivative::computeQpJacobian();
}

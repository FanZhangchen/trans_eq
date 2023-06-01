//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GNDDislocationDensity.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("TransEQApp", GNDDislocationDensity);

InputParameters
GNDDislocationDensity::validParams()
{
  InputParameters params = VectorAuxKernel::validParams();

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("rhoep", "The positive edge dislocation.");

  params.addRequiredCoupledVar("rhoen", "The negative edge dislocation.");

  return params;
}

GNDDislocationDensity::GNDDislocationDensity(const InputParameters & parameters)
  : VectorAuxKernel(parameters),

    // Get the value of the variables
    _rhoep(coupledValue("rhoep")),

    _rhoen(coupledValue("rhoen"))

{
}

RealVectorValue
GNDDislocationDensity::computeValue()
{
  // Access the gradient of the pressure at this quadrature point, then pull out the "component" of
  // it requested (x, y or z). Note, that getting a particular component of a gradient is done using
  // the parenthesis operator.
  return (_rhoep[_qp] - _rhoen[_qp])/16000;
}

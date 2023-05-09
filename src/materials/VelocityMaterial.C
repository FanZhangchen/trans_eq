//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VelocityMaterial.h"

registerMooseObject("TransEQApp", VelocityMaterial);

InputParameters
VelocityMaterial::validParams()
{
  InputParameters params = Material::validParams();

  // Allow users to specify vectors defining the points of a piecewise function formed via linear
  // interpolation.
  // params.addRequiredParam<std::vector<Real>>(
  //     "independent_vals",
  //     "The vector of z-coordinate values for a piecewise function's independent variable");
  // params.addRequiredParam<std::vector<Real>>(
  //     "dependent_vals", "The vector of diffusivity values for a piecewise function's dependent");
  // Allow the user to specify which independent variable's gradient to use for calculating the
  // convection velocity property:
  // params.addCoupledVar(
  //     "diffusion_gradient",
  //     "The gradient of this variable will be used to compute a velocity vector property.");
  params.addParam<Real>("initial_velocity", 7.5e-1, "The velocity of the flux");

  return params;
}

VelocityMaterial::VelocityMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _initial_velocity(getParam<Real>("initial_velocity")),

    _velocity(declareProperty<Real>("velocity")),

    _velocity_old(getMaterialPropertyOld<Real>("velocity"))
{
}

void
VelocityMaterial::initQpStatefulProperties()
{
  _velocity[_qp] = _initial_velocity;
}

void
VelocityMaterial::computeQpProperties()
{
  _velocity[_qp] = _velocity_old[_qp];
}

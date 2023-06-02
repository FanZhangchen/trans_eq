//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "LinearInterpolation.h"
#include "DerivativeMaterialInterface.h"

class DisloVelocityBC_ConstSlipRate : public DerivativeMaterialInterface<Material>
{
public:
  DisloVelocityBC_ConstSlipRate(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();

  const unsigned int _nss;

  std::vector<Real> _gssT;

private:
  /// member variable to hold the computed diffusivity coefficient
  MaterialProperty<std::vector<Real>> & _dislo_velocity;
  /// member variable to hold the computed convection velocity gradient term
  const MaterialProperty<std::vector<Real>> & _velocity_old;

  const VariableValue & _rhoep;

  const VariableValue & _rhoen;

  Real _initial_velocity;
};

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

class CDDSink : public Kernel
{
public:
  CDDSink(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  const MaterialProperty<std::vector<Real>> & _slip_rate;

  const VariableValue & _rho_edge_1;

  const VariableValue & _rho_edge_2;

  const VariableValue & _rho_edge_3;

  const VariableValue & _rho_edge_4;

  const VariableValue & _rho_screw_1;

  const VariableValue & _rho_screw_2;

  const VariableValue & _rho_screw_3;

  const VariableValue & _rho_screw_4;

  /// Determine which quandrant
  const enum class QuandrantType { qua_one, qua_two, qua_thr, qua_fou } _quandrant;

  // Character of dislocations (edge or screw)
  const enum class DisloCharacter { edge, screw } _dislo_character;
};

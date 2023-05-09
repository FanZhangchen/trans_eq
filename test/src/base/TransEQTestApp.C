//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "TransEQTestApp.h"
#include "TransEQApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
TransEQTestApp::validParams()
{
  InputParameters params = TransEQApp::validParams();
  return params;
}

TransEQTestApp::TransEQTestApp(InputParameters parameters) : MooseApp(parameters)
{
  TransEQTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

TransEQTestApp::~TransEQTestApp() {}

void
TransEQTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  TransEQApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"TransEQTestApp"});
    Registry::registerActionsTo(af, {"TransEQTestApp"});
  }
}

void
TransEQTestApp::registerApps()
{
  registerApp(TransEQApp);
  registerApp(TransEQTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
TransEQTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  TransEQTestApp::registerAll(f, af, s);
}
extern "C" void
TransEQTestApp__registerApps()
{
  TransEQTestApp::registerApps();
}

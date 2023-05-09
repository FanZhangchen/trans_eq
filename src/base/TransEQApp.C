#include "TransEQApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
TransEQApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

TransEQApp::TransEQApp(InputParameters parameters) : MooseApp(parameters)
{
  TransEQApp::registerAll(_factory, _action_factory, _syntax);
}

TransEQApp::~TransEQApp() {}

void
TransEQApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"TransEQApp"});
  Registry::registerActionsTo(af, {"TransEQApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
TransEQApp::registerApps()
{
  registerApp(TransEQApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
TransEQApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  TransEQApp::registerAll(f, af, s);
}
extern "C" void
TransEQApp__registerApps()
{
  TransEQApp::registerApps();
}

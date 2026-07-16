#include "palmApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
palmApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

palmApp::palmApp(InputParameters parameters) : MooseApp(parameters)
{
  palmApp::registerAll(_factory, _action_factory, _syntax);
}

palmApp::~palmApp() {}

void
palmApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<palmApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"palmApp"});
  Registry::registerActionsTo(af, {"palmApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
palmApp::registerApps()
{
  registerApp(palmApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
palmApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  palmApp::registerAll(f, af, s);
}
extern "C" void
palmApp__registerApps()
{
  palmApp::registerApps();
}

//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "palmTestApp.h"
#include "palmApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
palmTestApp::validParams()
{
  InputParameters params = palmApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

palmTestApp::palmTestApp(InputParameters parameters) : MooseApp(parameters)
{
  palmTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

palmTestApp::~palmTestApp() {}

void
palmTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  palmApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"palmTestApp"});
    Registry::registerActionsTo(af, {"palmTestApp"});
  }
}

void
palmTestApp::registerApps()
{
  registerApp(palmApp);
  registerApp(palmTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
palmTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  palmTestApp::registerAll(f, af, s);
}
extern "C" void
palmTestApp__registerApps()
{
  palmTestApp::registerApps();
}

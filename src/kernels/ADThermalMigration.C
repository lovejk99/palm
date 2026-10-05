#include "ADThermalMigration.h"

registerMooseObject("palmApp", ADThermalMigration);

InputParameters
ADThermalMigration::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription(
      "Thermal-migration term -div(M_T grad T / T) of the split Cahn-Hilliard equation, "
      "discretised as R = -int (M_T/T) grad T . grad psi, i.e. a flux J_T = +M_T grad T / T. "
      "The thermal mobility is read from an AD material property so the full Jacobian is "
      "obtained by automatic differentiation.");
  params.addRequiredCoupledVar("T", "Temperature variable (nonlinear or auxiliary)");
  params.addParam<MaterialPropertyName>(
      "thermal_mobility", "thermal_mobility", "AD material property holding M_T");
  return params;
}

ADThermalMigration::ADThermalMigration(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _T(adCoupledValue("T")),
    _grad_T(adCoupledGradient("T")),
    _MT(getADMaterialProperty<Real>("thermal_mobility"))
{
}

ADRealVectorValue
ADThermalMigration::precomputeQpResidual()
{
  return -_MT[_qp] / _T[_qp] * _grad_T[_qp];
}

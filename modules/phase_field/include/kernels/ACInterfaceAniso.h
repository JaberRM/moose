/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ACINTERFACEANISO_H
#define ACINTERFACEANISO_H

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

class ACInterfaceAniso;

template <>
InputParameters validParams<ACInterfaceAniso>();

/**
 * Compute the Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 */
class ACInterfaceAniso : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  ACInterfaceAniso(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  RealGradient gradL();
  RealGradient gradKappa();

  /// the \f$ \kappa\nabla(L\psi) \f$ term
  RealGradient kappaNablaLPsi();

  /// Mobility
  const MaterialProperty<Real> & _L;
  /// Interfacial parameter
  const MaterialProperty<RealTensorValue> & _kappa;

  /// flag set if L is a function of non-linear variables in args
  const bool _variable_L;

  /// @{ Mobility derivatives w.r.t. order parameter
  const MaterialProperty<Real> & _dLdop;
  const MaterialProperty<Real> & _d2Ldop2;
  /// @}

  /// kappa derivative w.r.t. order parameter
  const MaterialProperty<RealTensorValue> & _dkappadop;

  /// number of coupled variables
  const unsigned int _nvar;

  /// @{ Mobility derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dLdarg;
  std::vector<const MaterialProperty<Real> *> _d2Ldargdop;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2Ldarg2;
  /// @}

  /// kappa derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<RealTensorValue> *> _dkappadarg;

  /// Gradients for all coupled variables
  std::vector<const VariableGradient *> _gradarg;
};

#endif // ACINTERFACE_H

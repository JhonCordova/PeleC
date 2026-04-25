#ifndef PELEC_ODTLES_ODTDIFFUSION_H
#define PELEC_ODTLES_ODTDIFFUSION_H

#include <array>
#include <vector>

#include <AMReX_REAL.H>

#include "ODTLineGeometry.H"
#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTDiffusion
{
public:
  enum class MomentumMolecularForm : int
  {
    // Internal debug/test-only legacy alternative.
    ConservativeVariableNuOnRhoU = 0,
    // Production molecular momentum model.
    VelocityGradientMuFlux = 1
  };

  enum class Component : int
  {
    Rho = 0,
    RhoU = 1,
    RhoV = 2,
    RhoW = 3,
    RhoE = 4,
    NumComponents = 5
  };

  struct Controls
  {
    // MVP diffusion target: momentum components only by default.
    std::array<bool, static_cast<int>(Component::NumComponents)>
      diffuse_component{false, true, true, true, false};

    // Production default is to use host molecular viscosity for momentum.
    bool use_molecular_viscosity_for_momentum = true;
    // Select mathematical momentum operator when molecular mode is enabled.
    MomentumMolecularForm momentum_molecular_form =
      MomentumMolecularForm::VelocityGradientMuFlux;
    // Optional lower bound for molecular viscosity used by momentum operator.
    // Default 0 preserves strict host-transport value.
    amrex::Real molecular_viscosity_floor = 0.0;
    bool fail_on_molecular_viscosity_recovery_failure = true;

    // Constant line diffusivity per selected component.
    std::array<amrex::Real, static_cast<int>(Component::NumComponents)>
      diffusivity{0.0, 0.0, 0.0, 0.0, 0.0};
  };

  struct MomentumViscosityProfile
  {
    bool success = false;
    int failure_cell = -1;
    std::vector<amrex::Real> dynamic_viscosity_mu;
  };

  struct ValidationReport
  {
    bool uniform_invariant = false;
    bool nonuniform_evolved = false;
    amrex::Real uniform_max_abs_change = 0.0;
    amrex::Real nonuniform_max_abs_change = 0.0;
  };

  // Implicit 1D diffusion update on fixed, uniform line support.
  // For momentum components when molecular mode is enabled, the selected form
  // in Controls::momentum_molecular_form is applied:
  //   (A) d(rho*u_i)/dt = d/ds ( nu_eff * d(rho*u_i)/ds ), nu_eff = mu/rho
  //       (internal debug/test legacy form)
  //   (B) d(rho*u_i)/dt = d/ds ( mu * d(u_i)/ds ), u_i=(rho*u_i)/rho
  //       (production runtime form)
  // with mu recovered from host transport using local thermochemical state.
  // MVP boundary closure: zero-gradient (Neumann-like) treatment at line ends.
  // This operator is intentionally local to ODT line evolution and does not
  // perform PeleC SGS flux deposition.
  static void applyImplicitUniform(
    const ODTLineGeometry& geom,
    ODTLineState& line_state,
    amrex::Real dt,
    const Controls& controls);

  // Build local molecular viscosity profile used by momentum diffusion:
  // - mu from host transport at recovered (rho,T,Y)
  // Returns success=false if thermochemical recovery fails on any line cell.
  static MomentumViscosityProfile buildMomentumViscosityProfile(
    const ODTLineState& line_state);

  // Lightweight validation hook for MVP diffusion behavior:
  // - uniform profile invariance (up to roundoff)
  // - nonuniform profile evolution under diffusion
  static ValidationReport runMVPValidationHook();
};

} // namespace pelec::odtles

#endif

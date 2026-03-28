#ifndef PELEC_ODTLES_ODTDIFFUSION_H
#define PELEC_ODTLES_ODTDIFFUSION_H

#include <array>

#include <AMReX_REAL.H>

#include "ODTLineGeometry.H"
#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTDiffusion
{
public:
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

    // Constant line diffusivity per selected component.
    std::array<amrex::Real, static_cast<int>(Component::NumComponents)>
      diffusivity{0.0, 0.0, 0.0, 0.0, 0.0};
  };

  struct ValidationReport
  {
    bool uniform_invariant = false;
    bool nonuniform_evolved = false;
    amrex::Real uniform_max_abs_change = 0.0;
    amrex::Real nonuniform_max_abs_change = 0.0;
  };

  // Implicit 1D diffusion update on fixed, uniform line support.
  // MVP boundary closure: zero-gradient (Neumann-like) treatment at line ends.
  // This operator is intentionally local to ODT line evolution and does not
  // perform PeleC SGS flux deposition.
  static void applyImplicitUniform(
    const ODTLineGeometry& geom,
    ODTLineState& line_state,
    amrex::Real dt,
    const Controls& controls);

  // Lightweight validation hook for MVP diffusion behavior:
  // - uniform profile invariance (up to roundoff)
  // - nonuniform profile evolution under diffusion
  static ValidationReport runMVPValidationHook();
};

} // namespace pelec::odtles

#endif

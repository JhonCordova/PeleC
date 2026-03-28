#ifndef PELEC_ODTLES_ODTKERNEL_H
#define PELEC_ODTLES_ODTKERNEL_H

#include <vector>

#include <AMReX_REAL.H>

#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTKernel
{
public:
  struct VelocityKernelSpec
  {
    bool enabled = false;
    amrex::Real amp_u = 0.0;
    amrex::Real amp_v = 0.0;
    amrex::Real amp_w = 0.0;
    // Optional shape over the mapped interval. If empty, use uniform shape.
    std::vector<amrex::Real> shape;
  };

  // Apply a velocity-only redistribution hook over [i_lo, i_hi] with a
  // discrete zero-net momentum contribution per velocity component.
  // Phase scope: rho and rhoE are intentionally unchanged in this MVP hook.
  // Full variable-density/compressible kernel-energy coupling is deferred.
  static bool applyVelocityKernelZeroNet(
    ODTLineState& line_state,
    int i_lo,
    int i_hi,
    const VelocityKernelSpec& kernel);
};

} // namespace pelec::odtles

#endif

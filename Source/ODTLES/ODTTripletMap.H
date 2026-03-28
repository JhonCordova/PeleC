#ifndef PELEC_ODTLES_ODTTRIPLETMAP_H
#define PELEC_ODTLES_ODTTRIPLETMAP_H

#include "ODTKernel.H"
#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTTripletMap
{
public:
  struct EddyInterval
  {
    int i_lo = -1;
    int i_hi = -1;

    int size() const noexcept { return i_hi - i_lo + 1; }
    bool valid() const noexcept { return i_lo >= 0 && i_hi >= i_lo; }
  };

  struct ApplyReport
  {
    bool applied = false;
    bool interval_multiple_of_three = false;
    bool kernel_applied = false;
  };

  struct ValidationReport
  {
    bool uniform_invariant = false;
    bool nonuniform_modified = false;
    bool map_locality_preserved = false;
    bool kernel_zero_net_momentum = false;
    amrex::Real uniform_max_abs_change = 0.0;
    amrex::Real nonuniform_max_abs_change = 0.0;
    amrex::Real outside_interval_max_abs_change = 0.0;
    amrex::Real kernel_net_rhou_change = 0.0;
    amrex::Real kernel_net_rhov_change = 0.0;
    amrex::Real kernel_net_rhow_change = 0.0;
  };

  // MVP interval policy: apply only when interval length n is divisible by 3.
  static bool canApplyPermutation(
    const ODTLineState& line_state,
    const EddyInterval& interval);

  // Apply discrete triplet-map permutation on the selected interval only.
  // Outside [i_lo, i_hi], the line state is unchanged.
  static ApplyReport applyPermutation(
    ODTLineState& line_state,
    const EddyInterval& interval);

  // Same permutation with optional separated velocity-kernel hook.
  static ApplyReport applyPermutationWithKernel(
    ODTLineState& line_state,
    const EddyInterval& interval,
    const ODTKernel::VelocityKernelSpec& kernel);

  // Lightweight T5.2 validation hook for map/kernel invariants.
  static ValidationReport runMVPValidationHook();
};

} // namespace pelec::odtles

#endif

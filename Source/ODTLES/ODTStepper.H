#ifndef PELEC_ODTLES_ODTSTEPPER_H
#define PELEC_ODTLES_ODTSTEPPER_H

#include <cstdint>
#include <vector>

#include <AMReX_REAL.H>

#include "ODTDiffusion.H"
#include "ODTEventSampler.H"
#include "ODTKernel.H"
#include "ODTLineGeometry.H"
#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTStepper
{
public:
  enum class IncompatibleIntervalPolicy : int
  {
    AdjustToCompatibleElseReject = 0,
    RejectAndContinue = 1
  };

  struct Controls
  {
    ODTEventSampler::Controls sampler_controls{};
    ODTDiffusion::Controls diffusion_controls{};
    ODTKernel::VelocityKernelSpec kernel_spec{};
    IncompatibleIntervalPolicy incompatible_interval_policy =
      IncompatibleIntervalPolicy::AdjustToCompatibleElseReject;
    int max_internal_iterations = 100000;
  };

  struct SampleTrace
  {
    amrex::Real t_before = 0.0;
    amrex::Real wait_time = 0.0;
    amrex::Real sampled_interval_probability = 0.0;
    amrex::Real sampled_interval_physical_length = 0.0;
    amrex::Real state_energetic_admissibility = 0.0;
    amrex::Real state_turnover_time = 0.0;
    amrex::Real state_event_hazard = 0.0;
    amrex::Real state_acceptance_probability = 0.0;
    amrex::Real acceptance_draw = 1.0;
    ODTTripletMap::EddyInterval sampled_interval{};
    ODTTripletMap::EddyInterval applied_interval{};
    bool compatible_triplet_interval = false; // compatibility of sampled interval
    bool adjusted_to_compatible = false;
    bool accepted_by_state_model = false;
    bool rejected = false;
    bool applied = false;
  };

  struct StepReport
  {
    amrex::Real accumulated_internal_time = 0.0;
    int attempted_events = 0;
    int adjusted_events = 0;
    int resampled_events = 0;
    int applied_events = 0;
    int rejected_incompatible_events = 0;
    int rejected_state_model_events = 0;
    int rejected_events = 0;
    int diffusion_substeps = 0;
    bool closed_by_diffusion_only_catchup = false;
    bool reached_dt_les = false;
    bool deterministic_sampler = false;
    std::uint64_t sampler_seed = 0;
    std::vector<SampleTrace> trace;
  };

  struct ValidationReport
  {
    bool closes_exact_dt = false;
    bool no_overshoot = false;
    bool no_event_path_handled = false;
    bool event_path_attempted = false;
    bool compatibility_policy_exercised = false;
    bool deterministic_reproducible = false;
    bool metadata_consistent = false;
    amrex::Real closure_error = 0.0;
    amrex::Real overshoot_amount = 0.0;
  };

  // Advance one local ODT line over exactly one host LES step dt_les through
  // internal subcycling. This never modifies host timestep ownership.
  static StepReport advanceOneLESTimestep(
    const ODTLineGeometry& geom,
    ODTLineState& line_state,
    amrex::Real dt_les,
    const Controls& controls);

  // Lightweight T5.3 validation hook for local stepper invariants.
  static ValidationReport runMVPValidationHook();
};

} // namespace pelec::odtles

#endif

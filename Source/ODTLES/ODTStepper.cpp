#include "ODTStepper.H"

#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

#include "ODTTripletMap.H"

namespace pelec::odtles
{

namespace
{

bool
adjustIntervalToTripletCompatible(
  const ODTTripletMap::EddyInterval& sampled,
  int line_cells,
  ODTTripletMap::EddyInterval& adjusted)
{
  adjusted = {};
  if (!sampled.valid() || line_cells < 3) {
    return false;
  }

  int n_target = sampled.size();
  if (n_target >= 3) {
    n_target = 3 * (n_target / 3);
  } else {
    n_target = 3;
  }
  const int max_compatible_n = 3 * (line_cells / 3);
  if (max_compatible_n < 3) {
    return false;
  }
  n_target = std::min(n_target, max_compatible_n);
  if (n_target < 3 || (n_target % 3) != 0) {
    return false;
  }

  const amrex::Real center =
    0.5 * static_cast<amrex::Real>(sampled.i_lo + sampled.i_hi);
  int i_lo = static_cast<int>(std::round(center - 0.5 * (n_target - 1)));
  i_lo = std::max(0, std::min(i_lo, line_cells - n_target));

  adjusted.i_lo = i_lo;
  adjusted.i_hi = i_lo + n_target - 1;
  return adjusted.valid() && adjusted.i_hi < line_cells;
}

bool
statesEqual(const ODTLineState& a, const ODTLineState& b, amrex::Real tol)
{
  if (a.numCells() != b.numCells()) {
    return false;
  }
  for (int i = 0; i < a.numCells(); ++i) {
    const auto& ca = a.cell(i);
    const auto& cb = b.cell(i);
    if (std::abs(ca.rho - cb.rho) > tol) {
      return false;
    }
    if (std::abs(ca.rhou - cb.rhou) > tol) {
      return false;
    }
    if (std::abs(ca.rhov - cb.rhov) > tol) {
      return false;
    }
    if (std::abs(ca.rhow - cb.rhow) > tol) {
      return false;
    }
    if (std::abs(ca.rhoE - cb.rhoE) > tol) {
      return false;
    }
  }
  return true;
}

} // namespace

ODTStepper::StepReport
ODTStepper::advanceOneLESTimestep(
  const ODTLineGeometry& geom,
  ODTLineState& line_state,
  amrex::Real dt_les,
  const Controls& controls)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTStepper requires defined line geometry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTStepper requires initialized, valid line state");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.numCells() == geom.supportCellCount(),
    "ODTStepper requires state/geometry support-size consistency");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt_les >= 0.0, "ODTStepper requires dt_les >= 0");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    controls.max_internal_iterations > 0,
    "ODTStepper requires positive max_internal_iterations");

  StepReport report{};
  report.deterministic_sampler = controls.sampler_controls.deterministic;
  report.sampler_seed = controls.sampler_controls.seed;

  if (dt_les == 0.0) {
    report.accumulated_internal_time = 0.0;
    report.reached_dt_les = true;
    return report;
  }

  ODTEventSampler sampler(controls.sampler_controls);

  const amrex::Real tol =
    64.0 * std::numeric_limits<amrex::Real>::epsilon() *
    std::max<amrex::Real>(1.0, dt_les);

  amrex::Real t = 0.0;
  int iter = 0;
  while (t < dt_les - tol) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      iter < controls.max_internal_iterations,
      "ODTStepper exceeded max_internal_iterations before closing dt_les");
    ++iter;

    const amrex::Real remaining = dt_les - t;
    AMREX_ALWAYS_ASSERT(remaining > 0.0);

    const auto sample = sampler.sample(line_state.numCells());
    if (!sample.valid) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, remaining, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t = dt_les;
      report.closed_by_diffusion_only_catchup = true;
      break;
    }

    const amrex::Real wait_time = std::max<amrex::Real>(0.0, sample.wait_time);
    if (wait_time >= remaining - tol) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, remaining, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t = dt_les;
      report.closed_by_diffusion_only_catchup = true;
      break;
    }

    if (wait_time > 0.0) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, wait_time, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t += wait_time;
    }

    report.attempted_events += 1;
    SampleTrace tr{};
    tr.t_before = t;
    tr.wait_time = wait_time;
    tr.sampled_interval = sample.interval;
    tr.applied_interval = sample.interval;
    tr.compatible_triplet_interval =
      ODTTripletMap::canApplyPermutation(line_state, sample.interval);

    if (!tr.compatible_triplet_interval) {
      switch (controls.incompatible_interval_policy) {
      case IncompatibleIntervalPolicy::AdjustToCompatibleElseReject: {
        ODTTripletMap::EddyInterval adjusted{};
        if (
          adjustIntervalToTripletCompatible(
            sample.interval, line_state.numCells(), adjusted) &&
          ODTTripletMap::canApplyPermutation(line_state, adjusted)) {
          tr.adjusted_to_compatible = true;
          tr.applied_interval = adjusted;
          report.adjusted_events += 1;
        } else {
          report.rejected_incompatible_events += 1;
          report.rejected_events += 1;
          tr.rejected = true;
          report.trace.push_back(tr);
          continue;
        }
        break;
      }
      case IncompatibleIntervalPolicy::RejectAndContinue:
      default:
        report.rejected_incompatible_events += 1;
        report.rejected_events += 1;
        tr.rejected = true;
        report.trace.push_back(tr);
        continue;
      }
    }

    switch (controls.incompatible_interval_policy) {
    case IncompatibleIntervalPolicy::AdjustToCompatibleElseReject:
    case IncompatibleIntervalPolicy::RejectAndContinue: {
      const auto rep = ODTTripletMap::applyPermutationWithKernel(
        line_state, tr.applied_interval, controls.kernel_spec);
      tr.applied = rep.applied;
      if (rep.applied) {
        report.applied_events += 1;
      } else {
        tr.rejected = true;
        report.rejected_events += 1;
      }
      break;
    }
    default:
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "Unsupported interval policy");
    }

    report.trace.push_back(tr);
  }

  if (t < dt_les) {
    const amrex::Real remaining = dt_les - t;
    if (remaining > 0.0) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, remaining, controls.diffusion_controls);
      report.diffusion_substeps += 1;
    }
    t = dt_les;
  }

  if (std::abs(t - dt_les) <= tol) {
    t = dt_les;
  }

  report.accumulated_internal_time = t;
  report.reached_dt_les = (std::abs(t - dt_les) <= tol);
  return report;
}

ODTStepper::ValidationReport
ODTStepper::runMVPValidationHook()
{
  ValidationReport out{};

  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(2, 0, 0)));
  amrex::RealBox rb(
    AMREX_D_DECL(0.0, 0.0, 0.0), AMREX_D_DECL(3.0, 1.0, 1.0));
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(1, 0, 0));
  const ODTLineGeometry line_geom(geom, 0, owner, 0);

  ODTLineState base_state;
  base_state.initialize(line_geom);
  for (int i = 0; i < base_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0;
    c.rhou = static_cast<amrex::Real>(i + 1);
    c.rhov = 0.1 * static_cast<amrex::Real>(i);
    c.rhow = -0.2 * static_cast<amrex::Real>(i);
    // Keep manufactured conservative state thermochemically admissible for
    // the production diffusion path (thermochemical recovery + molecular mu).
    const amrex::Real kinetic =
      0.5 * (c.rhou * c.rhou + c.rhov * c.rhov + c.rhow * c.rhow) / c.rho;
    const amrex::Real eint = 5.0;
    c.rhoE = c.rho * eint + kinetic;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      c.rhoY[static_cast<std::size_t>(n)] =
        c.rho / static_cast<amrex::Real>(NUM_SPECIES);
    }
    base_state.setCell(i, c);
  }
  base_state.setValid(true);

  constexpr amrex::Real dt_les = 1.0e-2;
  constexpr amrex::Real tol = 1.0e-12;

  // No-event path: event_rate = 0 forces diffusion-only catch-up closure.
  Controls no_event_ctrl{};
  no_event_ctrl.sampler_controls.event_rate = 0.0;
  no_event_ctrl.sampler_controls.deterministic = true;
  no_event_ctrl.sampler_controls.seed = 11ULL;
  ODTLineState s_noevent = base_state;
  const auto no_event_rep =
    advanceOneLESTimestep(line_geom, s_noevent, dt_les, no_event_ctrl);
  out.no_event_path_handled =
    no_event_rep.closed_by_diffusion_only_catchup &&
    no_event_rep.attempted_events == 0 && no_event_rep.reached_dt_les;

  // Event-enabled path with deterministic seed.
  Controls event_ctrl{};
  event_ctrl.sampler_controls.event_rate = 1.0e3;
  event_ctrl.sampler_controls.min_interval_size = 2;
  event_ctrl.sampler_controls.deterministic = true;
  event_ctrl.sampler_controls.seed = 17ULL;
  event_ctrl.incompatible_interval_policy =
    IncompatibleIntervalPolicy::AdjustToCompatibleElseReject;

  ODTLineState s_event_a = base_state;
  ODTLineState s_event_b = base_state;
  const auto event_rep_a =
    advanceOneLESTimestep(line_geom, s_event_a, dt_les, event_ctrl);
  const auto event_rep_b =
    advanceOneLESTimestep(line_geom, s_event_b, dt_les, event_ctrl);

  out.event_path_attempted = event_rep_a.attempted_events > 0;
  out.compatibility_policy_exercised =
    (event_rep_a.adjusted_events > 0 || event_rep_a.rejected_incompatible_events > 0);

  const bool same_reports =
    event_rep_a.attempted_events == event_rep_b.attempted_events &&
    event_rep_a.adjusted_events == event_rep_b.adjusted_events &&
    event_rep_a.applied_events == event_rep_b.applied_events &&
    event_rep_a.rejected_incompatible_events ==
      event_rep_b.rejected_incompatible_events &&
    event_rep_a.rejected_events == event_rep_b.rejected_events &&
    event_rep_a.closed_by_diffusion_only_catchup ==
      event_rep_b.closed_by_diffusion_only_catchup &&
    std::abs(
      event_rep_a.accumulated_internal_time - event_rep_b.accumulated_internal_time) <=
      tol;
  out.deterministic_reproducible =
    same_reports && statesEqual(s_event_a, s_event_b, tol);

  const bool metadata_a_ok =
    event_rep_a.attempted_events >= 0 && event_rep_a.applied_events >= 0 &&
    event_rep_a.rejected_events >= 0 &&
    event_rep_a.applied_events + event_rep_a.rejected_events ==
      event_rep_a.attempted_events &&
    event_rep_a.adjusted_events <= event_rep_a.attempted_events &&
    event_rep_a.rejected_incompatible_events <= event_rep_a.rejected_events &&
    static_cast<int>(event_rep_a.trace.size()) == event_rep_a.attempted_events;
  out.metadata_consistent = metadata_a_ok;

  out.closure_error = std::abs(event_rep_a.accumulated_internal_time - dt_les);
  out.overshoot_amount =
    std::max<amrex::Real>(0.0, event_rep_a.accumulated_internal_time - dt_les);
  out.closes_exact_dt = event_rep_a.reached_dt_les && out.closure_error <= tol;
  out.no_overshoot = out.overshoot_amount <= tol;

  return out;
}

} // namespace pelec::odtles

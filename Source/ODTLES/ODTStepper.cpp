#include "ODTStepper.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

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

struct LocalEventModel
{
  amrex::Real energetic_admissibility = 0.0;
  amrex::Real turnover_time = 0.0;
  amrex::Real event_hazard = 0.0;
  amrex::Real acceptance_probability = 0.0;
};

struct IntervalProposal
{
  ODTTripletMap::EddyInterval interval{};
  amrex::Real probability = 0.0;
  bool valid = false;
};

amrex::Real
intervalPhysicalLength(
  const ODTLineGeometry& geom, const ODTTripletMap::EddyInterval& interval)
{
  if (!interval.valid()) {
    return 0.0;
  }
  const auto& segs = geom.supportCells();
  const int lo = std::max(0, interval.i_lo);
  const int hi = std::min(static_cast<int>(segs.size()) - 1, interval.i_hi);
  if (lo > hi) {
    return 0.0;
  }
  amrex::Real len = 0.0;
  for (int i = lo; i <= hi; ++i) {
    const auto& seg = segs[static_cast<std::size_t>(i)];
    len += seg.s_interval.hi - seg.s_interval.lo;
  }
  return std::max<amrex::Real>(0.0, len);
}

amrex::Real
supportPhysicalLength(const ODTLineGeometry& geom)
{
  const auto& segs = geom.supportCells();
  amrex::Real len = 0.0;
  for (const auto& seg : segs) {
    len += seg.s_interval.hi - seg.s_interval.lo;
  }
  return std::max<amrex::Real>(0.0, len);
}

IntervalProposal
sampleIntervalProposal(
  const ODTLineGeometry& geom,
  int line_cells,
  int min_interval_size,
  ODTEventSampler& sampler)
{
  IntervalProposal out{};
  if (line_cells <= 0) {
    return out;
  }
  const int min_n = std::max(3, min_interval_size);
  if (line_cells < min_n) {
    return out;
  }

  std::vector<int> admissible_sizes;
  std::vector<amrex::Real> weights;
  admissible_sizes.reserve(static_cast<std::size_t>(line_cells));
  weights.reserve(static_cast<std::size_t>(line_cells));
  amrex::Real w_sum = 0.0;
  for (int n = min_n; n <= line_cells; ++n) {
    if ((n % 3) != 0) {
      continue;
    }
    const int n_starts = line_cells - n + 1;
    if (n_starts <= 0) {
      continue;
    }
    ODTTripletMap::EddyInterval rep{};
    rep.i_lo = 0;
    rep.i_hi = n - 1;
    const amrex::Real length = intervalPhysicalLength(geom, rep);
    if (length <= 0.0) {
      continue;
    }

    // Proposal weights are inverse physical size with multiplicity by start
    // positions, so q(interval) is physically scaled and normalized.
    const amrex::Real w =
      static_cast<amrex::Real>(n_starts) / length;
    admissible_sizes.push_back(n);
    weights.push_back(w);
    w_sum += w;
  }
  if (admissible_sizes.empty() || w_sum <= 0.0) {
    return out;
  }

  const amrex::Real u = sampler.sampleUnitUniform() * w_sum;
  amrex::Real accum = 0.0;
  int n = admissible_sizes.back();
  std::size_t selected_k = weights.size() - 1;
  for (std::size_t k = 0; k < weights.size(); ++k) {
    accum += weights[k];
    if (u <= accum || k + 1 == weights.size()) {
      n = admissible_sizes[k];
      selected_k = k;
      break;
    }
  }

  const int n_starts = line_cells - n + 1;
  if (n_starts <= 0) {
    return out;
  }
  const int i_lo =
    std::min(n_starts - 1, static_cast<int>(sampler.sampleUnitUniform() * n_starts));

  out.interval.i_lo = i_lo;
  out.interval.i_hi = i_lo + n - 1;
  const amrex::Real p_n = weights[selected_k] / w_sum;
  const amrex::Real p_i = 1.0 / static_cast<amrex::Real>(n_starts);
  out.probability = p_n * p_i;
  out.valid = out.interval.valid();
  return out;
}

LocalEventModel
stateDependentEventModel(
  const ODTLineGeometry& geom,
  const ODTLineState& line_state,
  const ODTTripletMap::EddyInterval& interval)
{
  LocalEventModel out{};
  if (!interval.valid()) {
    return out;
  }

  const int i_lo = std::max(0, interval.i_lo);
  const int i_hi = std::min(line_state.numCells() - 1, interval.i_hi);
  if (i_lo > i_hi) {
    return out;
  }
  const auto& segs = geom.supportCells();
  if (static_cast<int>(segs.size()) != line_state.numCells()) {
    return out;
  }

  amrex::Real rho_sum = 0.0;
  amrex::Real rhou_sum = 0.0;
  amrex::Real rhov_sum = 0.0;
  amrex::Real rhow_sum = 0.0;
  for (int i = i_lo; i <= i_hi; ++i) {
    const auto& c = line_state.cell(i);
    if (c.rho <= 0.0) {
      return out;
    }
    rho_sum += c.rho;
    rhou_sum += c.rhou;
    rhov_sum += c.rhov;
    rhow_sum += c.rhow;
  }
  if (rho_sum <= 0.0) {
    return out;
  }

  const amrex::Real inv_rho_sum = 1.0 / rho_sum;
  const amrex::Real u_mean = rhou_sum * inv_rho_sum;
  const amrex::Real v_mean = rhov_sum * inv_rho_sum;
  const amrex::Real w_mean = rhow_sum * inv_rho_sum;

  amrex::Real var_mass_weighted = 0.0;
  for (int i = i_lo; i <= i_hi; ++i) {
    const auto& c = line_state.cell(i);
    const amrex::Real u = c.rhou / c.rho;
    const amrex::Real v = c.rhov / c.rho;
    const amrex::Real w = c.rhow / c.rho;
    const amrex::Real du = u - u_mean;
    const amrex::Real dv = v - v_mean;
    const amrex::Real dw = w - w_mean;
    var_mass_weighted += c.rho * (du * du + dv * dv + dw * dw);
  }
  var_mass_weighted *= inv_rho_sum;

  // Build a translation-invariant internal-gradient scale from velocity
  // differences across neighboring support cells in the candidate interval.
  // Scale by interval length so this term stays physically scaled (u^2) and
  // does not depend on raw segment count.
  amrex::Real grad_sq_sum = 0.0;
  amrex::Real grad_w_sum = 0.0;
  amrex::Real interval_length = 0.0;
  for (int i = i_lo; i <= i_hi; ++i) {
    const auto& seg = segs[static_cast<std::size_t>(i)];
    interval_length += seg.s_interval.hi - seg.s_interval.lo;
  }
  interval_length = std::max<amrex::Real>(0.0, interval_length);

  for (int i = i_lo; i < i_hi; ++i) {
    const auto& c_l = line_state.cell(i);
    const auto& c_r = line_state.cell(i + 1);
    const auto& seg_l = segs[static_cast<std::size_t>(i)];
    const auto& seg_r = segs[static_cast<std::size_t>(i + 1)];
    const amrex::Real ds = std::abs(seg_r.s_center - seg_l.s_center);
    if (ds <= 0.0) {
      continue;
    }
    const amrex::Real u_l = c_l.rhou / c_l.rho;
    const amrex::Real v_l = c_l.rhov / c_l.rho;
    const amrex::Real w_l = c_l.rhow / c_l.rho;
    const amrex::Real u_r = c_r.rhou / c_r.rho;
    const amrex::Real v_r = c_r.rhov / c_r.rho;
    const amrex::Real w_r = c_r.rhow / c_r.rho;
    const amrex::Real du = (u_r - u_l) / ds;
    const amrex::Real dv = (v_r - v_l) / ds;
    const amrex::Real dw = (w_r - w_l) / ds;
    const amrex::Real w_face = 0.5 * (c_l.rho + c_r.rho);
    grad_sq_sum += w_face * (du * du + dv * dv + dw * dw);
    grad_w_sum += w_face;
  }
  const amrex::Real grad_rate_mass_weighted =
    (grad_w_sum > 0.0) ? (grad_sq_sum / grad_w_sum) : 0.0;
  const amrex::Real grad_mass_weighted =
    interval_length * interval_length * grad_rate_mass_weighted;

  const amrex::Real eps = 1.0e-30;
  out.energetic_admissibility =
    grad_mass_weighted / (grad_mass_weighted + var_mass_weighted + eps);
  out.energetic_admissibility = std::max<amrex::Real>(
    0.0, std::min<amrex::Real>(1.0, out.energetic_admissibility));

  const amrex::Real u_char = std::sqrt(std::max<amrex::Real>(
    0.0, grad_mass_weighted + var_mass_weighted));
  out.turnover_time = interval_length / (u_char + eps);
  out.turnover_time = std::max<amrex::Real>(0.0, out.turnover_time);

  // Local accepted-event hazard [1/time], built from state alone.
  out.event_hazard =
    out.energetic_admissibility / (out.turnover_time + eps);
  out.event_hazard = std::max<amrex::Real>(0.0, out.event_hazard);
  return out;
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

    const amrex::Real event_rate =
      std::max<amrex::Real>(0.0, controls.sampler_controls.event_rate);
    const amrex::Real wait_time = sampler.sampleExponentialWait(event_rate);
    if (!std::isfinite(wait_time)) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, remaining, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t = dt_les;
      report.closed_by_diffusion_only_catchup = true;
      break;
    }

    const amrex::Real wait_time_clamped = std::max<amrex::Real>(0.0, wait_time);
    if (wait_time_clamped >= remaining - tol) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, remaining, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t = dt_les;
      report.closed_by_diffusion_only_catchup = true;
      break;
    }

    if (wait_time_clamped > 0.0) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, wait_time_clamped, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t += wait_time_clamped;
    }

    const auto proposal = sampleIntervalProposal(
      geom, line_state.numCells(), controls.sampler_controls.min_interval_size,
      sampler);
    if (!proposal.valid) {
      ODTDiffusion::applyImplicitUniform(
        geom, line_state, dt_les - t, controls.diffusion_controls);
      report.diffusion_substeps += 1;
      t = dt_les;
      report.closed_by_diffusion_only_catchup = true;
      break;
    }

    report.attempted_events += 1;
    SampleTrace tr{};
    tr.t_before = t;
    tr.wait_time = wait_time_clamped;
    tr.sampled_interval = proposal.interval;
    tr.applied_interval = proposal.interval;
    tr.sampled_interval_probability = proposal.probability;
    tr.sampled_interval_physical_length =
      intervalPhysicalLength(geom, proposal.interval);
    tr.applied_interval_physical_length = tr.sampled_interval_physical_length;
    tr.compatible_triplet_interval =
      ODTTripletMap::canApplyPermutation(line_state, proposal.interval);

    if (!tr.compatible_triplet_interval) {
      switch (controls.incompatible_interval_policy) {
      case IncompatibleIntervalPolicy::AdjustToCompatibleElseReject: {
        ODTTripletMap::EddyInterval adjusted{};
        if (
          adjustIntervalToTripletCompatible(
            proposal.interval, line_state.numCells(), adjusted) &&
          ODTTripletMap::canApplyPermutation(line_state, adjusted)) {
          tr.adjusted_to_compatible = true;
          tr.applied_interval = adjusted;
          tr.applied_interval_physical_length =
            intervalPhysicalLength(geom, adjusted);
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
      const auto local_model =
        stateDependentEventModel(geom, line_state, tr.applied_interval);
      tr.state_energetic_admissibility = local_model.energetic_admissibility;
      tr.state_turnover_time = local_model.turnover_time;
      const amrex::Real support_len = supportPhysicalLength(geom);
      const amrex::Real interval_len_frac =
        (support_len > 0.0)
          ? std::max<amrex::Real>(
              0.0, tr.applied_interval_physical_length / support_len)
          : 0.0;
      const amrex::Real sampled_interval_prob =
        std::max<amrex::Real>(0.0, tr.sampled_interval_probability);
      // The local state model provides an intrinsic interval turnover frequency
      // [1/time]. Convert this to a discrete interval-occurrence target rate by
      // accounting for proposal probability mass and physical support fraction.
      tr.state_event_hazard =
        local_model.event_hazard * sampled_interval_prob * interval_len_frac;
      const amrex::Real proposal_rate_density = std::max<amrex::Real>(
        0.0, event_rate * sampled_interval_prob);
      tr.state_proposal_rate_density = proposal_rate_density;
      const amrex::Real rate_eps = 1.0e-300;
      tr.state_acceptance_probability = std::max<amrex::Real>(
        0.0, std::min<amrex::Real>(
               1.0, tr.state_event_hazard / std::max<amrex::Real>(
                                            proposal_rate_density, rate_eps)));
      tr.acceptance_draw = sampler.sampleUnitUniform();
      tr.accepted_by_state_model =
        (tr.acceptance_draw <= tr.state_acceptance_probability);
      if (!tr.accepted_by_state_model) {
        tr.rejected = true;
        report.rejected_state_model_events += 1;
        report.rejected_events += 1;
        report.trace.push_back(tr);
        continue;
      }

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
    (event_rep_a.rejected_incompatible_events == 0);

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
    event_rep_a.rejected_state_model_events <= event_rep_a.rejected_events &&
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

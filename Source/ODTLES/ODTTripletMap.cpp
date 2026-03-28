#include "ODTTripletMap.H"

#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

namespace pelec::odtles
{

namespace
{

int
mappedSourceLocalIndex(int j, int n)
{
  AMREX_ALWAYS_ASSERT(n >= 3 && (n % 3) == 0);
  AMREX_ALWAYS_ASSERT(j >= 0 && j < n);
  const int m = n / 3;

  // Discrete triplet permutation over n=3m cells:
  // first third  : residues 0 (ascending)
  // middle third : residues 1 (descending)
  // final third  : residues 2 (ascending)
  // This is a pure reordering (measure-preserving permutation).
  if (j < m) {
    return 3 * j;
  }
  if (j < 2 * m) {
    const int i = j - m;
    return 3 * (m - 1 - i) + 1;
  }
  const int i = j - 2 * m;
  return 3 * i + 2;
}

amrex::Real
maxAbsStateChange(const ODTLineState& a, const ODTLineState& b)
{
  AMREX_ALWAYS_ASSERT(a.numCells() == b.numCells());
  amrex::Real max_abs = 0.0;
  for (int i = 0; i < a.numCells(); ++i) {
    const auto& ca = a.cell(i);
    const auto& cb = b.cell(i);
    max_abs = std::max(max_abs, std::abs(cb.rho - ca.rho));
    max_abs = std::max(max_abs, std::abs(cb.rhou - ca.rhou));
    max_abs = std::max(max_abs, std::abs(cb.rhov - ca.rhov));
    max_abs = std::max(max_abs, std::abs(cb.rhow - ca.rhow));
    max_abs = std::max(max_abs, std::abs(cb.rhoE - ca.rhoE));
  }
  return max_abs;
}

} // namespace

bool
ODTTripletMap::canApplyPermutation(
  const ODTLineState& line_state,
  const EddyInterval& interval)
{
  if (!line_state.initialized() || !line_state.valid()) {
    return false;
  }
  if (!interval.valid()) {
    return false;
  }
  if (interval.i_hi >= line_state.numCells()) {
    return false;
  }
  const int n = interval.size();
  // MVP triplet-map policy: interval size must be n=3m.
  return (n >= 3) && ((n % 3) == 0);
}

ODTTripletMap::ApplyReport
ODTTripletMap::applyPermutation(
  ODTLineState& line_state,
  const EddyInterval& interval)
{
  ApplyReport rep{};
  if (!interval.valid()) {
    return rep;
  }

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTTripletMap requires initialized, valid line state");
  AMREX_ALWAYS_ASSERT(interval.i_hi < line_state.numCells());

  const int n = interval.size();
  rep.interval_multiple_of_three = ((n % 3) == 0);
  if (n < 3 || !rep.interval_multiple_of_three) {
    return rep;
  }

  std::vector<ODTLineState::ConservativeCell> original(
    static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    original[static_cast<std::size_t>(j)] = line_state.cell(interval.i_lo + j);
  }

  for (int j = 0; j < n; ++j) {
    const int src = mappedSourceLocalIndex(j, n);
    line_state.setCell(
      interval.i_lo + j, original[static_cast<std::size_t>(src)]);
  }

  rep.applied = true;
  return rep;
}

ODTTripletMap::ApplyReport
ODTTripletMap::applyPermutationWithKernel(
  ODTLineState& line_state,
  const EddyInterval& interval,
  const ODTKernel::VelocityKernelSpec& kernel)
{
  auto rep = applyPermutation(line_state, interval);
  if (!rep.applied) {
    return rep;
  }

  rep.kernel_applied = ODTKernel::applyVelocityKernelZeroNet(
    line_state, interval.i_lo, interval.i_hi, kernel);
  return rep;
}

ODTTripletMap::ValidationReport
ODTTripletMap::runMVPValidationHook()
{
  ValidationReport report{};

  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(5, 0, 0)));
  amrex::RealBox rb(
    AMREX_D_DECL(0.0, 0.0, 0.0), AMREX_D_DECL(6.0, 1.0, 1.0));
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(3, 0, 0));
  ODTLineGeometry line_geom(geom, 0, owner, 0);

  ODTLineState base_state;
  base_state.initialize(line_geom);
  for (int i = 0; i < base_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0;
    c.rhou = 2.0;
    c.rhov = -1.0;
    c.rhow = 0.5;
    c.rhoE = 10.0;
    base_state.setCell(i, c);
  }
  base_state.setValid(true);

  constexpr amrex::Real tol = 1.0e-12;
  const EddyInterval interval{0, 2}; // n=3 satisfies MVP n=3m policy.

  // Uniform invariance under pure permutation.
  ODTLineState uniform_before = base_state;
  ODTLineState uniform_after = base_state;
  applyPermutation(uniform_after, interval);
  report.uniform_max_abs_change = maxAbsStateChange(uniform_before, uniform_after);
  report.uniform_invariant = (report.uniform_max_abs_change <= tol);

  // Nonuniform modification under permutation.
  ODTLineState nonuniform_before = base_state;
  for (int i = interval.i_lo; i <= interval.i_hi; ++i) {
    auto c = nonuniform_before.cell(i);
    c.rhou = static_cast<amrex::Real>(i + 1);
    nonuniform_before.setCell(i, c);
  }
  ODTLineState nonuniform_after = nonuniform_before;
  applyPermutation(nonuniform_after, interval);
  report.nonuniform_max_abs_change =
    maxAbsStateChange(nonuniform_before, nonuniform_after);
  report.nonuniform_modified = (report.nonuniform_max_abs_change > tol);

  // Locality check: outside eddy interval remains unchanged.
  report.outside_interval_max_abs_change = 0.0;
  for (int i = 0; i < nonuniform_before.numCells(); ++i) {
    if (i >= interval.i_lo && i <= interval.i_hi) {
      continue;
    }
    const auto& a = nonuniform_before.cell(i);
    const auto& b = nonuniform_after.cell(i);
    report.outside_interval_max_abs_change = std::max(
      report.outside_interval_max_abs_change, std::abs(b.rho - a.rho));
    report.outside_interval_max_abs_change = std::max(
      report.outside_interval_max_abs_change, std::abs(b.rhou - a.rhou));
    report.outside_interval_max_abs_change = std::max(
      report.outside_interval_max_abs_change, std::abs(b.rhov - a.rhov));
    report.outside_interval_max_abs_change = std::max(
      report.outside_interval_max_abs_change, std::abs(b.rhow - a.rhow));
    report.outside_interval_max_abs_change = std::max(
      report.outside_interval_max_abs_change, std::abs(b.rhoE - a.rhoE));
  }
  report.map_locality_preserved = (report.outside_interval_max_abs_change <= tol);

  // Kernel zero-net momentum check over mapped interval.
  ODTLineState kernel_state = nonuniform_after;
  amrex::Real sum_rhou_before = 0.0;
  amrex::Real sum_rhov_before = 0.0;
  amrex::Real sum_rhow_before = 0.0;
  for (int i = interval.i_lo; i <= interval.i_hi; ++i) {
    const auto& c = kernel_state.cell(i);
    sum_rhou_before += c.rhou;
    sum_rhov_before += c.rhov;
    sum_rhow_before += c.rhow;
  }

  ODTKernel::VelocityKernelSpec kernel{};
  kernel.enabled = true;
  kernel.amp_u = 0.3;
  kernel.amp_v = -0.2;
  kernel.amp_w = 0.1;
  kernel.shape = {0.0, 1.0, 0.0};

  applyPermutationWithKernel(kernel_state, interval, kernel);

  amrex::Real sum_rhou_after = 0.0;
  amrex::Real sum_rhov_after = 0.0;
  amrex::Real sum_rhow_after = 0.0;
  for (int i = interval.i_lo; i <= interval.i_hi; ++i) {
    const auto& c = kernel_state.cell(i);
    sum_rhou_after += c.rhou;
    sum_rhov_after += c.rhov;
    sum_rhow_after += c.rhow;
  }

  report.kernel_net_rhou_change = sum_rhou_after - sum_rhou_before;
  report.kernel_net_rhov_change = sum_rhov_after - sum_rhov_before;
  report.kernel_net_rhow_change = sum_rhow_after - sum_rhow_before;
  report.kernel_zero_net_momentum =
    (std::abs(report.kernel_net_rhou_change) <= tol) &&
    (std::abs(report.kernel_net_rhov_change) <= tol) &&
    (std::abs(report.kernel_net_rhow_change) <= tol);

  return report;
}

} // namespace pelec::odtles

#include "ODTMomentExtractor.H"

#include <algorithm>
#include <cmath>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

namespace pelec::odtles
{

namespace
{

amrex::Real
overlapLength(
  const ODTLineGeometry::Interval& a,
  const ODTLineGeometry::Interval& b)
{
  return std::max<amrex::Real>(
    0.0, std::min(a.hi, b.hi) - std::max(a.lo, b.lo));
}

amrex::Real
componentMomentum(const ODTLineState::ConservativeCell& c, int comp)
{
  AMREX_ALWAYS_ASSERT(comp >= 0 && comp < AMREX_SPACEDIM);
  if (comp == 0) {
    return c.rhou;
  }
  if (comp == 1) {
    return c.rhov;
  }
  return c.rhow;
}

amrex::Real
maxAbsDiff(const std::array<amrex::Real, AMREX_SPACEDIM>& a, const std::array<amrex::Real, AMREX_SPACEDIM>& b)
{
  amrex::Real m = 0.0;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    m = std::max(m, std::abs(a[static_cast<std::size_t>(i)] - b[static_cast<std::size_t>(i)]));
  }
  return m;
}

} // namespace

ODTMomentExtractor::DirectionalColumn
ODTMomentExtractor::extractDirectionalMomentumColumn(
  const ODTLineGeometry& geom,
  const ODTLineState& line_state)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTMomentExtractor requires defined line geometry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTMomentExtractor requires an initialized and valid line state");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.numCells() == geom.supportCellCount(),
    "ODTMomentExtractor requires state/geometry support-size consistency");

  DirectionalColumn out{};
  out.dir = geom.dir();
  AMREX_ALWAYS_ASSERT(out.dir >= 0 && out.dir < AMREX_SPACEDIM);

  const auto owner_interval = geom.ownerInterval();
  out.owner_interval_length = owner_interval.hi - owner_interval.lo;
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    out.owner_interval_length > 0.0,
    "ODTMomentExtractor requires positive owner-interval length");

  const auto& support = geom.supportCells();
  for (int k = 0; k < static_cast<int>(support.size()); ++k) {
    const auto& seg = support[static_cast<std::size_t>(k)];
    const amrex::Real overlap = overlapLength(owner_interval, seg.s_interval);
    if (overlap <= 0.0) {
      continue;
    }

    ++out.contributing_segments;
    out.overlap_measure += overlap;

    const amrex::Real wk = overlap / out.owner_interval_length;
    out.overlap_weight_sum += wk;

    const auto& c = line_state.cell(k);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      c.rho > 0.0,
      "ODTMomentExtractor requires strictly positive density in contributing "
      "segments");

    out.mean_rho += wk * c.rho;

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      const amrex::Real rho_ui = componentMomentum(c, i);
      out.mean_rho_ui[static_cast<std::size_t>(i)] += wk * rho_ui;
    }

    const amrex::Real rho_uj = componentMomentum(c, out.dir);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      const amrex::Real rho_ui = componentMomentum(c, i);
      const amrex::Real rho_ui_uj = (rho_ui * rho_uj) / c.rho;
      out.mean_rho_ui_uj[static_cast<std::size_t>(i)] += wk * rho_ui_uj;
    }
  }

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    out.contributing_segments > 0,
    "ODTMomentExtractor found no overlap between support and owner interval");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    out.mean_rho > 0.0,
    "ODTMomentExtractor computed non-positive central-interval mean density");

  const amrex::Real rho_uj_bar = out.mean_rho_ui[static_cast<std::size_t>(out.dir)];
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    out.tau_ij[static_cast<std::size_t>(i)] =
      out.mean_rho_ui_uj[static_cast<std::size_t>(i)] -
      (out.mean_rho_ui[static_cast<std::size_t>(i)] * rho_uj_bar) / out.mean_rho;
  }

  out.valid = true;
  return out;
}

ODTMomentExtractor::ValidationReport
ODTMomentExtractor::runMVPValidationHook()
{
  ValidationReport rep{};

  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(2, 0, 0)));
  amrex::RealBox rb(
    AMREX_D_DECL(0.0, 0.0, 0.0), AMREX_D_DECL(3.0, 1.0, 1.0));
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(1, 0, 0));

  ODTLineGeometry line_geom(geom, 0, owner, 0);
  ODTLineState line_state;
  line_state.initialize(line_geom);

  for (int i = 0; i < line_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0 + static_cast<amrex::Real>(i);
    c.rhou = 10.0 + static_cast<amrex::Real>(i);
    c.rhov = -2.0 + static_cast<amrex::Real>(0.5 * i);
    c.rhow = 0.25 * static_cast<amrex::Real>(i);
    c.rhoE = 20.0;
    line_state.setCell(i, c);
  }
  line_state.setValid(true);

  const auto col = extractDirectionalMomentumColumn(line_geom, line_state);
  AMREX_ALWAYS_ASSERT(col.valid);

  constexpr amrex::Real tol = 1.0e-12;
  rep.overlap_weight_sum_error = std::abs(col.overlap_weight_sum - 1.0);
  rep.overlap_weights_partition_unity = (rep.overlap_weight_sum_error <= tol);

  const auto owner_cell =
    line_state.cell(line_geom.ownerLocalOrdinal());
  std::array<amrex::Real, AMREX_SPACEDIM> owner_rho_ui{};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    owner_rho_ui[static_cast<std::size_t>(i)] = componentMomentum(owner_cell, i);
  }
  rep.owner_recovery_max_abs_error =
    std::max(std::abs(col.mean_rho - owner_cell.rho), maxAbsDiff(col.mean_rho_ui, owner_rho_ui));
  rep.owner_central_recovered = (rep.owner_recovery_max_abs_error <= tol);

  rep.tau_max_abs = 0.0;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    rep.tau_max_abs = std::max(rep.tau_max_abs, std::abs(col.tau_ij[i]));
  }
  rep.tau_zero_for_single_segment_central_interval = (rep.tau_max_abs <= tol);

  return rep;
}

} // namespace pelec::odtles

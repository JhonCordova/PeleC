#include "ODTReconcile.H"

#include <algorithm>
#include <array>
#include <cmath>

#include <AMReX.H>

namespace pelec::odtles
{

namespace
{

int
findOffsetIndex(const ODTLineGeometry& geom, int offset)
{
  const auto& cells = geom.supportCells();
  for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
    if (cells[static_cast<std::size_t>(i)].relative_offset == offset) {
      return i;
    }
  }
  return -1;
}

int
findNearestByOffsetSign(const ODTLineGeometry& geom, int sign)
{
  const auto& cells = geom.supportCells();
  int best = -1;
  int best_abs = 0;
  for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
    const int off = cells[static_cast<std::size_t>(i)].relative_offset;
    if ((sign < 0 && off < 0) || (sign > 0 && off > 0)) {
      const int aoff = std::abs(off);
      if (best < 0 || aoff < best_abs) {
        best = i;
        best_abs = aoff;
      }
    }
  }
  return best;
}

std::vector<int>
collectIndicesByOffsetSign(const ODTLineGeometry& geom, int sign)
{
  std::vector<int> idxs;
  const auto& cells = geom.supportCells();
  idxs.reserve(cells.size());
  for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
    const int off = cells[static_cast<std::size_t>(i)].relative_offset;
    if ((sign < 0 && off < 0) || (sign > 0 && off > 0)) {
      idxs.push_back(i);
    }
  }
  return idxs;
}

amrex::Real
computeLimitedSlope(
  amrex::Real ql,
  amrex::Real qc,
  amrex::Real qr,
  amrex::Real dx)
{
  const amrex::Real tol =
    1.0e-14 * (std::abs(ql) + std::abs(qc) + std::abs(qr) + 1.0);
  const bool uniform =
    (std::abs(qr - ql) <= tol && std::abs(qc - ql) <= tol);
  if (uniform) {
    return 0.0;
  }

  // Central slope from neighboring LES means.
  amrex::Real slope = (qr - ql) / (2.0 * dx);

  // If central slope vanishes but stencil is non-uniform, keep a directional
  // gradient from one side to avoid collapsing to a flat initialization.
  if (std::abs(slope) <= tol) {
    const amrex::Real dl = (qc - ql) / dx;
    const amrex::Real dr = (qr - qc) / dx;
    slope = (std::abs(dr) >= std::abs(dl)) ? dr : dl;
  }

  return slope;
}

void
setCellFromComponentArray(
  ODTLineState::ConservativeCell& cell,
  const std::vector<amrex::Real>& v)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(v.size()) == 5 + NUM_SPECIES,
    "ODTReconcile conservative component size mismatch");
  cell.rho = v[0];
  cell.rhou = v[1];
  cell.rhov = v[2];
  cell.rhow = v[3];
  cell.rhoE = v[4];
  if (cell.rhoY.size() != static_cast<std::size_t>(NUM_SPECIES)) {
    cell.rhoY.resize(static_cast<std::size_t>(NUM_SPECIES), 0.0);
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    cell.rhoY[static_cast<std::size_t>(n)] = v[5 + n];
  }
}

std::vector<amrex::Real>
componentArrayFromCell(const ODTLineState::ConservativeCell& c)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    c.rhoY.size() == static_cast<std::size_t>(NUM_SPECIES),
    "ODTReconcile expected rhoY container to match NUM_SPECIES");
  std::vector<amrex::Real> q(static_cast<std::size_t>(5 + NUM_SPECIES), 0.0);
  q[0] = c.rho;
  q[1] = c.rhou;
  q[2] = c.rhov;
  q[3] = c.rhow;
  q[4] = c.rhoE;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    q[static_cast<std::size_t>(5 + n)] = c.rhoY[static_cast<std::size_t>(n)];
  }
  return q;
}

amrex::Real
componentValue(const ODTLineState::ConservativeCell& c, int n)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    n >= 0 && n < 5 + NUM_SPECIES,
    "ODTReconcile conservative component index out of range");
  if (n == 0) {
    return c.rho;
  }
  if (n == 1) {
    return c.rhou;
  }
  if (n == 2) {
    return c.rhov;
  }
  if (n == 3) {
    return c.rhow;
  }
  if (n == 4) {
    return c.rhoE;
  }
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    c.rhoY.size() == static_cast<std::size_t>(NUM_SPECIES),
    "ODTReconcile expected rhoY container to match NUM_SPECIES");
  return c.rhoY[static_cast<std::size_t>(n - 5)];
}

amrex::Real
specificInternalEnergy(const std::vector<amrex::Real>& q)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(q.size()) == 5 + NUM_SPECIES,
    "ODTReconcile conservative component size mismatch");
  const amrex::Real rho = q[0];
  if (rho <= 0.0) {
    return -1.0;
  }
  const amrex::Real inv_rho = 1.0 / rho;
  const amrex::Real kinetic =
    0.5 * inv_rho * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  const amrex::Real rhoe = q[4] - kinetic;
  return rhoe * inv_rho;
}

bool
isAdmissible(
  const std::vector<amrex::Real>& q,
  const ODTReconcile::ReconcileControls& ctrl)
{
  if (q[0] <= ctrl.rho_floor) {
    return false;
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    if (q[static_cast<std::size_t>(5 + n)] < 0.0) {
      return false;
    }
  }
  const amrex::Real eint = specificInternalEnergy(q);
  return eint > ctrl.e_floor;
}

} // namespace

void
ODTReconcile::initializeLineStateFromLESSupportAverages(
  const ODTLineGeometry& geom,
  const std::vector<ConservativeCell>& local_support_cell_averages,
  ODTLineState& line_state)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTReconcile initialize requires defined geometry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(local_support_cell_averages.size()) == geom.supportCellCount(),
    "ODTReconcile initialize expects one LES average per support segment");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(geom.supportCells().size()) == geom.supportCellCount(),
    "ODTReconcile initialize requires geometry support metadata consistency");

  if (
    !line_state.initialized() ||
    line_state.numCells() != geom.supportCellCount()) {
    line_state.initialize(geom);
  }

  const int owner = geom.ownerLocalOrdinal();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    owner >= 0, "ODTReconcile initialize requires an owner segment in support");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    owner < static_cast<int>(local_support_cell_averages.size()),
    "ODTReconcile initialize owner index is out of range");

  // Input contract reminder: these are valid local LES support-cell averages
  // mapped to the line support ordering, not ghost-derived physical means.

  const int left_idx = findOffsetIndex(geom, -1);
  const int right_idx = findOffsetIndex(geom, 1);
  const int neg_idx = (left_idx >= 0) ? left_idx : findNearestByOffsetSign(geom, -1);
  const int pos_idx =
    (right_idx >= 0) ? right_idx : findNearestByOffsetSign(geom, 1);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    (neg_idx >= 0 || pos_idx >= 0),
    "ODTReconcile initialize requires at least one neighbor around owner");

  const amrex::Real dx = geom.deltaS();
  const auto owner_avg =
    local_support_cell_averages[static_cast<std::size_t>(owner)];
  constexpr int hydro_components = 5;
  const int ncomp = 5 + NUM_SPECIES;
  auto q0 = componentArrayFromCell(owner_avg);

  std::vector<amrex::Real> ql = q0;
  std::vector<amrex::Real> qr = q0;
  if (neg_idx >= 0) {
    ql = componentArrayFromCell(
      local_support_cell_averages[static_cast<std::size_t>(neg_idx)]);
  }
  if (pos_idx >= 0) {
    qr = componentArrayFromCell(
      local_support_cell_averages[static_cast<std::size_t>(pos_idx)]);
  }

  std::vector<amrex::Real> slopes(static_cast<std::size_t>(ncomp), 0.0);
  const auto neg_side = collectIndicesByOffsetSign(geom, -1);
  const auto pos_side = collectIndicesByOffsetSign(geom, 1);
  for (int n = 0; n < ncomp; ++n) {
    if (neg_idx >= 0 && pos_idx >= 0) {
      slopes[n] = computeLimitedSlope(ql[n], q0[n], qr[n], dx);
    } else {
      const std::vector<int>& side = (pos_idx >= 0) ? pos_side : neg_side;
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !side.empty(), "ODTReconcile initialize requires one-sided support cells");

      amrex::Real sum_s_dq = 0.0;
      amrex::Real sum_s2 = 0.0;
      amrex::Real slope_maxmag = 0.0;
      bool has_nonzero = false;
      bool has_nonuniform = false;
      const amrex::Real tol = 1.0e-14 * (std::abs(q0[n]) + 1.0);

      for (int i : side) {
        const auto& seg = geom.supportCells()[static_cast<std::size_t>(i)];
        const amrex::Real s = seg.s_center;
        if (s == 0.0) {
          continue;
        }
        const amrex::Real qi = componentValue(
          local_support_cell_averages[static_cast<std::size_t>(i)], n);
        const amrex::Real dqi = qi - q0[n];
        if (std::abs(dqi) > tol) {
          has_nonuniform = true;
        }
        sum_s_dq += s * dqi;
        sum_s2 += s * s;
        const amrex::Real mi = dqi / s;
        if (!has_nonzero || std::abs(mi) > std::abs(slope_maxmag)) {
          slope_maxmag = mi;
          has_nonzero = true;
        }
      }

      amrex::Real m = (sum_s2 > 0.0) ? (sum_s_dq / sum_s2) : 0.0;
      if (has_nonuniform && std::abs(m) <= tol && has_nonzero) {
        // Degeneracy guard: non-uniform one-sided support should not collapse
        // to a constant initialization.
        m = slope_maxmag;
      }
      slopes[n] = m;
    }
  }

  // Local overshoot control: clamp each reconstructed support average to the
  // envelope over the full local support stencil.
  std::vector<amrex::Real> qmin(static_cast<std::size_t>(ncomp), 0.0);
  std::vector<amrex::Real> qmax(static_cast<std::size_t>(ncomp), 0.0);
  for (int n = 0; n < ncomp; ++n) {
    qmin[n] = q0[n];
    qmax[n] = q0[n];
  }
  for (const auto& c : local_support_cell_averages) {
    for (int n = 0; n < ncomp; ++n) {
      const amrex::Real qi = componentValue(c, n);
      qmin[n] = std::min(qmin[n], qi);
      qmax[n] = std::max(qmax[n], qi);
    }
  }

  const auto& segs = geom.supportCells();
  for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
    const amrex::Real s_center = segs[static_cast<std::size_t>(i)].s_center;

    std::vector<amrex::Real> qrec(static_cast<std::size_t>(ncomp), 0.0);
    for (int n = 0; n < ncomp; ++n) {
      qrec[n] = q0[n] + slopes[n] * s_center;
      qrec[n] = std::clamp(qrec[n], qmin[n], qmax[n]);
      if (n >= hydro_components) {
        // Keep conservative species non-negative on the embedded line.
        qrec[n] = std::max<amrex::Real>(qrec[n], 0.0);
      }
    }

    ODTLineState::ConservativeCell cell{};
    setCellFromComponentArray(cell, qrec);
    line_state.setCell(i, cell);
  }

  // Owner interval is centered at s=0, so linear reconstruction preserves the
  // owner LES average exactly over the owner segment.
  line_state.setValid(true);
}

void
ODTReconcile::reconcileExistingLineStateToOwnerAverage(
  const ODTLineGeometry& geom,
  const ConservativeCell& owner_cell_average,
  ODTLineState& line_state)
{
  reconcileExistingLineStateToOwnerAverage(
    geom, owner_cell_average, line_state, defaultReconcileControls());
}

ODTReconcile::ReconcileControls
ODTReconcile::defaultReconcileControls()
{
  return ReconcileControls{};
}

void
ODTReconcile::reconcileExistingLineStateToOwnerAverage(
  const ODTLineGeometry& geom,
  const ConservativeCell& owner_cell_average,
  ODTLineState& line_state,
  const ReconcileControls& ctrl)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTReconcile reconcile requires defined geometry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTReconcile reconcile requires an existing valid line state");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.numCells() == geom.supportCellCount(),
    "ODTReconcile reconcile requires state/geometry support-size consistency");

  const int owner = geom.ownerLocalOrdinal();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    owner >= 0, "ODTReconcile reconcile requires owner segment in support");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    owner < line_state.numCells(),
    "ODTReconcile reconcile owner index out of range");

  constexpr int hydro_components = 5;
  const int ncomp = hydro_components + NUM_SPECIES;
  const auto q_target = componentArrayFromCell(owner_cell_average);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    isAdmissible(q_target, ctrl),
    "ODTReconcile reconcile target owner LES mean is not admissible");

  const auto q_owner_old = componentArrayFromCell(line_state.cell(owner));

  // Residual relative to old owner mean at each support cell.
  std::vector<std::vector<amrex::Real>> residuals(
    static_cast<std::size_t>(line_state.numCells()));
  for (int i = 0; i < line_state.numCells(); ++i) {
    const auto qi = componentArrayFromCell(line_state.cell(i));
    auto& ri = residuals[static_cast<std::size_t>(i)];
    ri.assign(static_cast<std::size_t>(ncomp), 0.0);
    for (int n = 0; n < ncomp; ++n) {
      ri[n] = qi[n] - q_owner_old[n];
    }
  }

  auto make_trial = [&](int i, amrex::Real alpha) -> std::vector<amrex::Real> {
    std::vector<amrex::Real> q(static_cast<std::size_t>(ncomp), 0.0);
    const auto& r = residuals[static_cast<std::size_t>(i)];
    for (int n = 0; n < ncomp; ++n) {
      q[n] = q_target[n] + alpha * r[n];
    }
    return q;
  };

  // Uniform residual scaling for admissibility: alpha starts at 1 and is
  // reduced only if needed.
  amrex::Real alpha = 1.0;
  bool admissible = false;
  for (int iter = 0; iter <= ctrl.max_alpha_reductions; ++iter) {
    admissible = true;
    for (int i = 0; i < line_state.numCells(); ++i) {
      if (!isAdmissible(make_trial(i, alpha), ctrl)) {
        admissible = false;
        break;
      }
    }
    if (admissible) {
      break;
    }
    alpha *= 0.5;
    if (alpha < ctrl.alpha_min) {
      alpha = 0.0;
      break;
    }
  }

  if (alpha == 0.0) {
    admissible = true;
    for (int i = 0; i < line_state.numCells(); ++i) {
      if (!isAdmissible(q_target, ctrl)) {
        admissible = false;
        break;
      }
    }
  }
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    admissible,
    "ODTReconcile reconcile failed admissibility even after residual scaling");

  for (int i = 0; i < line_state.numCells(); ++i) {
    const auto q = make_trial(i, alpha);
    ODTLineState::ConservativeCell c{};
    setCellFromComponentArray(c, q);
    line_state.setCell(i, c);
  }

  // Explicit mean treatment: enforce owner exactly after positivity control.
  ODTLineState::ConservativeCell owner_cell{};
  setCellFromComponentArray(owner_cell, q_target);
  line_state.setCell(owner, owner_cell);
}

} // namespace pelec::odtles

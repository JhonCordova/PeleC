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
  const std::array<amrex::Real, 5>& v)
{
  cell.rho = v[0];
  cell.rhou = v[1];
  cell.rhov = v[2];
  cell.rhow = v[3];
  cell.rhoE = v[4];
}

std::array<amrex::Real, 5>
componentArrayFromCell(const ODTLineState::ConservativeCell& c)
{
  return {c.rho, c.rhou, c.rhov, c.rhow, c.rhoE};
}

amrex::Real
specificInternalEnergy(const std::array<amrex::Real, 5>& q)
{
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
  const std::array<amrex::Real, 5>& q,
  const ODTReconcile::ReconcileControls& ctrl)
{
  if (q[0] <= ctrl.rho_floor) {
    return false;
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
  std::array<amrex::Real, 5> q0{
    owner_avg.rho, owner_avg.rhou, owner_avg.rhov, owner_avg.rhow, owner_avg.rhoE};

  std::array<amrex::Real, 5> ql = q0;
  std::array<amrex::Real, 5> qr = q0;
  if (neg_idx >= 0) {
    const auto l = local_support_cell_averages[static_cast<std::size_t>(neg_idx)];
    ql = {l.rho, l.rhou, l.rhov, l.rhow, l.rhoE};
  }
  if (pos_idx >= 0) {
    const auto r = local_support_cell_averages[static_cast<std::size_t>(pos_idx)];
    qr = {r.rho, r.rhou, r.rhov, r.rhow, r.rhoE};
  }

  std::array<amrex::Real, 5> slopes{};
  const auto neg_side = collectIndicesByOffsetSign(geom, -1);
  const auto pos_side = collectIndicesByOffsetSign(geom, 1);
  for (int n = 0; n < 5; ++n) {
    auto compVal = [&](const ConservativeCell& c) -> amrex::Real {
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
      return c.rhoE;
    };

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
        const amrex::Real qi =
          compVal(local_support_cell_averages[static_cast<std::size_t>(i)]);
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
  std::array<amrex::Real, 5> qmin{};
  std::array<amrex::Real, 5> qmax{};
  for (int n = 0; n < 5; ++n) {
    qmin[n] = q0[n];
    qmax[n] = q0[n];
  }
  for (const auto& c : local_support_cell_averages) {
    qmin[0] = std::min(qmin[0], c.rho);
    qmin[1] = std::min(qmin[1], c.rhou);
    qmin[2] = std::min(qmin[2], c.rhov);
    qmin[3] = std::min(qmin[3], c.rhow);
    qmin[4] = std::min(qmin[4], c.rhoE);

    qmax[0] = std::max(qmax[0], c.rho);
    qmax[1] = std::max(qmax[1], c.rhou);
    qmax[2] = std::max(qmax[2], c.rhov);
    qmax[3] = std::max(qmax[3], c.rhow);
    qmax[4] = std::max(qmax[4], c.rhoE);
  }

  const auto& segs = geom.supportCells();
  for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
    const amrex::Real s_center = segs[static_cast<std::size_t>(i)].s_center;

    std::array<amrex::Real, 5> qrec{};
    for (int n = 0; n < 5; ++n) {
      qrec[n] = q0[n] + slopes[n] * s_center;
      qrec[n] = std::clamp(qrec[n], qmin[n], qmax[n]);
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

  const auto q_target = componentArrayFromCell(owner_cell_average);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    isAdmissible(q_target, ctrl),
    "ODTReconcile reconcile target owner LES mean is not admissible");

  const auto q_owner_old = componentArrayFromCell(line_state.cell(owner));

  // Residual relative to old owner mean at each support cell.
  std::vector<std::array<amrex::Real, 5>> residuals(
    static_cast<std::size_t>(line_state.numCells()));
  for (int i = 0; i < line_state.numCells(); ++i) {
    const auto qi = componentArrayFromCell(line_state.cell(i));
    auto& ri = residuals[static_cast<std::size_t>(i)];
    for (int n = 0; n < 5; ++n) {
      ri[n] = qi[n] - q_owner_old[n];
    }
  }

  auto make_trial =
    [&](int i, amrex::Real alpha) -> std::array<amrex::Real, 5> {
    std::array<amrex::Real, 5> q{};
    const auto& r = residuals[static_cast<std::size_t>(i)];
    for (int n = 0; n < 5; ++n) {
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

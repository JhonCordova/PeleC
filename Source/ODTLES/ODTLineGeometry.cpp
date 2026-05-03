#include "ODTLineGeometry.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pelec::odtles
{

int
ODTLineGeometry::wrapToDomain(int idx, int ilo, int ihi)
{
  const int n = ihi - ilo + 1;
  if (n <= 0) {
    return idx;
  }
  int wrapped = (idx - ilo) % n;
  if (wrapped < 0) {
    wrapped += n;
  }
  return wrapped + ilo;
}

ODTLineGeometry::ODTLineGeometry(
  const amrex::Geometry& geom,
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  int subsegments_per_host_cell)
{
  define(geom, level, owner_cell, dir, subsegments_per_host_cell);
}

void
ODTLineGeometry::define(
  const amrex::Geometry& geom,
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  int subsegments_per_host_cell)
{
  if (dir < 0 || dir >= AMREX_SPACEDIM) {
    throw std::invalid_argument("ODTLineGeometry::define invalid line direction");
  }
  if (subsegments_per_host_cell < 3) {
    throw std::invalid_argument(
      "ODTLineGeometry::define requires subsegments_per_host_cell >= 3");
  }
  if ((subsegments_per_host_cell % 2) == 0) {
    throw std::invalid_argument(
      "ODTLineGeometry::define requires odd subsegments_per_host_cell");
  }

  m_defined = false;
  m_support_cells.clear();
  m_owner_local_ordinal = -1;
  m_subsegments_per_host_cell = subsegments_per_host_cell;

  m_key.level = level;
  m_key.owner_cell = owner_cell;
  m_key.dir = dir;

  const auto dx = geom.CellSizeArray();
  const auto prob_lo = geom.ProbLoArray();
  const auto domain = geom.Domain();
  const int ilo = domain.smallEnd(dir);
  const int ihi = domain.bigEnd(dir);
  const int owner_idx = owner_cell[dir];

  m_delta_s = dx[dir];
  m_owner_center_phys =
    prob_lo[dir] + (static_cast<amrex::Real>(owner_idx) + 0.5) * m_delta_s;
  m_owner_interval = {-0.5 * m_delta_s, 0.5 * m_delta_s};

  m_boundary.periodic_dir = geom.isPeriodic(dir);

  std::vector<int> host_offsets;
  std::vector<int> host_indices;
  std::vector<bool> wrapped_flags;

  if (m_boundary.periodic_dir) {
    host_offsets = {-1, 0, 1};
    host_indices.resize(host_offsets.size());
    wrapped_flags.resize(host_offsets.size(), false);
    for (int n = 0; n < static_cast<int>(host_offsets.size()); ++n) {
      const int raw = owner_idx + host_offsets[n];
      host_indices[n] = wrapToDomain(raw, ilo, ihi);
      wrapped_flags[n] = (raw < ilo || raw > ihi);
    }
  } else {
    int lo = owner_idx - 1;
    int hi = owner_idx + 1;
    if (lo < ilo) {
      hi += (ilo - lo);
      lo = ilo;
    }
    if (hi > ihi) {
      lo -= (hi - ihi);
      hi = ihi;
    }
    lo = std::max(lo, ilo);
    hi = std::min(hi, ihi);

    m_boundary.clipped_low = (owner_idx - 1 < ilo);
    m_boundary.clipped_high = (owner_idx + 1 > ihi);
    m_boundary.touches_low_boundary = (lo == ilo);
    m_boundary.touches_high_boundary = (hi == ihi);

    const int n_cells = hi - lo + 1;
    host_offsets.reserve(n_cells);
    host_indices.reserve(n_cells);
    wrapped_flags.assign(n_cells, false);
    for (int idx = lo; idx <= hi; ++idx) {
      host_offsets.push_back(idx - owner_idx);
      host_indices.push_back(idx);
    }
  }

  const int n_sub = m_subsegments_per_host_cell;
  const amrex::Real sub_width = m_delta_s / static_cast<amrex::Real>(n_sub);
  int local_ordinal = 0;
  for (int n = 0; n < static_cast<int>(host_offsets.size()); ++n) {
    const int host_offset = host_offsets[n];
    const amrex::Real host_center =
      static_cast<amrex::Real>(host_offset) * m_delta_s;
    for (int sub = 0; sub < n_sub; ++sub) {
      const amrex::Real lo_local =
        -0.5 * m_delta_s +
        static_cast<amrex::Real>(sub) * sub_width;
      const amrex::Real hi_local = lo_local + sub_width;

      CellSegment seg;
      seg.local_ordinal = local_ordinal++;
      seg.relative_offset = host_offset;
      seg.host_cell_index = host_indices[n];
      seg.wrapped = wrapped_flags[n];
      seg.s_interval = {host_center + lo_local, host_center + hi_local};
      seg.s_center = 0.5 * (seg.s_interval.lo + seg.s_interval.hi);
      m_support_cells.push_back(seg);

      const amrex::Real contains_origin_tol = 1.0e-14 * (m_delta_s + 1.0);
      if (
        seg.s_interval.lo <= 0.0 + contains_origin_tol &&
        seg.s_interval.hi >= 0.0 - contains_origin_tol) {
        m_owner_local_ordinal = seg.local_ordinal;
      }
    }
  }

  if (m_owner_local_ordinal < 0 && !m_support_cells.empty()) {
    int best = 0;
    amrex::Real best_abs = std::abs(m_support_cells[0].s_center);
    for (int i = 1; i < static_cast<int>(m_support_cells.size()); ++i) {
      const amrex::Real a = std::abs(m_support_cells[static_cast<std::size_t>(i)].s_center);
      if (a < best_abs) {
        best = i;
        best_abs = a;
      }
    }
    m_owner_local_ordinal = best;
  }

  if (!m_support_cells.empty()) {
    m_support_interval.lo = m_support_cells.front().s_interval.lo;
    m_support_interval.hi = m_support_cells.front().s_interval.hi;
    for (const auto& seg : m_support_cells) {
      m_support_interval.lo = std::min(m_support_interval.lo, seg.s_interval.lo);
      m_support_interval.hi = std::max(m_support_interval.hi, seg.s_interval.hi);
    }
  } else {
    m_support_interval = {0.0, 0.0};
  }

  m_support_length = m_support_interval.hi - m_support_interval.lo;
  m_boundary.one_sided_low =
    (!m_boundary.periodic_dir && !m_support_cells.empty() &&
     m_support_cells.front().relative_offset >= 0);
  m_boundary.one_sided_high =
    (!m_boundary.periodic_dir && !m_support_cells.empty() &&
     m_support_cells.back().relative_offset <= 0);

  m_defined = true;
}

} // namespace pelec::odtles

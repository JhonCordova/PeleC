#include "ODTLineGeometry.H"

#include <algorithm>
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
  int dir)
{
  define(geom, level, owner_cell, dir);
}

void
ODTLineGeometry::define(
  const amrex::Geometry& geom,
  int level,
  const amrex::IntVect& owner_cell,
  int dir)
{
  if (dir < 0 || dir >= AMREX_SPACEDIM) {
    throw std::invalid_argument("ODTLineGeometry::define invalid line direction");
  }

  m_defined = false;
  m_support_cells.clear();
  m_owner_local_ordinal = -1;

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

  std::vector<int> offsets;
  std::vector<int> host_indices;
  std::vector<bool> wrapped_flags;

  if (m_boundary.periodic_dir) {
    offsets = {-1, 0, 1};
    host_indices.resize(offsets.size());
    wrapped_flags.resize(offsets.size(), false);
    for (int n = 0; n < static_cast<int>(offsets.size()); ++n) {
      const int raw = owner_idx + offsets[n];
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
    offsets.reserve(n_cells);
    host_indices.reserve(n_cells);
    wrapped_flags.assign(n_cells, false);
    for (int idx = lo; idx <= hi; ++idx) {
      offsets.push_back(idx - owner_idx);
      host_indices.push_back(idx);
    }
  }

  for (int n = 0; n < static_cast<int>(offsets.size()); ++n) {
    CellSegment seg;
    seg.local_ordinal = n;
    seg.relative_offset = offsets[n];
    seg.host_cell_index = host_indices[n];
    seg.wrapped = wrapped_flags[n];
    seg.s_center = static_cast<amrex::Real>(seg.relative_offset) * m_delta_s;
    seg.s_interval = {
      seg.s_center - 0.5 * m_delta_s, seg.s_center + 0.5 * m_delta_s};
    m_support_cells.push_back(seg);
    if (seg.relative_offset == 0) {
      m_owner_local_ordinal = n;
    }
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

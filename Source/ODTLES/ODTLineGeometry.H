#ifndef PELEC_ODTLES_ODTLINEGEOMETRY_H
#define PELEC_ODTLES_ODTLINEGEOMETRY_H

#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

namespace pelec::odtles
{

class ODTLineGeometry
{
public:
  static constexpr int MVPNumSupportCells = 3;

  struct OwnershipKey
  {
    int level = -1;
    amrex::IntVect owner_cell{AMREX_D_DECL(0, 0, 0)};
    int dir = 0;
  };

  struct Interval
  {
    amrex::Real lo = 0.0;
    amrex::Real hi = 0.0;
  };

  struct CellSegment
  {
    int local_ordinal = -1;
    int relative_offset = 0;
    int host_cell_index = 0;
    bool wrapped = false;
    Interval s_interval;
    amrex::Real s_center = 0.0;
  };

  struct BoundaryMetadata
  {
    bool periodic_dir = false;
    bool touches_low_boundary = false;
    bool touches_high_boundary = false;
    bool clipped_low = false;
    bool clipped_high = false;
    bool one_sided_low = false;
    bool one_sided_high = false;
  };

  ODTLineGeometry() = default;
  ODTLineGeometry(
    const amrex::Geometry& geom,
    int level,
    const amrex::IntVect& owner_cell,
    int dir);

  void define(
    const amrex::Geometry& geom,
    int level,
    const amrex::IntVect& owner_cell,
    int dir);

  bool isDefined() const noexcept { return m_defined; }

  const OwnershipKey& key() const noexcept { return m_key; }
  int level() const noexcept { return m_key.level; }
  const amrex::IntVect& ownerCell() const noexcept { return m_key.owner_cell; }
  int dir() const noexcept { return m_key.dir; }

  amrex::Real deltaS() const noexcept { return m_delta_s; }
  amrex::Real ownerCenterPhysical() const noexcept { return m_owner_center_phys; }

  int supportCellCount() const noexcept
  {
    return static_cast<int>(m_support_cells.size());
  }
  amrex::Real supportLength() const noexcept { return m_support_length; }
  const Interval& supportInterval() const noexcept { return m_support_interval; }
  const Interval& ownerInterval() const noexcept { return m_owner_interval; }
  int ownerLocalOrdinal() const noexcept { return m_owner_local_ordinal; }

  const std::vector<CellSegment>& supportCells() const noexcept
  {
    return m_support_cells;
  }

  const BoundaryMetadata& boundary() const noexcept { return m_boundary; }

private:
  static int wrapToDomain(int idx, int ilo, int ihi);

  bool m_defined = false;
  OwnershipKey m_key{};
  amrex::Real m_delta_s = 0.0;
  amrex::Real m_owner_center_phys = 0.0;
  amrex::Real m_support_length = 0.0;
  Interval m_support_interval{};
  Interval m_owner_interval{};
  int m_owner_local_ordinal = -1;
  std::vector<CellSegment> m_support_cells;
  BoundaryMetadata m_boundary{};
};

} // namespace pelec::odtles

#endif

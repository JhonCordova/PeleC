#include "ODTManager.H"

#include <AMReX.H>

namespace pelec::odtles
{

ODTManager::ODTManager(const ODTParams& params) : m_params(params) {}

void
ODTManager::setParams(const ODTParams& params)
{
  m_params = params;
}

void
ODTManager::initializeLevel(
  int level,
  const amrex::BoxArray& grids,
  const amrex::DistributionMapping& dmap)
{
  static_cast<void>(grids);
  static_cast<void>(dmap);
  m_level = level;
  m_initialized = true;
}

void
ODTManager::clear()
{
  m_line_geometries.clear();
  m_line_states.clear();
  m_initialized = false;
  m_level = -1;
}

ODTLineGeometry&
ODTManager::getOrCreateLineGeometry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom)
{
  const LineKey key{level, owner_cell, dir};
  auto it = m_line_geometries.find(key);
  if (it == m_line_geometries.end()) {
    it =
      m_line_geometries.emplace(key, ODTLineGeometry(geom, level, owner_cell, dir))
        .first;
  } else {
    it->second.define(geom, level, owner_cell, dir);
  }
  return it->second;
}

ODTLineGeometry*
ODTManager::findLineGeometry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir)
{
  const LineKey key{level, owner_cell, dir};
  const auto it = m_line_geometries.find(key);
  return (it == m_line_geometries.end()) ? nullptr : &it->second;
}

const ODTLineGeometry*
ODTManager::findLineGeometry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir) const
{
  const LineKey key{level, owner_cell, dir};
  const auto it = m_line_geometries.find(key);
  return (it == m_line_geometries.end()) ? nullptr : &it->second;
}

ODTLineState&
ODTManager::getOrCreateLineState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const ODTLineGeometry& geom)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTManager::getOrCreateLineState requires defined geometry");
  AMREX_ALWAYS_ASSERT(geom.level() == level);
  AMREX_ALWAYS_ASSERT(geom.ownerCell() == owner_cell);
  AMREX_ALWAYS_ASSERT(geom.dir() == dir);

  const LineKey key{level, owner_cell, dir};
  auto it_inserted_pair = m_line_states.emplace(key, ODTLineState{});
  auto it = it_inserted_pair.first;
  if (
    it_inserted_pair.second ||
    !it->second.initialized() ||
    it->second.numCells() != geom.supportCellCount()) {
    it->second.initialize(geom);
  }
  return it->second;
}

ODTLineState*
ODTManager::findLineState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir)
{
  const LineKey key{level, owner_cell, dir};
  const auto it = m_line_states.find(key);
  return (it == m_line_states.end()) ? nullptr : &it->second;
}

const ODTLineState*
ODTManager::findLineState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir) const
{
  const LineKey key{level, owner_cell, dir};
  const auto it = m_line_states.find(key);
  return (it == m_line_states.end()) ? nullptr : &it->second;
}

bool
ODTManager::LineKey::operator<(const LineKey& other) const noexcept
{
  if (level != other.level) {
    return level < other.level;
  }
  if (owner_cell[0] != other.owner_cell[0]) {
    return owner_cell[0] < other.owner_cell[0];
  }
#if AMREX_SPACEDIM >= 2
  if (owner_cell[1] != other.owner_cell[1]) {
    return owner_cell[1] < other.owner_cell[1];
  }
#endif
#if AMREX_SPACEDIM == 3
  if (owner_cell[2] != other.owner_cell[2]) {
    return owner_cell[2] < other.owner_cell[2];
  }
#endif
  return dir < other.dir;
}

} // namespace pelec::odtles

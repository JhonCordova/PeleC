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
  if (m_initialized) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      level == m_level,
      "ODTManager::initializeLevel called with a different level than the "
      "already initialized level");
  }
  m_level = level;
  m_initialized = true;
}

void
ODTManager::clear()
{
  m_line_entries.clear();
  m_initialized = false;
  m_level = -1;
}

bool
ODTManager::hasLineEntry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir) const
{
  return findLineEntry(level, owner_cell, dir) != nullptr;
}

ODTManager::LineEntry*
ODTManager::findLineEntry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir)
{
  assertLevelScope(level);
  const LineKey key{level, owner_cell, dir};
  const auto it = m_line_entries.find(key);
  return (it == m_line_entries.end()) ? nullptr : &it->second;
}

const ODTManager::LineEntry*
ODTManager::findLineEntry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir) const
{
  assertLevelScope(level);
  const LineKey key{level, owner_cell, dir};
  const auto it = m_line_entries.find(key);
  return (it == m_line_entries.end()) ? nullptr : &it->second;
}

ODTManager::LineEntry&
ODTManager::getOrCreateLineEntry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom)
{
  assertLevelScope(level);
  const LineKey key{level, owner_cell, dir};
  auto it_inserted_pair = m_line_entries.emplace(key, LineEntry{});
  auto it = it_inserted_pair.first;
  if (it_inserted_pair.second) {
    it->second.geometry.define(geom, level, owner_cell, dir);
    it->second.state.initialize(it->second.geometry);
  } else {
    AMREX_ALWAYS_ASSERT(it->second.geometry.isDefined());
    AMREX_ALWAYS_ASSERT(it->second.geometry.level() == level);
    AMREX_ALWAYS_ASSERT(it->second.geometry.ownerCell() == owner_cell);
    AMREX_ALWAYS_ASSERT(it->second.geometry.dir() == dir);
    AMREX_ALWAYS_ASSERT(it->second.state.initialized());
    AMREX_ALWAYS_ASSERT(
      it->second.state.numCells() == it->second.geometry.supportCellCount());
  }
  return it->second;
}

ODTLineGeometry&
ODTManager::getOrCreateLineGeometry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom)
{
  return getOrCreateLineEntry(level, owner_cell, dir, geom).geometry;
}

ODTLineGeometry*
ODTManager::findLineGeometry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir)
{
  auto* entry = findLineEntry(level, owner_cell, dir);
  return entry == nullptr ? nullptr : &entry->geometry;
}

const ODTLineGeometry*
ODTManager::findLineGeometry(
  int level,
  const amrex::IntVect& owner_cell,
  int dir) const
{
  const auto* entry = findLineEntry(level, owner_cell, dir);
  return entry == nullptr ? nullptr : &entry->geometry;
}

ODTLineState&
ODTManager::getOrCreateLineState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom)
{
  return getOrCreateLineEntry(level, owner_cell, dir, geom).state;
}

ODTLineState*
ODTManager::findLineState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir)
{
  auto* entry = findLineEntry(level, owner_cell, dir);
  return entry == nullptr ? nullptr : &entry->state;
}

const ODTLineState*
ODTManager::findLineState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir) const
{
  const auto* entry = findLineEntry(level, owner_cell, dir);
  return entry == nullptr ? nullptr : &entry->state;
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

void
ODTManager::assertLevelScope(int level) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    m_initialized,
    "ODTManager persistent access requires initializeLevel() first");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    level == m_level,
    "ODTManager is level-local; mixed-level access is not allowed");
}

} // namespace pelec::odtles

#include "ODTManager.H"

#include <AMReX.H>

#include "IndexDefines.H"
#include "ODTReconcile.H"

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

ODTLineState&
ODTManager::initializeLineStateFromLESState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom,
  const amrex::MultiFab& state)
{
  auto& entry = getOrCreateLineEntry(level, owner_cell, dir, geom);
  const auto local_support_cell_averages =
    assembleLocalSupportCellAverages(entry.geometry, state);
  return initializeLineStateFromLESSupportAverages(
    level, owner_cell, dir, geom, local_support_cell_averages);
}

std::vector<ODTLineState::ConservativeCell>
ODTManager::assembleLocalSupportCellAverages(
  const ODTLineGeometry& line_geom,
  const amrex::MultiFab& state) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_geom.isDefined(),
    "ODTManager::assembleLocalSupportCellAverages requires defined geometry");

  const auto& support = line_geom.supportCells();
  std::vector<ODTLineState::ConservativeCell> local_support_cell_averages(
    support.size());
  const amrex::IntVect owner = line_geom.ownerCell();
  const int dir = line_geom.dir();

  for (int i = 0; i < static_cast<int>(support.size()); ++i) {
    const auto& seg = support[static_cast<std::size_t>(i)];
    amrex::IntVect iv = owner;
    iv[dir] = seg.host_cell_index;

    bool found = false;
    for (amrex::MFIter mfi(state, false); mfi.isValid(); ++mfi) {
      const amrex::Box vbx = mfi.validbox();
      if (!vbx.contains(iv)) {
        continue;
      }

      const auto arr = state.const_array(mfi);
      auto& c = local_support_cell_averages[static_cast<std::size_t>(i)];
      c.rho = arr(iv, URHO);
      c.rhou = arr(iv, UMX);
      c.rhov = arr(iv, UMY);
      c.rhow = arr(iv, UMZ);
      c.rhoE = arr(iv, UEDEN);
      found = true;
      break;
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      found,
      "ODTManager support-cell average assembly requires valid LES cell "
      "ownership (ghost values are not accepted as physical means)");
  }

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(local_support_cell_averages.size()) ==
      line_geom.supportCellCount(),
    "ODTManager support-cell average assembly size mismatch");

  return local_support_cell_averages;
}

ODTLineState&
ODTManager::initializeLineStateFromLESSupportAverages(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom,
  const std::vector<ODTLineState::ConservativeCell>&
    local_support_cell_averages)
{
  auto& entry = getOrCreateLineEntry(level, owner_cell, dir, geom);
  if (!entry.state.valid()) {
    ODTReconcile::initializeLineStateFromLESSupportAverages(
      entry.geometry, local_support_cell_averages, entry.state);
  }
  return entry.state;
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

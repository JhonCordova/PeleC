#include "ODTManager.H"

#include <algorithm>
#include <cmath>

#include <AMReX.H>

#include "IndexDefines.H"
#include "ODTReconcile.H"

namespace pelec::odtles
{

ODTManager::ODTManager(const ODTParams& params) : m_params(params) {}

int
ODTManager::SupportData::count(SupportProvenance p) const noexcept
{
  int n = 0;
  for (const auto& s : ordered_samples) {
    if (s.provenance == p) {
      ++n;
    }
  }
  return n;
}

bool
ODTManager::SupportData::allSameLevelAccepted() const noexcept
{
  for (const auto& s : ordered_samples) {
    if (
      s.provenance != SupportProvenance::SameLevelValidCell &&
      s.provenance != SupportProvenance::SameLevelFilledGhost) {
      return false;
    }
    if (!s.has_value) {
      return false;
    }
  }
  return true;
}

std::vector<ODTManager::ConservativeCell>
ODTManager::SupportData::conservativeAveragesOrdered() const
{
  std::vector<ConservativeCell> vals;
  vals.reserve(ordered_samples.size());
  for (const auto& s : ordered_samples) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      s.has_value,
      "SupportData cannot provide conservative averages for missing values");
    vals.push_back(s.value);
  }
  return vals;
}

ODTManager::ConservativeCell
ODTManager::SupportData::ownerAverage(const ODTLineGeometry& geom) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(ordered_samples.size()) == geom.supportCellCount(),
    "SupportData owner average requires sample/geometry support-size consistency");

  const auto owner_interval = geom.ownerInterval();
  const amrex::Real owner_len = owner_interval.hi - owner_interval.lo;
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    owner_len > 0.0, "SupportData owner interval length must be positive");

  ConservativeCell out{};
  amrex::Real sum_w = 0.0;
  const auto& support = geom.supportCells();
  for (int i = 0; i < static_cast<int>(ordered_samples.size()); ++i) {
    const auto& s = ordered_samples[static_cast<std::size_t>(i)];
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      s.has_value, "SupportData owner average requires populated support values");

    const auto& seg = support[static_cast<std::size_t>(i)];
    const amrex::Real overlap = std::max<amrex::Real>(
      0.0,
      std::min(owner_interval.hi, seg.s_interval.hi) -
        std::max(owner_interval.lo, seg.s_interval.lo));
    if (overlap <= 0.0) {
      continue;
    }

    const amrex::Real w = overlap / owner_len;
    sum_w += w;

    out.rho += w * s.value.rho;
    out.rhou += w * s.value.rhou;
    out.rhov += w * s.value.rhov;
    out.rhow += w * s.value.rhow;
    out.rhoE += w * s.value.rhoE;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      s.value.rhoY.size() == static_cast<std::size_t>(NUM_SPECIES),
      "SupportData owner average expected rhoY size to match NUM_SPECIES");
    for (int n = 0; n < NUM_SPECIES; ++n) {
      out.rhoY[static_cast<std::size_t>(n)] +=
        w * s.value.rhoY[static_cast<std::size_t>(n)];
    }
  }

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    sum_w > 0.0, "SupportData owner average found no owner-interval overlap");

  const amrex::Real tol = 1.0e-12;
  if (std::abs(sum_w - 1.0) > tol) {
    const amrex::Real inv_w = 1.0 / sum_w;
    out.rho *= inv_w;
    out.rhou *= inv_w;
    out.rhov *= inv_w;
    out.rhow *= inv_w;
    out.rhoE *= inv_w;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      out.rhoY[static_cast<std::size_t>(n)] *= inv_w;
    }
  }

  return out;
}

void
ODTManager::setParams(const ODTParams& params)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    params.subsegments_per_host_cell >= 3,
    "ODTManager::setParams requires subsegments_per_host_cell >= 3");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    (params.subsegments_per_host_cell % 2) == 1,
    "ODTManager::setParams requires odd subsegments_per_host_cell");
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
    it->second.geometry.define(
      geom, level, owner_cell, dir, m_params.subsegments_per_host_cell);
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

bool
ODTManager::findCellAverageInMF(
  const amrex::IntVect& iv,
  const amrex::MultiFab& state,
  bool use_validbox_only,
  ConservativeCell& out_cell)
{
  for (amrex::MFIter mfi(state, false); mfi.isValid(); ++mfi) {
    const amrex::Box bx = use_validbox_only ? mfi.validbox() : mfi.fabbox();
    if (!bx.contains(iv)) {
      continue;
    }

    const auto arr = state.const_array(mfi);
    out_cell.rho = arr(iv, URHO);
    out_cell.rhou = arr(iv, UMX);
    out_cell.rhov = arr(iv, UMY);
    out_cell.rhow = arr(iv, UMZ);
    out_cell.rhoE = arr(iv, UEDEN);
    if (out_cell.rhoY.size() != static_cast<std::size_t>(NUM_SPECIES)) {
      out_cell.rhoY.resize(static_cast<std::size_t>(NUM_SPECIES), 0.0);
    }
    for (int n = 0; n < NUM_SPECIES; ++n) {
      out_cell.rhoY[static_cast<std::size_t>(n)] = arr(iv, UFS + n);
    }
    return true;
  }
  return false;
}

ODTManager::SupportData
ODTManager::collectSupportData(
  const ODTLineGeometry& line_geom,
  const amrex::Geometry& geom,
  const amrex::MultiFab& state_valid,
  const amrex::MultiFab& state_same_level,
  const amrex::MultiFab& state_host_filled) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_geom.isDefined(), "ODTManager support collection requires geometry");

  SupportData data;
  const auto& support = line_geom.supportCells();
  data.ordered_samples.reserve(support.size());
  const amrex::IntVect owner = line_geom.ownerCell();
  const int dir = line_geom.dir();

  for (int i = 0; i < static_cast<int>(support.size()); ++i) {
    const auto& seg = support[static_cast<std::size_t>(i)];
    SupportSample sample;
    sample.support_ordinal = i;
    sample.iv = owner;
    sample.iv[dir] = seg.host_cell_index;

    ConservativeCell c{};
    if (!geom.Domain().contains(sample.iv)) {
      sample.provenance = SupportProvenance::PhysicalBoundaryGhost;
      sample.has_value = false;
    } else if (findCellAverageInMF(sample.iv, state_valid, true, c)) {
      sample.provenance = SupportProvenance::SameLevelValidCell;
      sample.value = c;
      sample.has_value = true;
    } else if (findCellAverageInMF(sample.iv, state_same_level, false, c)) {
      sample.provenance = SupportProvenance::SameLevelFilledGhost;
      sample.value = c;
      sample.has_value = true;
    } else if (findCellAverageInMF(sample.iv, state_host_filled, false, c)) {
      sample.provenance = SupportProvenance::AMRCoarseFineFilled;
      sample.value = c;
      sample.has_value = true;
    } else {
      sample.provenance = SupportProvenance::InvalidUnsupported;
      sample.has_value = false;
    }

    data.ordered_samples.push_back(sample);
  }

  return data;
}

ODTLineState&
ODTManager::initializeLineStateFromSupportData(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom,
  const SupportData& support_data)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    support_data.allSameLevelAccepted(),
    "ODTManager initialize requires same-level accepted support provenance");
  const auto local_support_cell_averages =
    support_data.conservativeAveragesOrdered();
  return initializeLineStateFromLESSupportAverages(
    level, owner_cell, dir, geom, local_support_cell_averages);
}

ODTLineState&
ODTManager::reconcileLineStateFromSupportData(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom,
  const SupportData& support_data)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    support_data.allSameLevelAccepted(),
    "ODTManager reconcile requires same-level accepted support provenance");

  assertLevelScope(level);
  auto* entry = findLineEntry(level, owner_cell, dir);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    entry != nullptr,
    "ODTManager reconcile requires an existing persistent line entry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    entry->state.initialized() && entry->state.valid(),
    "ODTManager reconcile requires an existing valid line state");
  AMREX_ALWAYS_ASSERT(entry->geometry.isDefined());
  AMREX_ALWAYS_ASSERT(entry->geometry.level() == level);
  AMREX_ALWAYS_ASSERT(entry->geometry.ownerCell() == owner_cell);
  AMREX_ALWAYS_ASSERT(entry->geometry.dir() == dir);

  const auto owner_avg = support_data.ownerAverage(entry->geometry);
  ODTReconcile::reconcileExistingLineStateToOwnerAverage(
    entry->geometry, owner_avg, entry->state);
  amrex::ignore_unused(geom);
  return entry->state;
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
    assembleLocalSupportCellAverages(entry.geometry, geom, state);
  return initializeLineStateFromLESSupportAverages(
    level, owner_cell, dir, geom, local_support_cell_averages);
}

ODTLineState&
ODTManager::reconcileLineStateToLESState(
  int level,
  const amrex::IntVect& owner_cell,
  int dir,
  const amrex::Geometry& geom,
  const amrex::MultiFab& state)
{
  assertLevelScope(level);
  auto* entry = findLineEntry(level, owner_cell, dir);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    entry != nullptr,
    "ODTManager reconcile requires an existing persistent line entry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    entry->state.initialized() && entry->state.valid(),
    "ODTManager reconcile requires an existing valid line state");
  AMREX_ALWAYS_ASSERT(entry->geometry.isDefined());
  AMREX_ALWAYS_ASSERT(entry->geometry.level() == level);
  AMREX_ALWAYS_ASSERT(entry->geometry.ownerCell() == owner_cell);
  AMREX_ALWAYS_ASSERT(entry->geometry.dir() == dir);

  const auto local_support_cell_averages =
    assembleLocalSupportCellAverages(entry->geometry, geom, state);
  SupportData support_data{};
  support_data.ordered_samples.resize(local_support_cell_averages.size());
  for (int i = 0; i < static_cast<int>(local_support_cell_averages.size()); ++i) {
    auto& sample = support_data.ordered_samples[static_cast<std::size_t>(i)];
    sample.support_ordinal = i;
    sample.value = local_support_cell_averages[static_cast<std::size_t>(i)];
    sample.has_value = true;
    sample.provenance = SupportProvenance::SameLevelValidCell;
  }
  const auto owner_avg = support_data.ownerAverage(entry->geometry);
  ODTReconcile::reconcileExistingLineStateToOwnerAverage(
    entry->geometry, owner_avg, entry->state);
  return entry->state;
}

std::vector<ODTLineState::ConservativeCell>
ODTManager::assembleLocalSupportCellAverages(
  const ODTLineGeometry& line_geom,
  const amrex::Geometry& geom,
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
    local_support_cell_averages[static_cast<std::size_t>(i)] =
      getConservativeCellAverageAt(iv, geom, state);
  }

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    static_cast<int>(local_support_cell_averages.size()) ==
      line_geom.supportCellCount(),
    "ODTManager support-cell average assembly size mismatch");

  return local_support_cell_averages;
}

ODTLineState::ConservativeCell
ODTManager::getConservativeCellAverageAt(
  const amrex::IntVect& iv,
  const amrex::Geometry& geom,
  const amrex::MultiFab& state) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.Domain().contains(iv),
    "ODTManager support-cell collection requires physical LES support indices "
    "inside the level domain");

  ODTLineState::ConservativeCell c{};
  const bool found = findCellAverageInMF(iv, state, false, c);

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    found,
    "ODTManager could not locate support cell in same-level host state; "
    "AMR/coarse-fine support provenance is not accepted in this phase");
  return c;
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

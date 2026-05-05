#ifndef PELEC_ODTLES_ODTMANAGER_H
#define PELEC_ODTLES_ODTMANAGER_H

#include <cstdint>
#include <cstddef>
#include <map>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include "ODTLineGeometry.H"
#include "ODTLineState.H"
#include "ODTParams.H"

namespace pelec::odtles
{

class ODTManager
{
public:
  using ConservativeCell = ODTLineState::ConservativeCell;

  struct LineEntry
  {
    ODTLineGeometry geometry;
    ODTLineState state;
    std::uint64_t runtime_event_seed = 0ULL;
    bool runtime_event_seed_initialized = false;
  };

  enum class SupportProvenance : int
  {
    SameLevelValidCell = 0,
    SameLevelFilledGhost = 1,
    PhysicalBoundaryGhost = 2,
    AMRCoarseFineFilled = 3,
    InvalidUnsupported = 4
  };

  struct SupportSample
  {
    int support_ordinal = -1;
    amrex::IntVect iv{AMREX_D_DECL(0, 0, 0)};
    SupportProvenance provenance = SupportProvenance::InvalidUnsupported;
    ConservativeCell value{};
    bool has_value = false;
  };

  struct SupportData
  {
    std::vector<SupportSample> ordered_samples;

    int count(SupportProvenance p) const noexcept;
    bool allSameLevelAccepted() const noexcept;
    bool allRuntimeAccepted(bool allow_amr_coarse_fine_support) const noexcept;
    std::vector<ConservativeCell> conservativeAveragesOrdered() const;
    ConservativeCell ownerAverage(const ODTLineGeometry& geom) const;
  };

  ODTManager() = default;
  explicit ODTManager(const ODTParams& params);

  void setParams(const ODTParams& params);
  const ODTParams& params() const noexcept { return m_params; }

  void initializeLevel(
    int level,
    const amrex::BoxArray& grids,
    const amrex::DistributionMapping& dmap);
  bool isInitialized() const noexcept { return m_initialized; }
  int level() const noexcept { return m_level; }

  void clear();

  bool hasLineEntry(int level, const amrex::IntVect& owner_cell, int dir) const;
  LineEntry*
  findLineEntry(int level, const amrex::IntVect& owner_cell, int dir);
  const LineEntry*
  findLineEntry(int level, const amrex::IntVect& owner_cell, int dir) const;
  LineEntry& getOrCreateLineEntry(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom);

  std::size_t numLineEntries() const noexcept { return m_line_entries.size(); }

  ODTLineGeometry& getOrCreateLineGeometry(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom);
  ODTLineGeometry*
  findLineGeometry(int level, const amrex::IntVect& owner_cell, int dir);
  const ODTLineGeometry* findLineGeometry(
    int level,
    const amrex::IntVect& owner_cell,
    int dir) const;

  std::size_t numLineGeometries() const noexcept
  {
    return m_line_entries.size();
  }

  ODTLineState&
  getOrCreateLineState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom);
  SupportData collectSupportData(
    const ODTLineGeometry& line_geom,
    const amrex::Geometry& geom,
    const amrex::MultiFab& state_valid,
    const amrex::MultiFab& state_same_level,
    const amrex::MultiFab& state_host_filled) const;
  ODTLineState& initializeLineStateFromSupportData(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const SupportData& support_data);
  ODTLineState& reconcileLineStateFromSupportData(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const SupportData& support_data);
  // Host-side input-contract path: `state` must carry same-level LES valid
  // cells plus same-level ghost copies (for inter-rank/periodic coverage),
  // not AMR coarse-fine ghost surrogates for support means.
  ODTLineState& initializeLineStateFromLESState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const amrex::MultiFab& state);
  ODTLineState& reconcileLineStateToLESState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const amrex::MultiFab& state);
  ODTLineState* findLineState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir);
  const ODTLineState* findLineState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir) const;

  std::size_t numLineStates() const noexcept { return m_line_entries.size(); }

private:
  struct LineKey
  {
    int level = -1;
    amrex::IntVect owner_cell{AMREX_D_DECL(0, 0, 0)};
    int dir = 0;

    bool operator<(const LineKey& other) const noexcept;
  };

  void assertLevelScope(int level) const;
  static bool findCellAverageInMF(
    const amrex::IntVect& iv,
    const amrex::MultiFab& state,
    bool use_validbox_only,
    ConservativeCell& out_cell);
  std::vector<ODTLineState::ConservativeCell> assembleLocalSupportCellAverages(
    const ODTLineGeometry& line_geom,
    const amrex::Geometry& geom,
    const amrex::MultiFab& state) const;
  ODTLineState::ConservativeCell getConservativeCellAverageAt(
    const amrex::IntVect& iv,
    const amrex::Geometry& geom,
    const amrex::MultiFab& state) const;
  ODTLineState& initializeLineStateFromLESSupportAverages(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const std::vector<ODTLineState::ConservativeCell>&
      local_support_cell_averages);

  ODTParams m_params{};
  bool m_initialized = false;
  int m_level = -1;
  std::map<LineKey, LineEntry> m_line_entries;
};

} // namespace pelec::odtles

#endif

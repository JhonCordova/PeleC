#ifndef PELEC_ODTLES_ODTRUNTIMELESBRIDGE_H
#define PELEC_ODTLES_ODTRUNTIMELESBRIDGE_H

#include <functional>

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include "ODTManager.H"

namespace pelec::odtles
{

struct RuntimeDepositionStats
{
  long rejected_boundary_entries = 0;
  long rejected_amr_entries = 0;
  long rejected_invalid_entries = 0;
  long rejected_mixed_entries = 0;
  long stepped_entries = 0;
  long stepper_attempted_events = 0;
  long stepper_applied_events = 0;
  long stepper_rejected_events = 0;
  long stepper_diffusion_only_catchup = 0;
  long moment_columns_built = 0;
};

using FluxRegisterUpdater = std::function<void(
  const amrex::MFIter&,
  amrex::FabType,
  const amrex::Array<amrex::FArrayBox const*, AMREX_SPACEDIM>&)>;

RuntimeDepositionStats
accumulateRuntimeODTMomentumLESTerm(
  int level,
  const amrex::Geometry& geom,
  const amrex::BoxArray& grids,
  const amrex::DistributionMapping& dmap,
  const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& area,
  const amrex::MultiFab& volume,
  ODTManager& odt_manager,
  const amrex::MultiFab& state_valid,
  const amrex::MultiFab& state_same_level,
  const amrex::MultiFab& state_host_filled,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  const FluxRegisterUpdater& reflux_updater = FluxRegisterUpdater{});

} // namespace pelec::odtles

#endif

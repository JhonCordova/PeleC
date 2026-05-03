#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <array>
#include <fstream>
#include <iomanip>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include "Diffterm.H"
#include "IndexDefines.H"
#include "ODTDiffusion.H"
#include "ODTMomentExtractor.H"
#include "ODTReconcile.H"
#include "ODTRuntimeLESBridge.H"
#include "ODTStepper.H"
#include "ODTThermoBridge.H"
#include "ODTTripletMap.H"
#include "PeleC.H"
#include "PelePhysics.H"

namespace pelec_tests
{
namespace
{

using ConservativeCell = pelec::odtles::ODTLineState::ConservativeCell;

ConservativeCell
makeCell(amrex::Real rho, amrex::Real rhou, amrex::Real rhov, amrex::Real rhow, amrex::Real rhoE)
{
  ConservativeCell c{};
  c.rho = rho;
  c.rhou = rhou;
  c.rhov = rhov;
  c.rhow = rhow;
  c.rhoE = rhoE;
  return c;
}

void
setSpeciesSequence(ConservativeCell& c, amrex::Real scale, amrex::Real shift = 0.0)
{
  if (c.rhoY.size() != static_cast<std::size_t>(NUM_SPECIES)) {
    c.rhoY.resize(static_cast<std::size_t>(NUM_SPECIES), 0.0);
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    c.rhoY[static_cast<std::size_t>(n)] =
      shift + scale * static_cast<amrex::Real>(n + 1);
  }
}

std::array<amrex::Real, 5>
asArray(const ConservativeCell& c)
{
  return {c.rho, c.rhou, c.rhov, c.rhow, c.rhoE};
}

int
supportIndexForOffset(const pelec::odtles::ODTLineGeometry& geom, int offset)
{
  const auto& segs = geom.supportCells();
  const amrex::Real target =
    static_cast<amrex::Real>(offset) * geom.deltaS();
  int best = -1;
  amrex::Real best_dist = 0.0;
  for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
    const auto& seg = segs[static_cast<std::size_t>(i)];
    if (seg.relative_offset != offset) {
      continue;
    }
    const amrex::Real dist = std::abs(seg.s_center - target);
    if (best < 0 || dist < best_dist) {
      best = i;
      best_dist = dist;
    }
  }
  return best;
}

std::vector<int>
supportIndicesForOffset(const pelec::odtles::ODTLineGeometry& geom, int offset)
{
  std::vector<int> idxs;
  const auto& segs = geom.supportCells();
  for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
    if (segs[static_cast<std::size_t>(i)].relative_offset == offset) {
      idxs.push_back(i);
    }
  }
  return idxs;
}

void
setSupportCellsForOffset(
  const pelec::odtles::ODTLineGeometry& geom,
  std::vector<ConservativeCell>& support,
  int offset,
  const ConservativeCell& c)
{
  const auto idxs = supportIndicesForOffset(geom, offset);
  for (int i : idxs) {
    support[static_cast<std::size_t>(i)] = c;
  }
}

ConservativeCell
ownerIntervalMean(
  const pelec::odtles::ODTLineGeometry& geom,
  const pelec::odtles::ODTLineState& line_state)
{
  const auto owner_interval = geom.ownerInterval();
  const amrex::Real owner_len = owner_interval.hi - owner_interval.lo;
  AMREX_ALWAYS_ASSERT(owner_len > 0.0);
  ConservativeCell out{};
  out.rho = 0.0;
  out.rhou = 0.0;
  out.rhov = 0.0;
  out.rhow = 0.0;
  out.rhoE = 0.0;
  if (out.rhoY.size() != static_cast<std::size_t>(NUM_SPECIES)) {
    out.rhoY.resize(static_cast<std::size_t>(NUM_SPECIES), 0.0);
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    out.rhoY[static_cast<std::size_t>(n)] = 0.0;
  }

  amrex::Real sum_w = 0.0;
  const auto& segs = geom.supportCells();
  for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
    const auto& seg = segs[static_cast<std::size_t>(i)];
    const amrex::Real overlap = std::max<amrex::Real>(
      0.0,
      std::min(owner_interval.hi, seg.s_interval.hi) -
        std::max(owner_interval.lo, seg.s_interval.lo));
    if (overlap <= 0.0) {
      continue;
    }
    const amrex::Real w = overlap / owner_len;
    sum_w += w;
    const auto& c = line_state.cell(i);
    out.rho += w * c.rho;
    out.rhou += w * c.rhou;
    out.rhov += w * c.rhov;
    out.rhow += w * c.rhow;
    out.rhoE += w * c.rhoE;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      out.rhoY[static_cast<std::size_t>(n)] +=
        w * c.rhoY[static_cast<std::size_t>(n)];
    }
  }
  AMREX_ALWAYS_ASSERT(sum_w > 0.0);

  return out;
}

void
setConservativeAt(
  amrex::MultiFab& mf, const amrex::IntVect& iv, const ConservativeCell& c)
{
  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
    if (!mfi.fabbox().contains(iv)) {
      continue;
    }
    auto const arr = mf.array(mfi);
    arr(iv, URHO) = c.rho;
    arr(iv, UMX) = c.rhou;
    arr(iv, UMY) = c.rhov;
    arr(iv, UMZ) = c.rhow;
    arr(iv, UEDEN) = c.rhoE;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      arr(iv, UFS + n) = c.rhoY[static_cast<std::size_t>(n)];
    }
    return;
  }
}

amrex::Real
specificInternalEnergy(const ConservativeCell& c)
{
  if (c.rho <= 0.0) {
    return -1.0;
  }
  const amrex::Real inv_rho = 1.0 / c.rho;
  const amrex::Real kinetic =
    0.5 * inv_rho * (c.rhou * c.rhou + c.rhov * c.rhov + c.rhow * c.rhow);
  return (c.rhoE - kinetic) * inv_rho;
}

amrex::Real
hostDynamicViscosityFromRecoveredState(
  const pelec::odtles::ODTThermoBridge::RecoveredState& rec)
{
  const bool get_xi = false;
  const bool get_mu = true;
  const bool get_lam = false;
  const bool get_Ddiag = false;
  const bool get_chi = false;
  amrex::Real mu = 0.0;
  amrex::Real xi_dummy = 0.0;
  amrex::Real lam_dummy = 0.0;
  std::vector<amrex::Real> Y = rec.Y;
  auto trans = pele::physics::PhysicsType::transport();
  auto const* tparm = &PeleC::trans_parms.host_parm();
  trans.transport(
    get_xi, get_mu, get_lam, get_Ddiag, get_chi, rec.temperature, rec.rho,
    Y.data(), nullptr, nullptr, mu, xi_dummy, lam_dummy, tparm);
  return mu;
}

amrex::Real
maxAbsMomentumComponentChange(
  const pelec::odtles::ODTLineState& before,
  const pelec::odtles::ODTLineState& after,
  int comp)
{
  AMREX_ALWAYS_ASSERT(before.numCells() == after.numCells());
  amrex::Real max_abs = 0.0;
  for (int i = 0; i < before.numCells(); ++i) {
    amrex::Real q_before = 0.0;
    amrex::Real q_after = 0.0;
    if (comp == static_cast<int>(pelec::odtles::ODTDiffusion::Component::RhoU)) {
      q_before = before.cell(i).rhou;
      q_after = after.cell(i).rhou;
    } else if (
      comp == static_cast<int>(pelec::odtles::ODTDiffusion::Component::RhoV)) {
      q_before = before.cell(i).rhov;
      q_after = after.cell(i).rhov;
    } else if (
      comp == static_cast<int>(pelec::odtles::ODTDiffusion::Component::RhoW)) {
      q_before = before.cell(i).rhow;
      q_after = after.cell(i).rhow;
    } else {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "Momentum component expected");
    }
    const amrex::Real d = std::abs(q_after - q_before);
    if (d > max_abs) {
      max_abs = d;
    }
  }
  return max_abs;
}

void
buildOdtDirectionalFluxes(
  const amrex::Box& cbox,
  const amrex::Box& domain,
  const amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& odt_dir_cc,
  const amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& area,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& flux_ec)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    const auto ebox = amrex::surroundingNodes(cbox, dir);
    auto const qcc = odt_dir_cc[dir].const_array();
    auto const a = area[dir].const_array();
    auto const flx = flux_ec[dir].array();

    for (amrex::IntVect iv = ebox.smallEnd(); iv <= ebox.bigEnd();
         ebox.next(iv)) {
      const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);
      const bool has_hi = domain.contains(iv);
      const bool has_lo = domain.contains(ivm);

      const int i = iv[0];
      const int j = iv[1];
      const int k = iv[2];
      const int il = i - (dir == 0 ? 1 : 0);
      const int jl = j - (dir == 1 ? 1 : 0);
      const int kl = k - (dir == 2 ? 1 : 0);

      auto centered_face_avg = [&](int comp) -> amrex::Real {
        const amrex::Real q_hi = has_hi ? qcc(i, j, k, comp) : 0.0;
        const amrex::Real q_lo = has_lo ? qcc(il, jl, kl, comp) : 0.0;
        if (has_hi && has_lo) {
          return 0.5 * (q_hi + q_lo);
        }
        if (has_hi) {
          return q_hi;
        }
        if (has_lo) {
          return q_lo;
        }
        return 0.0;
      };

      flx(i, j, k, UMX) = a(i, j, k) * centered_face_avg(
        pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomX);
      flx(i, j, k, UMY) = a(i, j, k) * centered_face_avg(
        pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomY);
      flx(i, j, k, UMZ) = a(i, j, k) * centered_face_avg(
        pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomZ);
      flx(i, j, k, UEDEN) = 0.0;
    }
  }
}

void
computeDivFromFluxes(
  const amrex::Box& cbox,
  const amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& flux_ec,
  const amrex::FArrayBox& volume,
  amrex::FArrayBox& lterm)
{
  const auto fx = flux_ec[0].const_array();
  const auto fy = flux_ec[1].const_array();
  const auto fz = flux_ec[2].const_array();
  const auto vol = volume.const_array();
  const auto D = lterm.array();

  for (amrex::IntVect iv = cbox.smallEnd(); iv <= cbox.bigEnd(); cbox.next(iv)) {
    for (int n = 0; n < NVAR; ++n) {
      pc_flux_div(iv[0], iv[1], iv[2], n, fx, fy, fz, vol, D);
    }
  }
}

void
computeLtermFromCenteredContributionMF(
  const amrex::BoxArray& ba,
  const amrex::DistributionMapping& dm,
  const amrex::Box& domain,
  amrex::MultiFab& lterm_out)
{
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> odt_dir_cc;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    odt_dir_cc[dir].define(ba, dm, 4, 1);
    odt_dir_cc[dir].setVal(0.0);
  }

  // Controlled nonuniform directional contribution in dir=0 only.
  for (amrex::MFIter mfi(odt_dir_cc[0], false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    auto const q = odt_dir_cc[0].array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      const amrex::Real x = static_cast<amrex::Real>(iv[0]);
      q(iv, 0) = x;
      q(iv, 1) = -x;
      q(iv, 2) = 2.0 * x;
      q(iv, 3) = 0.0;
    }
  }
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    odt_dir_cc[dir].FillBoundary();
  }

  amrex::MultiFab volume(ba, dm, 1, 0);
  volume.setVal(1.0);
  lterm_out.define(ba, dm, NVAR, 0);
  lterm_out.setVal(0.0);

  for (amrex::MFIter mfi(lterm_out, false); mfi.isValid(); ++mfi) {
    const amrex::Box vbox = mfi.validbox();
    amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> area{
      AMREX_D_DECL(
        amrex::FArrayBox(amrex::surroundingNodes(vbox, 0), 1),
        amrex::FArrayBox(amrex::surroundingNodes(vbox, 1), 1),
        amrex::FArrayBox(amrex::surroundingNodes(vbox, 2), 1))};
    amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec{
      AMREX_D_DECL(
        amrex::FArrayBox(amrex::surroundingNodes(vbox, 0), NVAR),
        amrex::FArrayBox(amrex::surroundingNodes(vbox, 1), NVAR),
        amrex::FArrayBox(amrex::surroundingNodes(vbox, 2), NVAR))};

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      area[dir].setVal(1.0);
      flux_ec[dir].setVal(0.0);
      const auto ebox = amrex::surroundingNodes(vbox, dir);
      auto const qcc = odt_dir_cc[dir].const_array(mfi);
      auto const a = area[dir].const_array();
      auto const flx = flux_ec[dir].array();

      for (amrex::IntVect iv = ebox.smallEnd(); iv <= ebox.bigEnd();
           ebox.next(iv)) {
        const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);
        const bool has_hi = domain.contains(iv);
        const bool has_lo = domain.contains(ivm);

        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];
        const int il = i - (dir == 0 ? 1 : 0);
        const int jl = j - (dir == 1 ? 1 : 0);
        const int kl = k - (dir == 2 ? 1 : 0);

        auto centered_face_avg = [&](int comp) -> amrex::Real {
          const amrex::Real q_hi = has_hi ? qcc(i, j, k, comp) : 0.0;
          const amrex::Real q_lo = has_lo ? qcc(il, jl, kl, comp) : 0.0;
          if (has_hi && has_lo) {
            return 0.5 * (q_hi + q_lo);
          }
          if (has_hi) {
            return q_hi;
          }
          if (has_lo) {
            return q_lo;
          }
          return 0.0;
        };

        flx(i, j, k, UMX) = a(i, j, k) * centered_face_avg(
          pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomX);
        flx(i, j, k, UMY) = a(i, j, k) * centered_face_avg(
          pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomY);
        flx(i, j, k, UMZ) = a(i, j, k) * centered_face_avg(
          pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomZ);
        flx(i, j, k, UEDEN) = 0.0;
      }
    }

    const auto fx = flux_ec[0].const_array();
    const auto fy = flux_ec[1].const_array();
    const auto fz = flux_ec[2].const_array();
    const auto vol = volume.const_array(mfi);
    const auto D = lterm_out.array(mfi);
    for (amrex::IntVect iv = vbox.smallEnd(); iv <= vbox.bigEnd(); vbox.next(iv)) {
      for (int n = 0; n < NVAR; ++n) {
        pc_flux_div(iv[0], iv[1], iv[2], n, fx, fy, fz, vol, D);
      }
    }
  }
}

amrex::Real
valueAt(const amrex::MultiFab& mf, const amrex::IntVect& iv, int n)
{
  for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
    if (mfi.validbox().contains(iv)) {
      return mf.const_array(mfi)(iv, n);
    }
  }
  return 0.0;
}

void
fillLinearDirectionalContribution(
  amrex::FArrayBox& qfab,
  const amrex::Box& bx,
  amrex::Real dx)
{
  auto const q = qfab.array();
  for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
    const amrex::Real x = (static_cast<amrex::Real>(iv[0]) + 0.5) * dx;
    q(iv, 0) = x;
    q(iv, 1) = -x;
    q(iv, 2) = 2.0 * x;
    q(iv, 3) = 0.0;
  }
}

void
computeSingleLevelLinearDeposition(
  const amrex::Box& cbox,
  const amrex::Box& domain,
  amrex::Real dx,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& odt_dir_cc,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& flux_ec,
  amrex::FArrayBox& lterm)
{
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> area{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 0), 1),
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 1), 1),
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 2), 1))};
  amrex::FArrayBox volume(cbox, 1);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    area[dir].setVal(1.0);
    flux_ec[dir].setVal(0.0);
    odt_dir_cc[dir].setVal(0.0);
  }

  // Build only j=0 contribution profile; j=1,2 stay zero.
  fillLinearDirectionalContribution(odt_dir_cc[0], odt_dir_cc[0].box(), dx);
  volume.setVal(dx);
  lterm.setVal(0.0);

  buildOdtDirectionalFluxes(cbox, domain, odt_dir_cc, area, flux_ec);
  computeDivFromFluxes(cbox, flux_ec, volume, lterm);
}

struct RuntimeSignature
{
  amrex::Real sum_umx = 0.0;
  amrex::Real sum_umy = 0.0;
  amrex::Real sum_umz = 0.0;
  amrex::Real sum_ueden = 0.0;
  amrex::Real l1_umx = 0.0;
  amrex::Real l1_umy = 0.0;
  amrex::Real l1_umz = 0.0;
  amrex::Real l1_ueden = 0.0;
  amrex::Real l1_others = 0.0;
  long stepped_entries = 0;
  long moment_columns_built = 0;
};

RuntimeSignature
computeRuntimeODTLESTermSignature(
  amrex::MultiFab& lterm_out, const pelec::odtles::ODTParams& runtime_params_in)
{
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(31, 3, 3)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(32.0, 4.0, 4.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);

  amrex::BoxArray ba(domain);
  ba.maxSize(4);
  amrex::DistributionMapping dm(ba);

  constexpr int support_ng = pelec::odtles::ODTLineGeometry::supportGhostCells();
  amrex::MultiFab state_valid(ba, dm, NVAR, 0);
  amrex::MultiFab state_same_level(ba, dm, NVAR, support_ng);
  amrex::MultiFab state_host_filled(ba, dm, NVAR, support_ng);
  auto initializeBackgroundConservativeState = [](amrex::MultiFab& mf) {
    mf.setVal(0.0);
    mf.setVal(1.0, URHO, 1, mf.nGrow());
    mf.setVal(3.0, UEDEN, 1, mf.nGrow());
    for (int n = 0; n < NUM_SPECIES; ++n) {
      mf.setVal(
        1.0 / static_cast<amrex::Real>(NUM_SPECIES), UFS + n, 1, mf.nGrow());
    }
  };
  initializeBackgroundConservativeState(state_valid);
  initializeBackgroundConservativeState(state_same_level);
  initializeBackgroundConservativeState(state_host_filled);

  for (amrex::MFIter mfi(state_valid, false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    const auto s = state_valid.array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      const amrex::Real x = static_cast<amrex::Real>(iv[0]);
      const amrex::Real y = static_cast<amrex::Real>(iv[1]);
      const amrex::Real z = static_cast<amrex::Real>(iv[2]);
      s(iv, URHO) = 1.0 + 0.01 * x;
      s(iv, UMX) = 0.15 * x + 0.02 * y + 0.5;
      s(iv, UMY) = -0.10 * x + 0.03 * z + 0.25;
      s(iv, UMZ) = 0.30 * x - 0.02 * y + 0.01 * z - 0.75;
      // Keep manufactured conservative state thermodynamically admissible
      // now that runtime molecular momentum diffusion always recovers (rho,T,Y).
      s(iv, UEDEN) = 80.0 + 0.20 * x + 0.10 * y;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s(iv, UFS + n) =
          (1.0 + 0.01 * x) * (0.01 * static_cast<amrex::Real>(n + 1));
      }
    }
  }
  amrex::MultiFab::Copy(state_same_level, state_valid, 0, 0, NVAR, 0);
  state_same_level.FillBoundary(geom.periodicity());
  amrex::MultiFab::Copy(
    state_host_filled, state_same_level, 0, 0, NVAR, state_same_level.nGrow());

  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> area{
    AMREX_D_DECL(
      amrex::MultiFab(amrex::convert(ba, amrex::IntVect::TheDimensionVector(0)), dm, 1, 0),
      amrex::MultiFab(amrex::convert(ba, amrex::IntVect::TheDimensionVector(1)), dm, 1, 0),
      amrex::MultiFab(amrex::convert(ba, amrex::IntVect::TheDimensionVector(2)), dm, 1, 0))};
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    area[dir].setVal(1.0);
  }
  amrex::MultiFab volume(ba, dm, 1, 0);
  volume.setVal(1.0);

  lterm_out.define(ba, dm, NVAR, 0);
  lterm_out.setVal(0.0);

  pelec::odtles::ODTManager odt_manager;
  pelec::odtles::ODTParams params = runtime_params_in;
  params.enabled = true;
  if (params.max_local_substeps <= 0) {
    params.max_local_substeps = 1;
  }
  odt_manager.setParams(params);

  const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> area_ptr = {
    AMREX_D_DECL(&area[0], &area[1], &area[2])};
  const auto stats = pelec::odtles::accumulateRuntimeODTMomentumLESTerm(
    0, geom, ba, dm, area_ptr, volume, odt_manager, state_valid, state_same_level,
    state_host_filled, 1.0e-2, lterm_out);

  RuntimeSignature sig{};
  for (amrex::MFIter mfi(lterm_out, false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    const auto arr = lterm_out.const_array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      const amrex::Real umx = arr(iv, UMX);
      const amrex::Real umy = arr(iv, UMY);
      const amrex::Real umz = arr(iv, UMZ);
      const amrex::Real ueden = arr(iv, UEDEN);
      sig.sum_umx += umx;
      sig.sum_umy += umy;
      sig.sum_umz += umz;
      sig.sum_ueden += ueden;
      sig.l1_umx += std::abs(umx);
      sig.l1_umy += std::abs(umy);
      sig.l1_umz += std::abs(umz);
      sig.l1_ueden += std::abs(ueden);
      sig.l1_others += std::abs(arr(iv, URHO));
      for (int n = UFS; n < UFS + NUM_SPECIES; ++n) {
        sig.l1_others += std::abs(arr(iv, n));
      }
#if NUM_AUX > 0
      for (int n = UFX; n < UFX + NUM_AUX; ++n) {
        sig.l1_others += std::abs(arr(iv, n));
      }
#endif
#if NUM_ADV > 0
      for (int n = UFA; n < UFA + NUM_ADV; ++n) {
        sig.l1_others += std::abs(arr(iv, n));
      }
#endif
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(sig.sum_umx);
  amrex::ParallelDescriptor::ReduceRealSum(sig.sum_umy);
  amrex::ParallelDescriptor::ReduceRealSum(sig.sum_umz);
  amrex::ParallelDescriptor::ReduceRealSum(sig.sum_ueden);
  amrex::ParallelDescriptor::ReduceRealSum(sig.l1_umx);
  amrex::ParallelDescriptor::ReduceRealSum(sig.l1_umy);
  amrex::ParallelDescriptor::ReduceRealSum(sig.l1_umz);
  amrex::ParallelDescriptor::ReduceRealSum(sig.l1_ueden);
  amrex::ParallelDescriptor::ReduceRealSum(sig.l1_others);
  sig.stepped_entries = stats.stepped_entries;
  sig.moment_columns_built = stats.moment_columns_built;
  return sig;
}

RuntimeSignature
computeRuntimeODTLESTermSignature(amrex::MultiFab& lterm_out)
{
  const pelec::odtles::ODTParams params{};
  return computeRuntimeODTLESTermSignature(lterm_out, params);
}

void
writeRuntimeSignatureIfRequested(const RuntimeSignature& sig)
{
  const char* out_path = std::getenv("ODTLES_MPI_SIGNATURE_FILE");
  if (out_path == nullptr || !amrex::ParallelDescriptor::IOProcessor()) {
    return;
  }
  std::ofstream out(out_path);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    out.good(), "Unable to open ODTLES_MPI_SIGNATURE_FILE for writing");
  out << std::setprecision(17) << std::scientific;
  out << "sum_umx " << sig.sum_umx << "\n";
  out << "sum_umy " << sig.sum_umy << "\n";
  out << "sum_umz " << sig.sum_umz << "\n";
  out << "sum_ueden " << sig.sum_ueden << "\n";
  out << "l1_umx " << sig.l1_umx << "\n";
  out << "l1_umy " << sig.l1_umy << "\n";
  out << "l1_umz " << sig.l1_umz << "\n";
  out << "l1_ueden " << sig.l1_ueden << "\n";
  out << "l1_others " << sig.l1_others << "\n";
}

amrex::Real
sumWeightedComponent(
  const amrex::FArrayBox& lterm,
  const amrex::FArrayBox& volume,
  const amrex::Box& bx,
  int comp)
{
  const auto d = lterm.const_array();
  const auto vol = volume.const_array();
  amrex::Real sum = 0.0;
  for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
    sum += vol(iv, 0) * d(iv, comp);
  }
  return sum;
}

} // namespace

TEST(ODTLESLocalEngine, DiffusionMVPValidation)
{
  constexpr amrex::Real tol = 1.0e-12;
  const auto rep = pelec::odtles::ODTDiffusion::runMVPValidationHook();

  EXPECT_TRUE(rep.uniform_invariant);
  EXPECT_TRUE(rep.nonuniform_evolved);
  EXPECT_LE(rep.uniform_max_abs_change, tol);
  EXPECT_GT(rep.nonuniform_max_abs_change, tol);
}

TEST(ODTLESLocalEngine, TripletMapKernelMVPValidation)
{
  constexpr amrex::Real tol = 1.0e-12;
  const auto rep = pelec::odtles::ODTTripletMap::runMVPValidationHook();

  EXPECT_TRUE(rep.uniform_invariant);
  EXPECT_TRUE(rep.nonuniform_modified);
  EXPECT_TRUE(rep.map_locality_preserved);
  EXPECT_TRUE(rep.kernel_zero_net_momentum);

  EXPECT_LE(rep.uniform_max_abs_change, tol);
  EXPECT_GT(rep.nonuniform_max_abs_change, tol);
  EXPECT_LE(rep.outside_interval_max_abs_change, tol);
  EXPECT_LE(std::abs(rep.kernel_net_rhou_change), tol);
  EXPECT_LE(std::abs(rep.kernel_net_rhov_change), tol);
  EXPECT_LE(std::abs(rep.kernel_net_rhow_change), tol);
}

TEST(ODTLESLocalEngine, StepperMVPValidation)
{
  constexpr amrex::Real tol = 1.0e-12;
  const auto rep = pelec::odtles::ODTStepper::runMVPValidationHook();

  EXPECT_TRUE(rep.closes_exact_dt);
  EXPECT_TRUE(rep.no_overshoot);
  EXPECT_TRUE(rep.no_event_path_handled);
  EXPECT_TRUE(rep.event_path_attempted);
  EXPECT_TRUE(rep.compatibility_policy_exercised);
  EXPECT_TRUE(rep.deterministic_reproducible);
  EXPECT_TRUE(rep.metadata_consistent);

  EXPECT_LE(rep.closure_error, tol);
  EXPECT_LE(rep.overshoot_amount, tol);
}

TEST(ODTLESLocalEngine, MomentExtractorMVPValidation)
{
  constexpr amrex::Real tol = 1.0e-12;
  const auto rep = pelec::odtles::ODTMomentExtractor::runMVPValidationHook();

  EXPECT_TRUE(rep.overlap_weights_partition_unity);
  EXPECT_TRUE(rep.owner_central_recovered);
  EXPECT_FALSE(rep.tau_zero_for_single_segment_central_interval);
  EXPECT_LE(rep.overlap_weight_sum_error, tol);
  EXPECT_LE(rep.owner_recovery_max_abs_error, tol);
  EXPECT_GT(rep.tau_max_abs, tol);
}

TEST(ODTLESLocalEngine, OwnerIntervalHasMultipleContributingSegments)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  pelec::odtles::ODTLineState line_state;
  line_state.initialize(line_geom);
  for (int i = 0; i < line_state.numCells(); ++i) {
    line_state.setCell(i, makeCell(1.0, 0.25, -0.10, 0.05, 3.0));
  }
  line_state.setValid(true);

  const auto col =
    pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
      line_geom, line_state);
  ASSERT_TRUE(col.valid);
  EXPECT_EQ(
    col.contributing_segments,
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);
  EXPECT_GT(col.contributing_segments, 1);
  EXPECT_NEAR(col.overlap_weight_sum, 1.0, tol);
}

TEST(ODTLESLocalEngine, InteriorGeometryHasNineSupportAndThreeOwnerSubsegments)
{
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(6, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(7.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(3, 0, 0));

  const pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  EXPECT_EQ(
    line_geom.supportCellCount(),
    pelec::odtles::ODTLineGeometry::HostSupportCells *
      pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);
  EXPECT_EQ(
    static_cast<int>(supportIndicesForOffset(line_geom, -1).size()),
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);
  EXPECT_EQ(
    static_cast<int>(supportIndicesForOffset(line_geom, 0).size()),
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);
  EXPECT_EQ(
    static_cast<int>(supportIndicesForOffset(line_geom, 1).size()),
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);
}

TEST(ODTLESLocalEngine, InteriorGeometryWithFiveSubsegmentsHasExpectedSupport)
{
  constexpr int n_sub = 5;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(6, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(7.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(3, 0, 0));

  const pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0, n_sub);
  EXPECT_EQ(line_geom.subsegmentsPerHostCell(), n_sub);
  EXPECT_EQ(
    line_geom.supportCellCount(),
    pelec::odtles::ODTLineGeometry::HostSupportCells * n_sub);
  EXPECT_EQ(static_cast<int>(supportIndicesForOffset(line_geom, -1).size()), n_sub);
  EXPECT_EQ(static_cast<int>(supportIndicesForOffset(line_geom, 0).size()), n_sub);
  EXPECT_EQ(static_cast<int>(supportIndicesForOffset(line_geom, 1).size()), n_sub);
}

TEST(
  ODTLESLocalEngine,
  ManufacturedOwnerSubcellVariationChangesExtractedMoments)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  const auto owner_idxs = supportIndicesForOffset(line_geom, 0);
  ASSERT_EQ(
    static_cast<int>(owner_idxs.size()),
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);

  pelec::odtles::ODTLineState uniform_state;
  uniform_state.initialize(line_geom);
  for (int i = 0; i < uniform_state.numCells(); ++i) {
    uniform_state.setCell(i, makeCell(1.0, 1.0, 0.0, 0.0, 3.0));
  }
  uniform_state.setValid(true);

  auto varied_state = uniform_state;
  const int n_owner = static_cast<int>(owner_idxs.size());
  const amrex::Real center = 0.5 * static_cast<amrex::Real>(n_owner - 1);
  for (int k = 0; k < n_owner; ++k) {
    // Preserve owner-interval mean rho*u while introducing subcell variance.
    varied_state.cell(owner_idxs[static_cast<std::size_t>(k)]).rhou =
      1.0 + 0.2 * (static_cast<amrex::Real>(k) - center);
  }

  const auto col_uniform =
    pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
      line_geom, uniform_state);
  const auto col_varied =
    pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
      line_geom, varied_state);
  ASSERT_TRUE(col_uniform.valid);
  ASSERT_TRUE(col_varied.valid);
  EXPECT_EQ(
    col_varied.contributing_segments,
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);

  EXPECT_NEAR(col_uniform.mean_rho, col_varied.mean_rho, tol);
  EXPECT_NEAR(col_uniform.mean_rho_ui[0], col_varied.mean_rho_ui[0], tol);
  EXPECT_NEAR(col_uniform.tau_ij[0], 0.0, tol);
  EXPECT_GT(std::abs(col_varied.tau_ij[0]), tol);
  EXPECT_GT(
    std::abs(col_varied.mean_rho_ui_uj[0] - col_uniform.mean_rho_ui_uj[0]),
    tol);
}

TEST(
  ODTLESLocalEngine,
  ManufacturedOwnerSubcellVariationChangesExtractedMomentsWithFiveSubsegments)
{
  constexpr amrex::Real tol = 1.0e-12;
  constexpr int n_sub = 5;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0, n_sub);
  const auto owner_idxs = supportIndicesForOffset(line_geom, 0);
  ASSERT_EQ(static_cast<int>(owner_idxs.size()), n_sub);

  pelec::odtles::ODTLineState uniform_state;
  uniform_state.initialize(line_geom);
  for (int i = 0; i < uniform_state.numCells(); ++i) {
    uniform_state.setCell(i, makeCell(1.0, 1.0, 0.0, 0.0, 3.0));
  }
  uniform_state.setValid(true);

  auto varied_state = uniform_state;
  const amrex::Real center = 0.5 * static_cast<amrex::Real>(n_sub - 1);
  for (int k = 0; k < n_sub; ++k) {
    varied_state.cell(owner_idxs[static_cast<std::size_t>(k)]).rhou =
      1.0 + 0.15 * (static_cast<amrex::Real>(k) - center);
  }

  const auto col_uniform =
    pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
      line_geom, uniform_state);
  const auto col_varied =
    pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
      line_geom, varied_state);
  ASSERT_TRUE(col_uniform.valid);
  ASSERT_TRUE(col_varied.valid);
  EXPECT_EQ(col_varied.contributing_segments, n_sub);
  EXPECT_NEAR(col_uniform.mean_rho, col_varied.mean_rho, tol);
  EXPECT_NEAR(col_uniform.mean_rho_ui[0], col_varied.mean_rho_ui[0], tol);
  EXPECT_NEAR(col_uniform.tau_ij[0], 0.0, tol);
  EXPECT_GT(std::abs(col_varied.tau_ij[0]), tol);
}

TEST(ODTLESLocalEngine, ODTLineStateExtendedSpeciesConstructionResetCopy)
{
  pelec::odtles::ODTLineState line_state;
  line_state.initializeForValidation(2);
  ASSERT_TRUE(line_state.initialized());
  ASSERT_EQ(line_state.numCells(), 2);

  auto c0 = line_state.cell(0);
  c0.rho = 2.0;
  c0.rhou = 1.0;
  c0.rhov = -0.5;
  c0.rhow = 0.25;
  c0.rhoE = 8.0;
  ASSERT_EQ(c0.rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    c0.rhoY[static_cast<std::size_t>(n)] = 0.1 * static_cast<amrex::Real>(n + 1);
  }
  line_state.setCell(0, c0);
  line_state.setValid(true);

  const auto vel = line_state.velocities(0);
  EXPECT_DOUBLE_EQ(vel[0], 0.5);
#if AMREX_SPACEDIM >= 2
  EXPECT_DOUBLE_EQ(vel[1], -0.25);
#endif
#if AMREX_SPACEDIM == 3
  EXPECT_DOUBLE_EQ(vel[2], 0.125);
#endif

  pelec::odtles::ODTLineState copied = line_state;
  ASSERT_TRUE(copied.initialized());
  ASSERT_TRUE(copied.valid());
  ASSERT_EQ(copied.numCells(), line_state.numCells());
  ASSERT_EQ(
    copied.cell(0).rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_DOUBLE_EQ(
      copied.cell(0).rhoY[static_cast<std::size_t>(n)],
      line_state.cell(0).rhoY[static_cast<std::size_t>(n)]);
  }

  line_state.reset();
  EXPECT_FALSE(line_state.initialized());
  EXPECT_FALSE(line_state.valid());
  EXPECT_EQ(line_state.numCells(), 0);
}

TEST(ODTLESLocalEngine, ODTLineStateSpeciesMassFractionRecovery)
{
  constexpr amrex::Real tol = 1.0e-12;
  pelec::odtles::ODTLineState line_state;
  line_state.initializeForValidation(1);

  ConservativeCell c{};
  c.rho = 2.5;
  c.rhou = 0.75;
  c.rhov = -0.25;
  c.rhow = 0.5;
  c.rhoE = 9.0;
  ASSERT_EQ(c.rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    c.rhoY[static_cast<std::size_t>(n)] =
      0.05 * static_cast<amrex::Real>(n + 1);
  }
  line_state.setCell(0, c);

  const auto Y = line_state.speciesMassFractions(0);
  ASSERT_EQ(Y.size(), static_cast<std::size_t>(NUM_SPECIES));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real expected =
      c.rhoY[static_cast<std::size_t>(n)] / c.rho;
    EXPECT_NEAR(line_state.speciesMassFraction(0, n), expected, tol);
    EXPECT_NEAR(Y[static_cast<std::size_t>(n)], expected, tol);
  }
}

TEST(ODTLESLocalEngine, ODTLineStateAdmissibilityChecksDensityAndSpecies)
{
  pelec::odtles::ODTLineState line_state;
  line_state.initializeForValidation(3);

  ConservativeCell admissible{};
  admissible.rho = 1.0;
  admissible.rhou = 0.2;
  admissible.rhov = -0.1;
  admissible.rhow = 0.05;
  admissible.rhoE = 4.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    admissible.rhoY[static_cast<std::size_t>(n)] =
      0.01 * static_cast<amrex::Real>(n + 1);
  }
  line_state.setCell(0, admissible);
  EXPECT_TRUE(line_state.hasExpectedSpeciesContainerSize(0));
  EXPECT_TRUE(line_state.isCellAdmissible(0));

  ConservativeCell neg_rho = admissible;
  neg_rho.rho = -1.0;
  line_state.setCell(1, neg_rho);
  EXPECT_FALSE(line_state.isCellAdmissible(1));

  ConservativeCell bad_species = admissible;
  if (NUM_SPECIES > 0) {
    bad_species.rhoY[0] = -1.0e-6;
  } else {
    bad_species.rhoY.push_back(0.0);
  }
  line_state.setCell(2, bad_species);
  EXPECT_FALSE(line_state.isCellAdmissible(2));

  ConservativeCell bad_size = admissible;
  bad_size.rhoY.resize(static_cast<std::size_t>(NUM_SPECIES + 1), 0.0);
  line_state.setCell(1, bad_size);
  EXPECT_FALSE(line_state.hasExpectedSpeciesContainerSize(1));
  EXPECT_FALSE(line_state.isCellAdmissible(1));
}

TEST(ODTLESLocalEngine, ODTThermoBridgeRecoversNormalizedSpeciesAndDerivedState)
{
  constexpr amrex::Real tol = 1.0e-12;
  pelec::odtles::ODTLineState line_state;
  line_state.initializeForValidation(1);

  ConservativeCell c{};
  c.rho = 2.0;
  c.rhou = 0.8;
  c.rhov = -0.2;
  c.rhow = 0.4;
  c.rhoE = 10.0;
  amrex::Real species_sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    c.rhoY[static_cast<std::size_t>(n)] = 0.05 * static_cast<amrex::Real>(n + 1);
    species_sum += c.rhoY[static_cast<std::size_t>(n)];
  }
  line_state.setCell(0, c);

  const auto rec = pelec::odtles::ODTThermoBridge::recoverCellThermoState(
    line_state, 0);
  ASSERT_TRUE(rec.success)
    << pelec::odtles::ODTThermoBridge::failureModeString(rec.failure_mode);

  EXPECT_NEAR(rec.velocity[0], c.rhou / c.rho, tol);
#if AMREX_SPACEDIM >= 2
  EXPECT_NEAR(rec.velocity[1], c.rhov / c.rho, tol);
#endif
#if AMREX_SPACEDIM == 3
  EXPECT_NEAR(rec.velocity[2], c.rhow / c.rho, tol);
#endif

  amrex::Real sum_rhoY = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real expected_Y =
      c.rhoY[static_cast<std::size_t>(n)] / species_sum;
    EXPECT_NEAR(rec.Y[static_cast<std::size_t>(n)], expected_Y, tol);
    EXPECT_NEAR(rec.rhoY[static_cast<std::size_t>(n)], c.rho * expected_Y, tol);
    sum_rhoY += rec.rhoY[static_cast<std::size_t>(n)];
  }
  EXPECT_NEAR(sum_rhoY, c.rho, tol);

  const amrex::Real u = c.rhou / c.rho;
  const amrex::Real v = c.rhov / c.rho;
  const amrex::Real w = c.rhow / c.rho;
  const amrex::Real expected_kinetic =
    0.5 * AMREX_D_TERM(u * u, +v * v, +w * w);
  const amrex::Real expected_e = c.rhoE / c.rho - expected_kinetic;
  EXPECT_NEAR(rec.kinetic_energy, expected_kinetic, tol);
  EXPECT_NEAR(rec.specific_internal_energy, expected_e, tol);
}

TEST(ODTLESLocalEngine, ODTThermoBridgeTemperatureMatchesHostEOS)
{
  constexpr amrex::Real tol = 1.0e-8;
  auto eos = pele::physics::PhysicsType::eos();

  ConservativeCell c{};
  c.rho = 1.7;
  constexpr amrex::Real T_ref = 420.0;
  amrex::Real Y_ref[NUM_SPECIES] = {0.0};
  amrex::Real ysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y_ref[n] = static_cast<amrex::Real>(n + 1);
    ysum += Y_ref[n];
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y_ref[n] /= ysum;
    c.rhoY[static_cast<std::size_t>(n)] = c.rho * Y_ref[n];
  }

  amrex::Real e_ref = 0.0;
  eos.RTY2E(c.rho, T_ref, Y_ref, e_ref);
  const amrex::Real u = 0.2;
  const amrex::Real v = -0.1;
  const amrex::Real w = 0.05;
  c.rhou = c.rho * u;
  c.rhov = c.rho * v;
  c.rhow = c.rho * w;
  c.rhoE = c.rho * (e_ref + 0.5 * (u * u + v * v + w * w));

  const auto rec = pelec::odtles::ODTThermoBridge::recoverCellThermoState(c);
  ASSERT_TRUE(rec.success)
    << pelec::odtles::ODTThermoBridge::failureModeString(rec.failure_mode);
  EXPECT_NEAR(rec.temperature, T_ref, tol);

  amrex::Real T_eos = 0.0;
  eos.REY2T(c.rho, rec.specific_internal_energy, rec.Y.data(), T_eos);
  EXPECT_NEAR(rec.temperature, T_eos, tol);
}

TEST(ODTLESLocalEngine, ODTThermoBridgeNegativeSpeciesRepairOrFail)
{
  ConservativeCell c{};
  c.rho = 1.5;
  c.rhou = 0.15;
  c.rhov = -0.03;
  c.rhow = 0.06;
  c.rhoE = 4.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    c.rhoY[static_cast<std::size_t>(n)] = 0.0;
  }
  if (NUM_SPECIES > 0) {
    c.rhoY[0] = -0.2;
  }
  if (NUM_SPECIES > 1) {
    c.rhoY[1] = 0.4;
  }
  for (int n = 2; n < NUM_SPECIES; ++n) {
    c.rhoY[static_cast<std::size_t>(n)] = 0.1;
  }

  const auto rec = pelec::odtles::ODTThermoBridge::recoverCellThermoState(c);
  if (NUM_SPECIES > 1) {
    ASSERT_TRUE(rec.success)
      << pelec::odtles::ODTThermoBridge::failureModeString(rec.failure_mode);
    EXPECT_GE(rec.rhoY[0], 0.0);
    EXPECT_EQ(rec.rhoY[0], 0.0);
  } else {
    EXPECT_FALSE(rec.success);
    EXPECT_EQ(
      rec.failure_mode,
      pelec::odtles::ODTThermoBridge::FailureMode::NonPositiveSpeciesSum);
  }
}

TEST(ODTLESLocalEngine, ODTThermoBridgeFailsOnNonPositiveSpeciesSum)
{
  ConservativeCell c{};
  c.rho = 1.1;
  c.rhou = 0.0;
  c.rhov = 0.0;
  c.rhow = 0.0;
  c.rhoE = 2.5;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    c.rhoY[static_cast<std::size_t>(n)] = -1.0e-3;
  }

  const auto rec = pelec::odtles::ODTThermoBridge::recoverCellThermoState(c);
  EXPECT_FALSE(rec.success);
  EXPECT_EQ(
    rec.failure_mode,
    pelec::odtles::ODTThermoBridge::FailureMode::NonPositiveSpeciesSum);
}

TEST(
  ODTLESLocalEngine,
  ODTDiffusionMolecularViscosityProfileMatchesHostTransportAndIsNonNegative)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  pelec::odtles::ODTLineState line_state;
  line_state.initialize(line_geom);

  auto eos = pele::physics::PhysicsType::eos();
  const int n_cells = line_state.numCells();
  for (int i = 0; i < n_cells; ++i) {
    ConservativeCell c{};
    c.rho = 1.0 + 0.2 * static_cast<amrex::Real>(i);

    amrex::Real Y_ref[NUM_SPECIES] = {0.0};
    amrex::Real ysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_ref[n] = static_cast<amrex::Real>(n + 1);
      ysum += Y_ref[n];
    }
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_ref[n] /= ysum;
      c.rhoY[static_cast<std::size_t>(n)] = c.rho * Y_ref[n];
    }

    const amrex::Real T_ref = 360.0 + 20.0 * static_cast<amrex::Real>(i);
    amrex::Real e_ref = 0.0;
    eos.RTY2E(c.rho, T_ref, Y_ref, e_ref);

    const amrex::Real u = 0.05 * static_cast<amrex::Real>(i + 1);
    const amrex::Real v = -0.03 * static_cast<amrex::Real>(i + 1);
    const amrex::Real w = 0.02 * static_cast<amrex::Real>(i + 1);
    c.rhou = c.rho * u;
    c.rhov = c.rho * v;
    c.rhow = c.rho * w;
    c.rhoE = c.rho * (e_ref + 0.5 * (u * u + v * v + w * w));
    line_state.setCell(i, c);
  }
  line_state.setValid(true);

  const auto profile =
    pelec::odtles::ODTDiffusion::buildMomentumViscosityProfile(line_state);
  ASSERT_TRUE(profile.success);
  ASSERT_EQ(
    static_cast<int>(profile.dynamic_viscosity_mu.size()), line_state.numCells());

  for (int i = 0; i < line_state.numCells(); ++i) {
    const auto rec = pelec::odtles::ODTThermoBridge::recoverCellThermoState(
      line_state, i);
    ASSERT_TRUE(rec.success)
      << pelec::odtles::ODTThermoBridge::failureModeString(rec.failure_mode);
    const amrex::Real mu_ref = hostDynamicViscosityFromRecoveredState(rec);
    EXPECT_NEAR(
      profile.dynamic_viscosity_mu[static_cast<std::size_t>(i)], mu_ref, tol);
    EXPECT_GE(profile.dynamic_viscosity_mu[static_cast<std::size_t>(i)], 0.0);
  }
}

TEST(
  ODTLESLocalEngine,
  ODTDiffusionDefaultControlsUseProductionMolecularMomentumPath)
{
  const pelec::odtles::ODTDiffusion::Controls ctrl{};
  EXPECT_TRUE(ctrl.use_molecular_viscosity_for_momentum);
  EXPECT_EQ(
    ctrl.momentum_molecular_form,
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::VelocityGradientMuFlux);
  EXPECT_TRUE(ctrl.fail_on_molecular_viscosity_recovery_failure);
}

TEST(
  ODTLESLocalEngine,
  RuntimeStepperControlPlumbingEnforcesProductionMolecularModel)
{
  pelec::odtles::ODTParams params{};
  pelec::odtles::ODTStepper::Controls step_ctrl{};
  pelec::odtles::configureRuntimeODTStepperControls(params, step_ctrl);

  EXPECT_EQ(step_ctrl.max_internal_iterations, 1);
  EXPECT_TRUE(step_ctrl.diffusion_controls.use_molecular_viscosity_for_momentum);
  EXPECT_EQ(
    step_ctrl.diffusion_controls.momentum_molecular_form,
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::VelocityGradientMuFlux);
  EXPECT_DOUBLE_EQ(step_ctrl.diffusion_controls.molecular_viscosity_floor, 0.0);
  EXPECT_TRUE(
    step_ctrl.diffusion_controls.fail_on_molecular_viscosity_recovery_failure);
}

TEST(
  ODTLESLocalEngine,
  RuntimeStepperControlPlumbingOnlyExposesSubstepCountTuning)
{
  pelec::odtles::ODTParams params{};
  params.max_local_substeps = 7;

  pelec::odtles::ODTStepper::Controls step_ctrl{};
  pelec::odtles::configureRuntimeODTStepperControls(params, step_ctrl);

  EXPECT_EQ(step_ctrl.max_internal_iterations, 7);
  EXPECT_TRUE(step_ctrl.diffusion_controls.use_molecular_viscosity_for_momentum);
  EXPECT_EQ(
    step_ctrl.diffusion_controls.momentum_molecular_form,
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::VelocityGradientMuFlux);
  EXPECT_DOUBLE_EQ(step_ctrl.diffusion_controls.molecular_viscosity_floor, 0.0);
  EXPECT_TRUE(
    step_ctrl.diffusion_controls.fail_on_molecular_viscosity_recovery_failure);
}

TEST(
  ODTLESLocalEngine,
  ODTDiffusionNonUniformDensityDistinguishesConservativeAndVelocityGradientForms)
{
  constexpr amrex::Real tol_quiet = 1.0e-12;
  constexpr amrex::Real tol_active = 1.0e-10;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  pelec::odtles::ODTLineState base_state;
  base_state.initialize(line_geom);
  auto eos = pele::physics::PhysicsType::eos();

  const amrex::Real u_const = 0.7;
  amrex::Real Y_ref[NUM_SPECIES] = {0.0};
  amrex::Real ysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y_ref[n] = static_cast<amrex::Real>(n + 1);
    ysum += Y_ref[n];
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y_ref[n] /= ysum;
  }

  for (int i = 0; i < base_state.numCells(); ++i) {
    ConservativeCell c{};
    const amrex::Real ii = static_cast<amrex::Real>(i);
    c.rho = 1.0 + 0.25 * ii * ii;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      c.rhoY[static_cast<std::size_t>(n)] = c.rho * Y_ref[n];
    }
    amrex::Real e_ref = 0.0;
    eos.RTY2E(c.rho, 430.0, Y_ref, e_ref);
    c.rhou = c.rho * u_const;
    c.rhov = 0.0;
    c.rhow = 0.0;
    c.rhoE = c.rho * (e_ref + 0.5 * u_const * u_const);
    base_state.setCell(i, c);
  }
  base_state.setValid(true);

  pelec::odtles::ODTDiffusion::Controls ctrl{};
  ctrl.diffuse_component.fill(false);
  ctrl.diffuse_component[static_cast<int>(pelec::odtles::ODTDiffusion::Component::RhoU)] =
    true;
  ctrl.use_molecular_viscosity_for_momentum = true;
  ctrl.molecular_viscosity_floor = 1.0;
  ctrl.fail_on_molecular_viscosity_recovery_failure = true;

  const auto mu_profile =
    pelec::odtles::ODTDiffusion::buildMomentumViscosityProfile(base_state);
  ASSERT_TRUE(mu_profile.success);
  amrex::Real max_mu = 0.0;
  for (const auto mu : mu_profile.dynamic_viscosity_mu) {
    if (mu > max_mu) {
      max_mu = mu;
    }
  }

  pelec::odtles::ODTLineState conservative_form_state = base_state;
  ctrl.momentum_molecular_form =
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::ConservativeVariableNuOnRhoU;
  pelec::odtles::ODTDiffusion::applyImplicitUniform(
    line_geom, conservative_form_state, 1.0e-2, ctrl);

  pelec::odtles::ODTLineState velocity_flux_state = base_state;
  ctrl.momentum_molecular_form =
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::VelocityGradientMuFlux;
  pelec::odtles::ODTDiffusion::applyImplicitUniform(
    line_geom, velocity_flux_state, 1.0e-2, ctrl);

  const int comp_rhou = static_cast<int>(pelec::odtles::ODTDiffusion::Component::RhoU);
  const amrex::Real conservative_change = maxAbsMomentumComponentChange(
    base_state, conservative_form_state, comp_rhou);
  const amrex::Real velocity_flux_change = maxAbsMomentumComponentChange(
    base_state, velocity_flux_state, comp_rhou);

  // u is uniform: velocity-gradient viscous flux should be quiescent.
  EXPECT_LE(velocity_flux_change, tol_quiet);
  // Conservative-variable nu_eff diffusion operates on rho*u and changes state.
  EXPECT_GT(conservative_change, tol_active)
    << "max_mu=" << max_mu << " (zero molecular viscosity gives zero update)";
}

TEST(ODTLESLocalEngine, SupportDataSamplingTransfersConservativeSpeciesFromLES)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  amrex::BoxArray ba(domain);
  amrex::DistributionMapping dm(ba);

  constexpr int support_ng = pelec::odtles::ODTLineGeometry::supportGhostCells();
  amrex::MultiFab state_valid(ba, dm, NVAR, 0);
  amrex::MultiFab state_same_level(ba, dm, NVAR, support_ng);
  amrex::MultiFab state_host_filled(ba, dm, NVAR, support_ng);
  state_valid.setVal(0.0);
  state_same_level.setVal(0.0);
  state_host_filled.setVal(0.0);

  for (amrex::MFIter mfi(state_valid, false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    auto const s = state_valid.array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      const amrex::Real x = static_cast<amrex::Real>(iv[0]);
      s(iv, URHO) = 1.0 + 0.1 * x;
      s(iv, UMX) = 0.2 + 0.03 * x;
      s(iv, UMY) = -0.1 + 0.02 * x;
      s(iv, UMZ) = 0.05 + 0.01 * x;
      s(iv, UEDEN) = 3.0 + 0.2 * x;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s(iv, UFS + n) =
          0.01 * static_cast<amrex::Real>(n + 1) * (1.0 + x);
      }
    }
  }
  amrex::MultiFab::Copy(state_same_level, state_valid, 0, 0, NVAR, 0);
  state_same_level.FillBoundary(geom.periodicity());
  amrex::MultiFab::Copy(
    state_host_filled, state_same_level, 0, 0, NVAR, state_same_level.nGrow());

  pelec::odtles::ODTManager odt_manager;
  odt_manager.initializeLevel(0, ba, dm);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));
  const pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  const auto support_data = odt_manager.collectSupportData(
    line_geom, geom, state_valid, state_same_level, state_host_filled);

  ASSERT_TRUE(support_data.allSameLevelAccepted());
  ASSERT_EQ(
    static_cast<int>(support_data.ordered_samples.size()),
    line_geom.supportCellCount());
  for (const auto& sample : support_data.ordered_samples) {
    ASSERT_TRUE(sample.has_value);
    ASSERT_EQ(sample.value.rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
    const amrex::Real x = static_cast<amrex::Real>(sample.iv[0]);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      const amrex::Real expected =
        0.01 * static_cast<amrex::Real>(n + 1) * (1.0 + x);
      EXPECT_NEAR(sample.value.rhoY[static_cast<std::size_t>(n)], expected, tol);
    }
  }
}

TEST(
  ODTLESLocalEngine,
  ODTLineInitializationFromLocalSupportIncludesSpeciesConservativeState)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  std::vector<ConservativeCell> support(
    static_cast<std::size_t>(line_geom.supportCellCount()));
  const int i_m = supportIndexForOffset(line_geom, -1);
  const int i_0 = supportIndexForOffset(line_geom, 0);
  const int i_p = supportIndexForOffset(line_geom, 1);
  ASSERT_GE(i_m, 0);
  ASSERT_GE(i_0, 0);
  ASSERT_GE(i_p, 0);

  ConservativeCell c_m_in = makeCell(1.0, 0.40, -0.15, 0.30, 4.2);
  ConservativeCell c_0_in = makeCell(1.5, 0.60, -0.05, 0.45, 4.8);
  ConservativeCell c_p_in = makeCell(2.2, 0.95, 0.10, 0.62, 5.9);
  setSpeciesSequence(c_m_in, 0.01);
  setSpeciesSequence(c_0_in, 0.02);
  setSpeciesSequence(c_p_in, 0.03);
  setSupportCellsForOffset(line_geom, support, -1, c_m_in);
  setSupportCellsForOffset(line_geom, support, 0, c_0_in);
  setSupportCellsForOffset(line_geom, support, 1, c_p_in);

  pelec::odtles::ODTLineState line_state;
  pelec::odtles::ODTReconcile::initializeLineStateFromLESSupportAverages(
    line_geom, support, line_state);

  ASSERT_TRUE(line_state.initialized());
  ASSERT_TRUE(line_state.valid());
  const auto& c_m = line_state.cell(i_m);
  const auto& c_0 = line_state.cell(i_0);
  const auto& c_p = line_state.cell(i_p);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_NEAR(
      c_m.rhoY[static_cast<std::size_t>(n)],
      support[static_cast<std::size_t>(i_m)].rhoY[static_cast<std::size_t>(n)],
      tol);
    EXPECT_NEAR(
      c_0.rhoY[static_cast<std::size_t>(n)],
      support[static_cast<std::size_t>(i_0)].rhoY[static_cast<std::size_t>(n)],
      tol);
    EXPECT_NEAR(
      c_p.rhoY[static_cast<std::size_t>(n)],
      support[static_cast<std::size_t>(i_p)].rhoY[static_cast<std::size_t>(n)],
      tol);
  }

  for (int i = 0; i < line_state.numCells(); ++i) {
    ASSERT_EQ(
      line_state.cell(i).rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
    for (int n = 0; n < NUM_SPECIES; ++n) {
      EXPECT_GE(line_state.cell(i).rhoY[static_cast<std::size_t>(n)], 0.0);
    }
  }
}

TEST(ODTLESLocalEngine, ODTLineInitializationFromLocalSupportIsNonConstant)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  std::vector<ConservativeCell> support(
    static_cast<std::size_t>(line_geom.supportCellCount()));
  const int i_m = supportIndexForOffset(line_geom, -1);
  const int i_0 = supportIndexForOffset(line_geom, 0);
  const int i_p = supportIndexForOffset(line_geom, 1);
  ASSERT_GE(i_m, 0);
  ASSERT_GE(i_0, 0);
  ASSERT_GE(i_p, 0);

  setSupportCellsForOffset(
    line_geom, support, -1, makeCell(1.0, 0.40, -0.15, 0.30, 4.2));
  setSupportCellsForOffset(
    line_geom, support, 0, makeCell(1.5, 0.60, -0.05, 0.45, 4.8));
  setSupportCellsForOffset(
    line_geom, support, 1, makeCell(2.2, 0.95, 0.10, 0.62, 5.9));

  pelec::odtles::ODTLineState line_state;
  pelec::odtles::ODTReconcile::initializeLineStateFromLESSupportAverages(
    line_geom, support, line_state);

  ASSERT_TRUE(line_state.initialized());
  ASSERT_TRUE(line_state.valid());
  EXPECT_EQ(line_state.numCells(), line_geom.supportCellCount());

  const auto& c_m = line_state.cell(i_m);
  const auto& c_0 = line_state.cell(i_0);
  const auto& c_p = line_state.cell(i_p);
  EXPECT_GT(std::abs(c_m.rho - c_0.rho), tol);
  EXPECT_GT(std::abs(c_p.rho - c_0.rho), tol);
  EXPECT_GT(std::abs(c_m.rhou - c_0.rhou), tol);
  EXPECT_GT(std::abs(c_p.rhou - c_0.rhou), tol);
}

TEST(ODTLESLocalEngine, ODTLineInitializationPreservesOwnerMeanExactly)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  std::vector<ConservativeCell> support(
    static_cast<std::size_t>(line_geom.supportCellCount()));
  const int i_m = supportIndexForOffset(line_geom, -1);
  const int i_0 = supportIndexForOffset(line_geom, 0);
  const int i_p = supportIndexForOffset(line_geom, 1);
  ASSERT_GE(i_m, 0);
  ASSERT_GE(i_0, 0);
  ASSERT_GE(i_p, 0);

  setSupportCellsForOffset(
    line_geom, support, -1, makeCell(1.1, 0.30, -0.08, 0.20, 3.8));
  setSupportCellsForOffset(
    line_geom, support, 1, makeCell(2.4, 0.90, 0.11, 0.50, 6.1));

  const auto owner_idxs = supportIndicesForOffset(line_geom, 0);
  ASSERT_EQ(
    static_cast<int>(owner_idxs.size()),
    pelec::odtles::ODTLineGeometry::DefaultSubsegmentsPerHostCell);
  support[static_cast<std::size_t>(owner_idxs[0])] =
    makeCell(1.4, 0.45, 0.00, 0.28, 4.5);
  support[static_cast<std::size_t>(owner_idxs[1])] =
    makeCell(1.7, 0.55, 0.02, 0.33, 4.9);
  support[static_cast<std::size_t>(owner_idxs[2])] =
    makeCell(2.0, 0.65, 0.04, 0.38, 5.3);

  ConservativeCell owner_target{};
  owner_target.rho = 0.0;
  owner_target.rhou = 0.0;
  owner_target.rhov = 0.0;
  owner_target.rhow = 0.0;
  owner_target.rhoE = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    owner_target.rhoY[static_cast<std::size_t>(n)] = 0.0;
  }
  for (int i : owner_idxs) {
    const auto& c = support[static_cast<std::size_t>(i)];
    owner_target.rho += c.rho / static_cast<amrex::Real>(owner_idxs.size());
    owner_target.rhou += c.rhou / static_cast<amrex::Real>(owner_idxs.size());
    owner_target.rhov += c.rhov / static_cast<amrex::Real>(owner_idxs.size());
    owner_target.rhow += c.rhow / static_cast<amrex::Real>(owner_idxs.size());
    owner_target.rhoE += c.rhoE / static_cast<amrex::Real>(owner_idxs.size());
    for (int n = 0; n < NUM_SPECIES; ++n) {
      owner_target.rhoY[static_cast<std::size_t>(n)] +=
        c.rhoY[static_cast<std::size_t>(n)] /
        static_cast<amrex::Real>(owner_idxs.size());
    }
  }

  pelec::odtles::ODTLineState line_state;
  pelec::odtles::ODTReconcile::initializeLineStateFromLESSupportAverages(
    line_geom, support, line_state);

  const auto owner_mean = ownerIntervalMean(line_geom, line_state);
  EXPECT_NEAR(owner_mean.rho, owner_target.rho, tol);
  EXPECT_NEAR(owner_mean.rhou, owner_target.rhou, tol);
  EXPECT_NEAR(owner_mean.rhov, owner_target.rhov, tol);
  EXPECT_NEAR(owner_mean.rhow, owner_target.rhow, tol);
  EXPECT_NEAR(owner_mean.rhoE, owner_target.rhoE, tol);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_NEAR(
      owner_mean.rhoY[static_cast<std::size_t>(n)],
      owner_target.rhoY[static_cast<std::size_t>(n)], tol);
  }
}

TEST(
  ODTLESLocalEngine,
  ODTLineInitializationAndReconcilePreserveOwnerMeanWithFiveSubsegments)
{
  constexpr amrex::Real tol = 1.0e-12;
  constexpr int n_sub = 5;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0, n_sub);
  std::vector<ConservativeCell> support(
    static_cast<std::size_t>(line_geom.supportCellCount()));
  setSupportCellsForOffset(
    line_geom, support, -1, makeCell(1.1, 0.30, -0.08, 0.20, 3.8));
  setSupportCellsForOffset(
    line_geom, support, 1, makeCell(2.4, 0.90, 0.11, 0.50, 6.1));

  const auto owner_idxs = supportIndicesForOffset(line_geom, 0);
  ASSERT_EQ(static_cast<int>(owner_idxs.size()), n_sub);

  ConservativeCell owner_target{};
  owner_target.rho = 0.0;
  owner_target.rhou = 0.0;
  owner_target.rhov = 0.0;
  owner_target.rhow = 0.0;
  owner_target.rhoE = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    owner_target.rhoY[static_cast<std::size_t>(n)] = 0.0;
  }

  for (int k = 0; k < n_sub; ++k) {
    ConservativeCell ck = makeCell(
      1.3 + 0.2 * static_cast<amrex::Real>(k),
      0.40 + 0.05 * static_cast<amrex::Real>(k),
      -0.02 + 0.01 * static_cast<amrex::Real>(k),
      0.24 + 0.03 * static_cast<amrex::Real>(k),
      4.3 + 0.4 * static_cast<amrex::Real>(k));
    support[static_cast<std::size_t>(owner_idxs[static_cast<std::size_t>(k)])] = ck;
    owner_target.rho += ck.rho / static_cast<amrex::Real>(n_sub);
    owner_target.rhou += ck.rhou / static_cast<amrex::Real>(n_sub);
    owner_target.rhov += ck.rhov / static_cast<amrex::Real>(n_sub);
    owner_target.rhow += ck.rhow / static_cast<amrex::Real>(n_sub);
    owner_target.rhoE += ck.rhoE / static_cast<amrex::Real>(n_sub);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      owner_target.rhoY[static_cast<std::size_t>(n)] +=
        ck.rhoY[static_cast<std::size_t>(n)] / static_cast<amrex::Real>(n_sub);
    }
  }

  pelec::odtles::ODTLineState line_state;
  pelec::odtles::ODTReconcile::initializeLineStateFromLESSupportAverages(
    line_geom, support, line_state);

  const auto owner_mean_init = ownerIntervalMean(line_geom, line_state);
  EXPECT_NEAR(owner_mean_init.rho, owner_target.rho, tol);
  EXPECT_NEAR(owner_mean_init.rhou, owner_target.rhou, tol);
  EXPECT_NEAR(owner_mean_init.rhov, owner_target.rhov, tol);
  EXPECT_NEAR(owner_mean_init.rhow, owner_target.rhow, tol);
  EXPECT_NEAR(owner_mean_init.rhoE, owner_target.rhoE, tol);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_NEAR(
      owner_mean_init.rhoY[static_cast<std::size_t>(n)],
      owner_target.rhoY[static_cast<std::size_t>(n)], tol);
  }

  ConservativeCell reconcile_target = makeCell(1.9, 0.52, -0.01, 0.33, 5.1);
  setSpeciesSequence(reconcile_target, 0.025, 0.002);
  pelec::odtles::ODTReconcile::reconcileExistingLineStateToOwnerAverage(
    line_geom, reconcile_target, line_state);

  const auto owner_mean_reconcile = ownerIntervalMean(line_geom, line_state);
  EXPECT_NEAR(owner_mean_reconcile.rho, reconcile_target.rho, tol);
  EXPECT_NEAR(owner_mean_reconcile.rhou, reconcile_target.rhou, tol);
  EXPECT_NEAR(owner_mean_reconcile.rhov, reconcile_target.rhov, tol);
  EXPECT_NEAR(owner_mean_reconcile.rhow, reconcile_target.rhow, tol);
  EXPECT_NEAR(owner_mean_reconcile.rhoE, reconcile_target.rhoE, tol);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_NEAR(
      owner_mean_reconcile.rhoY[static_cast<std::size_t>(n)],
      reconcile_target.rhoY[static_cast<std::size_t>(n)], tol);
  }
}

TEST(
  ODTLESLocalEngine,
  RuntimeODTLinePreparationDistinguishesInitializeAndReconcile)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  amrex::BoxArray ba(domain);
  amrex::DistributionMapping dm(ba);

  constexpr int support_ng = pelec::odtles::ODTLineGeometry::supportGhostCells();
  amrex::MultiFab state_valid(ba, dm, NVAR, 0);
  amrex::MultiFab state_same_level(ba, dm, NVAR, support_ng);
  amrex::MultiFab state_host_filled(ba, dm, NVAR, support_ng);
  state_valid.setVal(0.0);
  state_same_level.setVal(0.0);
  state_host_filled.setVal(0.0);

  for (amrex::MFIter mfi(state_valid, false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    auto const s = state_valid.array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      const amrex::Real x = static_cast<amrex::Real>(iv[0]);
      s(iv, URHO) = 1.2 + 0.2 * x;
      s(iv, UMX) = 0.24 + 0.04 * x;
      s(iv, UMY) = -0.12 + 0.02 * x;
      s(iv, UMZ) = 0.06 + 0.01 * x;
      s(iv, UEDEN) = 4.0 + 0.3 * x;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s(iv, UFS + n) =
          0.01 * static_cast<amrex::Real>(n + 1) * (2.0 + x);
      }
    }
  }
  amrex::MultiFab::Copy(state_same_level, state_valid, 0, 0, NVAR, 0);
  state_same_level.FillBoundary(geom.periodicity());
  amrex::MultiFab::Copy(
    state_host_filled, state_same_level, 0, 0, NVAR, state_same_level.nGrow());

  pelec::odtles::ODTManager odt_manager;
  odt_manager.initializeLevel(0, ba, dm);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::RuntimeDepositionStats stats{};
  const auto prep0 = pelec::odtles::prepareRuntimeODTLineFromLESSupport(
    0, owner, 0, geom, odt_manager, state_valid, state_same_level,
    state_host_filled, stats);
  ASSERT_TRUE(prep0.accepted);
  EXPECT_TRUE(prep0.initialized);
  EXPECT_FALSE(prep0.reconciled);

  auto* entry0 = odt_manager.findLineEntry(0, owner, 0);
  ASSERT_NE(entry0, nullptr);
  const int owner_local = entry0->geometry.ownerLocalOrdinal();
  const int left_local = supportIndexForOffset(entry0->geometry, -1);
  ASSERT_GE(left_local, 0);
  ASSERT_GE(owner_local, 0);
  const amrex::Real old_left_rho = entry0->state.cell(left_local).rho;
  const amrex::Real old_owner_rho = entry0->state.cell(owner_local).rho;
  EXPECT_GT(std::abs(old_left_rho - old_owner_rho), tol);

  ConservativeCell target = makeCell(1.85, 0.37, -0.01, 0.095, 5.4);
  setSpeciesSequence(target, 0.015, 0.001);
  setConservativeAt(state_valid, owner, target);
  setConservativeAt(state_same_level, owner, target);
  setConservativeAt(state_host_filled, owner, target);

  const auto prep1 = pelec::odtles::prepareRuntimeODTLineFromLESSupport(
    0, owner, 0, geom, odt_manager, state_valid, state_same_level,
    state_host_filled, stats);
  ASSERT_TRUE(prep1.accepted);
  EXPECT_FALSE(prep1.initialized);
  EXPECT_TRUE(prep1.reconciled);

  auto* entry1 = odt_manager.findLineEntry(0, owner, 0);
  ASSERT_NE(entry1, nullptr);
  const auto owner_mean = ownerIntervalMean(entry1->geometry, entry1->state);
  EXPECT_NEAR(owner_mean.rho, target.rho, tol);
  EXPECT_NEAR(owner_mean.rhou, target.rhou, tol);
  EXPECT_NEAR(owner_mean.rhov, target.rhov, tol);
  EXPECT_NEAR(owner_mean.rhow, target.rhow, tol);
  EXPECT_NEAR(owner_mean.rhoE, target.rhoE, tol);
  ASSERT_EQ(owner_mean.rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_NEAR(
      owner_mean.rhoY[static_cast<std::size_t>(n)],
      target.rhoY[static_cast<std::size_t>(n)], tol);
  }

  // Reconcile path preserves residual structure away from owner instead of
  // reinitializing all support cells to owner target.
  EXPECT_GT(std::abs(entry1->state.cell(left_local).rho - target.rho), tol);
}

TEST(ODTLESLocalEngine, ODTLineReconcilePreservesResidualWithAlphaReduction)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  pelec::odtles::ODTLineState line_state;
  line_state.initialize(line_geom);

  const int i_m = supportIndexForOffset(line_geom, -1);
  const int i_0 = supportIndexForOffset(line_geom, 0);
  const int i_p = supportIndexForOffset(line_geom, 1);
  ASSERT_GE(i_m, 0);
  ASSERT_GE(i_0, 0);
  ASSERT_GE(i_p, 0);

  ConservativeCell old_owner = makeCell(1.0, 0.40, 0.20, -0.10, 4.0);
  ConservativeCell old_left = makeCell(0.2, 0.08, 0.04, -0.02, 0.8);
  ConservativeCell old_right = makeCell(1.6, 0.64, 0.32, -0.16, 6.4);
  setSpeciesSequence(old_owner, 0.12);
  setSpeciesSequence(old_left, 0.02);
  setSpeciesSequence(old_right, 0.24);
  for (int i : supportIndicesForOffset(line_geom, -1)) {
    line_state.setCell(i, old_left);
  }
  for (int i : supportIndicesForOffset(line_geom, 0)) {
    line_state.setCell(i, old_owner);
  }
  for (int i : supportIndicesForOffset(line_geom, 1)) {
    line_state.setCell(i, old_right);
  }
  line_state.setValid(true);

  ConservativeCell target = makeCell(0.3, 0.12, 0.06, -0.03, 1.2);
  setSpeciesSequence(target, 0.03);
  pelec::odtles::ODTReconcile::ReconcileControls ctrl{};
  ctrl.rho_floor = 0.15;
  ctrl.e_floor = 1.0e-4;
  ctrl.alpha_min = 1.0e-6;
  ctrl.max_alpha_reductions = 60;

  const auto old_owner_v = asArray(old_owner);
  const auto old_left_v = asArray(old_left);
  std::array<amrex::Real, 5> old_left_res{};
  for (int n = 0; n < 5; ++n) {
    old_left_res[static_cast<std::size_t>(n)] =
      old_left_v[static_cast<std::size_t>(n)] -
      old_owner_v[static_cast<std::size_t>(n)];
  }

  pelec::odtles::ODTReconcile::reconcileExistingLineStateToOwnerAverage(
    line_geom, target, line_state, ctrl);

  const auto new_owner_mean = ownerIntervalMean(line_geom, line_state);
  EXPECT_NEAR(new_owner_mean.rho, target.rho, tol);
  EXPECT_NEAR(new_owner_mean.rhou, target.rhou, tol);
  EXPECT_NEAR(new_owner_mean.rhov, target.rhov, tol);
  EXPECT_NEAR(new_owner_mean.rhow, target.rhow, tol);
  EXPECT_NEAR(new_owner_mean.rhoE, target.rhoE, tol);
  ASSERT_EQ(new_owner_mean.rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    EXPECT_NEAR(
      new_owner_mean.rhoY[static_cast<std::size_t>(n)],
      target.rhoY[static_cast<std::size_t>(n)], tol);
  }

  const auto target_v = asArray(target);
  const auto new_left_v = asArray(line_state.cell(i_m));
  std::array<amrex::Real, 5> new_left_res{};
  for (int n = 0; n < 5; ++n) {
    new_left_res[static_cast<std::size_t>(n)] =
      new_left_v[static_cast<std::size_t>(n)] -
      target_v[static_cast<std::size_t>(n)];
  }

  const amrex::Real alpha =
    new_left_res[0] / old_left_res[0];
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 1.0);
  EXPECT_GT(std::abs(new_left_res[0]), tol);
  EXPECT_LT(std::abs(new_left_res[0]), std::abs(old_left_res[0]));
  for (int n = 1; n < 5; ++n) {
    EXPECT_NEAR(
      new_left_res[static_cast<std::size_t>(n)],
      alpha * old_left_res[static_cast<std::size_t>(n)], 5.0e-12);
  }

  for (int i = 0; i < line_state.numCells(); ++i) {
    const auto c = line_state.cell(i);
    EXPECT_GT(c.rho, ctrl.rho_floor);
    EXPECT_GT(specificInternalEnergy(c), ctrl.e_floor);
    ASSERT_EQ(c.rhoY.size(), static_cast<std::size_t>(NUM_SPECIES));
    for (int n = 0; n < NUM_SPECIES; ++n) {
      EXPECT_GE(c.rhoY[static_cast<std::size_t>(n)], 0.0);
    }
  }
}

TEST(ODTLESLocalEngine, RuntimeODTLinePreparationRejectsUnsupportedProvenance)
{
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);

  const amrex::Box patch(
    amrex::IntVect(AMREX_D_DECL(1, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(3, 0, 0)));
  amrex::BoxArray ba(patch);
  amrex::DistributionMapping dm(ba);

  amrex::MultiFab state_valid(ba, dm, NVAR, 0);
  amrex::MultiFab state_same_level(ba, dm, NVAR, 0);
  amrex::MultiFab state_host_filled(ba, dm, NVAR, 1);
  state_valid.setVal(0.0);
  state_same_level.setVal(0.0);
  state_host_filled.setVal(7.0);
  amrex::MultiFab::Copy(state_same_level, state_valid, 0, 0, NVAR, 0);
  amrex::MultiFab::Copy(state_host_filled, state_valid, 0, 0, NVAR, 0);

  pelec::odtles::ODTManager odt_manager;
  odt_manager.initializeLevel(0, ba, dm);
  pelec::odtles::RuntimeDepositionStats stats{};
  const amrex::IntVect owner(AMREX_D_DECL(1, 0, 0));

  const auto prep = pelec::odtles::prepareRuntimeODTLineFromLESSupport(
    0, owner, 0, geom, odt_manager, state_valid, state_same_level,
    state_host_filled, stats);

  EXPECT_FALSE(prep.accepted);
  EXPECT_FALSE(prep.initialized);
  EXPECT_FALSE(prep.reconciled);
  EXPECT_EQ(prep.line_entry, nullptr);
  EXPECT_EQ(stats.rejected_amr_entries, 1);
  EXPECT_EQ(stats.rejected_boundary_entries, 0);
  EXPECT_EQ(stats.rejected_invalid_entries, 0);
  EXPECT_EQ(stats.rejected_mixed_entries, 0);
  EXPECT_FALSE(odt_manager.hasLineEntry(0, owner, 0));
}

TEST(
  ODTLESLocalEngine,
  RuntimeODTLinePreparationUsesConfiguredSubsegmentsPerHostCell)
{
  constexpr amrex::Real tol = 1.0e-12;
  constexpr int n_sub = 5;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  amrex::BoxArray ba(domain);
  amrex::DistributionMapping dm(ba);

  constexpr int support_ng = pelec::odtles::ODTLineGeometry::supportGhostCells();
  amrex::MultiFab state_valid(ba, dm, NVAR, 0);
  amrex::MultiFab state_same_level(ba, dm, NVAR, support_ng);
  amrex::MultiFab state_host_filled(ba, dm, NVAR, support_ng);
  state_valid.setVal(0.0);
  state_same_level.setVal(0.0);
  state_host_filled.setVal(0.0);

  for (amrex::MFIter mfi(state_valid, false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    auto const s = state_valid.array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      const amrex::Real x = static_cast<amrex::Real>(iv[0]);
      s(iv, URHO) = 1.2 + 0.1 * x;
      s(iv, UMX) = 0.5 + 0.05 * x;
      s(iv, UMY) = 0.0;
      s(iv, UMZ) = 0.0;
      s(iv, UEDEN) = 4.0 + 0.2 * x;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s(iv, UFS + n) = (1.2 + 0.1 * x) / static_cast<amrex::Real>(NUM_SPECIES);
      }
    }
  }
  amrex::MultiFab::Copy(state_same_level, state_valid, 0, 0, NVAR, 0);
  state_same_level.FillBoundary(geom.periodicity());
  amrex::MultiFab::Copy(
    state_host_filled, state_same_level, 0, 0, NVAR, state_same_level.nGrow());

  pelec::odtles::ODTManager odt_manager;
  pelec::odtles::ODTParams params{};
  params.subsegments_per_host_cell = n_sub;
  odt_manager.setParams(params);
  odt_manager.initializeLevel(0, ba, dm);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::RuntimeDepositionStats stats{};
  const auto prep = pelec::odtles::prepareRuntimeODTLineFromLESSupport(
    0, owner, 0, geom, odt_manager, state_valid, state_same_level,
    state_host_filled, stats);
  ASSERT_TRUE(prep.accepted);
  ASSERT_NE(prep.line_entry, nullptr);

  const auto& line_geom = prep.line_entry->geometry;
  EXPECT_EQ(line_geom.subsegmentsPerHostCell(), n_sub);
  EXPECT_EQ(
    line_geom.supportCellCount(),
    pelec::odtles::ODTLineGeometry::HostSupportCells * n_sub);
  EXPECT_EQ(static_cast<int>(supportIndicesForOffset(line_geom, 0).size()), n_sub);

  const auto col = pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
    line_geom, prep.line_entry->state);
  ASSERT_TRUE(col.valid);
  EXPECT_EQ(col.contributing_segments, n_sub);
  EXPECT_NEAR(col.overlap_weight_sum, 1.0, tol);
}

TEST(ODTLESLocalEngine, DirectionalMomentumContributionLayout)
{
  pelec::odtles::ODTMomentExtractor::DirectionalColumn col{};
  col.dir = 2;
  col.valid = true;
  col.tau_ij = {1.5, -2.0, 0.25};

  const auto q =
    pelec::odtles::ODTMomentExtractor::buildDirectionalMomentumContribution(col);

  EXPECT_TRUE(q.valid);
  EXPECT_EQ(q.dir, 2);
  EXPECT_DOUBLE_EQ(q.q[0], -1.5);
  EXPECT_DOUBLE_EQ(q.q[1], 2.0);
  EXPECT_DOUBLE_EQ(q.q[2], -0.25);
  EXPECT_DOUBLE_EQ(q.q[3], 0.0);
}

TEST(ODTLESLocalEngine, ConservativeMomentumDepositionValidation)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box cbox(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(3, 0, 0)));
  const amrex::Box domain = cbox;

  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> odt_dir_cc{
    AMREX_D_DECL(
      amrex::FArrayBox(cbox, 4),
      amrex::FArrayBox(cbox, 4),
      amrex::FArrayBox(cbox, 4))};
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> area{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 0), 1),
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 1), 1),
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 2), 1))};
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 0), NVAR),
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 1), NVAR),
      amrex::FArrayBox(amrex::surroundingNodes(cbox, 2), NVAR))};
  amrex::FArrayBox volume(cbox, 1);
  amrex::FArrayBox lterm(cbox, NVAR);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    odt_dir_cc[dir].setVal(0.0);
    area[dir].setVal(1.0);
    flux_ec[dir].setVal(0.0);
  }
  volume.setVal(1.0);
  lterm.setVal(0.0);

  // Case 1: uniform centered ODT directional momentum contribution should
  // produce zero divergence after conservative flux differencing.
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    auto const q = odt_dir_cc[dir].array();
    for (amrex::IntVect iv = cbox.smallEnd(); iv <= cbox.bigEnd(); cbox.next(iv)) {
      q(iv, 0) = 1.0;
      q(iv, 1) = -2.0;
      q(iv, 2) = 0.5;
      q(iv, 3) = 0.0;
    }
  }
  buildOdtDirectionalFluxes(cbox, domain, odt_dir_cc, area, flux_ec);
  computeDivFromFluxes(cbox, flux_ec, volume, lterm);

  for (amrex::IntVect iv = cbox.smallEnd(); iv <= cbox.bigEnd(); cbox.next(iv)) {
    for (int n = 0; n < NVAR; ++n) {
      EXPECT_LE(std::abs(lterm(iv, n)), tol);
    }
  }

  // Case 2: controlled nonuniform x-direction contribution (j=0 only) with
  // known sign pattern in momentum divergence. Keep UEDEN inactive and all
  // non-momentum components untouched.
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    odt_dir_cc[dir].setVal(0.0);
    flux_ec[dir].setVal(0.0);
  }
  lterm.setVal(0.0);

  auto const qx = odt_dir_cc[0].array();
  for (amrex::IntVect iv = cbox.smallEnd(); iv <= cbox.bigEnd(); cbox.next(iv)) {
    const amrex::Real x = static_cast<amrex::Real>(iv[0]);
    qx(iv, 0) = x;       // slope +1 -> interior divergence -1
    qx(iv, 1) = -x;      // slope -1 -> interior divergence +1
    qx(iv, 2) = 2.0 * x; // slope +2 -> interior divergence -2
    qx(iv, 3) = 0.0;     // UEDEN inactive
  }

  buildOdtDirectionalFluxes(cbox, domain, odt_dir_cc, area, flux_ec);
  computeDivFromFluxes(cbox, flux_ec, volume, lterm);

  const amrex::IntVect iv1(AMREX_D_DECL(1, 0, 0));
  const amrex::IntVect iv2(AMREX_D_DECL(2, 0, 0));
  const amrex::IntVect iv0(AMREX_D_DECL(0, 0, 0));
  const amrex::IntVect iv3(AMREX_D_DECL(3, 0, 0));

  // Boundary-sensitive one-sided edge behavior.
  EXPECT_NEAR(lterm(iv0, UMX), -0.5, tol);
  EXPECT_NEAR(lterm(iv1, UMX), -1.0, tol);
  EXPECT_NEAR(lterm(iv2, UMX), -1.0, tol);
  EXPECT_NEAR(lterm(iv3, UMX), -0.5, tol);
  EXPECT_NEAR(lterm(iv0, UMY), 0.5, tol);
  EXPECT_NEAR(lterm(iv1, UMY), 1.0, tol);
  EXPECT_NEAR(lterm(iv2, UMY), 1.0, tol);
  EXPECT_NEAR(lterm(iv3, UMY), 0.5, tol);
  EXPECT_NEAR(lterm(iv0, UMZ), -1.0, tol);
  EXPECT_NEAR(lterm(iv1, UMZ), -2.0, tol);
  EXPECT_NEAR(lterm(iv2, UMZ), -2.0, tol);
  EXPECT_NEAR(lterm(iv3, UMZ), -1.0, tol);

  // Global conservative budget over the box:
  // sum(V*D) = -(F_right - F_left) per component.
  const auto fx = flux_ec[0].const_array();
  const amrex::Real rhs_umx =
    -(fx(4, 0, 0, UMX) - fx(0, 0, 0, UMX));
  const amrex::Real rhs_umy =
    -(fx(4, 0, 0, UMY) - fx(0, 0, 0, UMY));
  const amrex::Real rhs_umz =
    -(fx(4, 0, 0, UMZ) - fx(0, 0, 0, UMZ));
  amrex::Real sum_umx = 0.0;
  amrex::Real sum_umy = 0.0;
  amrex::Real sum_umz = 0.0;
  for (amrex::IntVect iv = cbox.smallEnd(); iv <= cbox.bigEnd(); cbox.next(iv)) {
    sum_umx += lterm(iv, UMX);
    sum_umy += lterm(iv, UMY);
    sum_umz += lterm(iv, UMZ);
  }
  EXPECT_NEAR(sum_umx, rhs_umx, tol);
  EXPECT_NEAR(sum_umy, rhs_umy, tol);
  EXPECT_NEAR(sum_umz, rhs_umz, tol);

  for (amrex::IntVect iv = cbox.smallEnd(); iv <= cbox.bigEnd(); cbox.next(iv)) {
    EXPECT_LE(std::abs(lterm(iv, UEDEN)), tol);
    EXPECT_LE(std::abs(lterm(iv, URHO)), tol);
    for (int n = UFS; n < UFS + NUM_SPECIES; ++n) {
      EXPECT_LE(std::abs(lterm(iv, n)), tol);
    }
#if NUM_AUX > 0
    for (int n = UFX; n < UFX + NUM_AUX; ++n) {
      EXPECT_LE(std::abs(lterm(iv, n)), tol);
    }
#endif
#if NUM_ADV > 0
    for (int n = UFA; n < UFA + NUM_ADV; ++n) {
      EXPECT_LE(std::abs(lterm(iv, n)), tol);
    }
#endif
  }
}

TEST(ODTLESLocalEngine, ConservativeMomentumDepositionPartitionConsistency)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(7, 0, 0)));

  amrex::BoxArray ba_single(domain);
  amrex::DistributionMapping dm_single(ba_single);
  amrex::MultiFab lterm_single;
  computeLtermFromCenteredContributionMF(ba_single, dm_single, domain, lterm_single);

  amrex::BoxArray ba_split(domain);
  ba_split.maxSize(4);
  amrex::DistributionMapping dm_split(ba_split);
  amrex::MultiFab lterm_split;
  computeLtermFromCenteredContributionMF(ba_split, dm_split, domain, lterm_split);

  for (amrex::IntVect iv = domain.smallEnd(); iv <= domain.bigEnd();
       domain.next(iv)) {
    EXPECT_NEAR(valueAt(lterm_single, iv, UMX), valueAt(lterm_split, iv, UMX), tol);
    EXPECT_NEAR(valueAt(lterm_single, iv, UMY), valueAt(lterm_split, iv, UMY), tol);
    EXPECT_NEAR(valueAt(lterm_single, iv, UMZ), valueAt(lterm_split, iv, UMZ), tol);
    EXPECT_NEAR(valueAt(lterm_single, iv, UEDEN), valueAt(lterm_split, iv, UEDEN), tol);
  }
}

TEST(
  ODTLESLocalEngine,
  RuntimeODTMomentumDepositionDefaultConfigMatchesExplicitCanonicalControls)
{
  constexpr amrex::Real tol = 1.0e-12;
  amrex::MultiFab lterm_default;
  amrex::MultiFab lterm_explicit;

  const pelec::odtles::ODTParams params_default{};
  const RuntimeSignature sig_default =
    computeRuntimeODTLESTermSignature(lterm_default, params_default);

  pelec::odtles::ODTParams params_explicit{};
  params_explicit.max_local_substeps = 1;
  const RuntimeSignature sig_explicit =
    computeRuntimeODTLESTermSignature(lterm_explicit, params_explicit);

  EXPECT_NEAR(sig_default.sum_umx, sig_explicit.sum_umx, tol);
  EXPECT_NEAR(sig_default.sum_umy, sig_explicit.sum_umy, tol);
  EXPECT_NEAR(sig_default.sum_umz, sig_explicit.sum_umz, tol);
  EXPECT_NEAR(sig_default.sum_ueden, sig_explicit.sum_ueden, tol);
  EXPECT_NEAR(sig_default.l1_umx, sig_explicit.l1_umx, tol);
  EXPECT_NEAR(sig_default.l1_umy, sig_explicit.l1_umy, tol);
  EXPECT_NEAR(sig_default.l1_umz, sig_explicit.l1_umz, tol);
  EXPECT_NEAR(sig_default.l1_ueden, sig_explicit.l1_ueden, tol);
  EXPECT_NEAR(sig_default.l1_others, sig_explicit.l1_others, tol);
}

TEST(
  ODTLESLocalEngine,
  RuntimeMolecularMomentumPathUsesProductionFormAndAllowsDebugComparisonOverride)
{
  constexpr amrex::Real tol_quiet = 1.0e-12;
  constexpr amrex::Real tol_active = 1.0e-10;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(4, 0, 0)));
  const amrex::RealBox rb(
    {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(5.0, 1.0, 1.0)});
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(2, 0, 0));

  pelec::odtles::ODTLineGeometry line_geom(geom, 0, owner, 0);
  pelec::odtles::ODTLineState base_state;
  base_state.initialize(line_geom);
  auto eos = pele::physics::PhysicsType::eos();

  const amrex::Real u_const = 0.65;
  amrex::Real Y_ref[NUM_SPECIES] = {0.0};
  amrex::Real ysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y_ref[n] = static_cast<amrex::Real>(n + 1);
    ysum += Y_ref[n];
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y_ref[n] /= ysum;
  }
  for (int i = 0; i < base_state.numCells(); ++i) {
    ConservativeCell c{};
    const amrex::Real ii = static_cast<amrex::Real>(i);
    c.rho = 1.0 + 0.3 * ii * ii;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      c.rhoY[static_cast<std::size_t>(n)] = c.rho * Y_ref[n];
    }
    amrex::Real e_ref = 0.0;
    eos.RTY2E(c.rho, 410.0, Y_ref, e_ref);
    c.rhou = c.rho * u_const;
    c.rhov = 0.0;
    c.rhow = 0.0;
    c.rhoE = c.rho * (e_ref + 0.5 * u_const * u_const);
    base_state.setCell(i, c);
  }
  base_state.setValid(true);

  pelec::odtles::ODTParams params{};
  pelec::odtles::ODTStepper::Controls step_ctrl{};
  pelec::odtles::configureRuntimeODTStepperControls(params, step_ctrl);
  EXPECT_EQ(
    step_ctrl.diffusion_controls.momentum_molecular_form,
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::VelocityGradientMuFlux);
  step_ctrl.diffusion_controls.molecular_viscosity_floor = 1.0;

  pelec::odtles::ODTLineState production_state = base_state;
  pelec::odtles::ODTDiffusion::applyImplicitUniform(
    line_geom, production_state, 1.0e-2, step_ctrl.diffusion_controls);

  auto debug_diff_ctrl = step_ctrl.diffusion_controls;
  debug_diff_ctrl.momentum_molecular_form =
    pelec::odtles::ODTDiffusion::MomentumMolecularForm::ConservativeVariableNuOnRhoU;
  pelec::odtles::ODTLineState debug_state = base_state;
  pelec::odtles::ODTDiffusion::applyImplicitUniform(
    line_geom, debug_state, 1.0e-2, debug_diff_ctrl);

  const int comp_rhou = static_cast<int>(pelec::odtles::ODTDiffusion::Component::RhoU);
  const amrex::Real change_production = maxAbsMomentumComponentChange(
    base_state, production_state, comp_rhou);
  const amrex::Real change_debug = maxAbsMomentumComponentChange(
    base_state, debug_state, comp_rhou);

  EXPECT_LE(change_production, tol_quiet);
  EXPECT_GT(change_debug, tol_active);
}

TEST(ODTLESLocalEngine, RuntimeODTMomentumDepositionMPIBudget)
{
  constexpr amrex::Real tol = 1.0e-12;
  amrex::MultiFab lterm_runtime;
  const RuntimeSignature sig = computeRuntimeODTLESTermSignature(lterm_runtime);

  // Ensure the real runtime path was exercised.
  EXPECT_GT(sig.stepped_entries, 0);
  EXPECT_GT(sig.moment_columns_built, 0);

  // Explicit variable targeting.
  EXPECT_LE(std::abs(sig.sum_ueden), tol);
  EXPECT_LE(sig.l1_ueden, tol);
  EXPECT_LE(sig.l1_others, tol);
}

TEST(ODTLESLocalEngine, RuntimeODTMomentumDepositionMPIParitySignature)
{
  constexpr amrex::Real tol = 1.0e-12;
  amrex::MultiFab lterm_runtime;
  const RuntimeSignature sig = computeRuntimeODTLESTermSignature(lterm_runtime);

  // In-test guardrails; cross-run parity is done by CTest wrapper.
  EXPECT_GT(sig.stepped_entries, 0);
  EXPECT_GT(sig.moment_columns_built, 0);
  EXPECT_LE(sig.l1_ueden, tol);
  EXPECT_LE(sig.l1_others, tol);

  writeRuntimeSignatureIfRequested(sig);
}

TEST(ODTLESLocalEngine, ConservativeMomentumDepositionMultiRankStyleBudget)
{
  constexpr amrex::Real tol = 1.0e-12;
  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(63, 0, 0)));

  amrex::BoxArray ba(domain);
  ba.maxSize(4);
  amrex::DistributionMapping dm(ba);
  amrex::MultiFab lterm;
  computeLtermFromCenteredContributionMF(ba, dm, domain, lterm);

  amrex::Real sum_umx = 0.0;
  amrex::Real sum_umy = 0.0;
  amrex::Real sum_umz = 0.0;
  amrex::Real sum_ueden = 0.0;
  for (amrex::MFIter mfi(lterm, false); mfi.isValid(); ++mfi) {
    const auto vbx = mfi.validbox();
    const auto arr = lterm.const_array(mfi);
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv)) {
      sum_umx += arr(iv, UMX);
      sum_umy += arr(iv, UMY);
      sum_umz += arr(iv, UMZ);
      sum_ueden += arr(iv, UEDEN);
      EXPECT_LE(std::abs(arr(iv, URHO)), tol);
      for (int n = UFS; n < UFS + NUM_SPECIES; ++n) {
        EXPECT_LE(std::abs(arr(iv, n)), tol);
      }
#if NUM_AUX > 0
      for (int n = UFX; n < UFX + NUM_AUX; ++n) {
        EXPECT_LE(std::abs(arr(iv, n)), tol);
      }
#endif
#if NUM_ADV > 0
      for (int n = UFA; n < UFA + NUM_ADV; ++n) {
        EXPECT_LE(std::abs(arr(iv, n)), tol);
      }
#endif
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(sum_umx);
  amrex::ParallelDescriptor::ReduceRealSum(sum_umy);
  amrex::ParallelDescriptor::ReduceRealSum(sum_umz);
  amrex::ParallelDescriptor::ReduceRealSum(sum_ueden);

  // For qx=[x,-x,2x,0] and one-sided physical boundaries used in the current directional momentum deposition mapping:
  // sum(D) = -(F_right - F_left) = [-(63-0), +(63-0), -2*(63-0)].
  EXPECT_NEAR(sum_umx, -63.0, tol);
  EXPECT_NEAR(sum_umy, 63.0, tol);
  EXPECT_NEAR(sum_umz, -126.0, tol);
  EXPECT_NEAR(sum_ueden, 0.0, tol);
}

TEST(ODTLESLocalEngine, ConservativeMomentumDepositionAMRCoarseFineInterface)
{
  constexpr amrex::Real tol = 1.0e-12;

  // Coarse level: 8 cells over [0,8], dx=1.0.
  const amrex::Box coarse_domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(7, 0, 0)));
  // Fine level (r=2) logical full domain [0,16), but only right-half patch active.
  const amrex::Box fine_domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(15, 0, 0)));
  const amrex::Box fine_patch(
    amrex::IntVect(AMREX_D_DECL(8, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(15, 0, 0)));

  // Coarse buffers.
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> coarse_q{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::grow(coarse_domain, 1), 4),
      amrex::FArrayBox(amrex::grow(coarse_domain, 1), 4),
      amrex::FArrayBox(amrex::grow(coarse_domain, 1), 4))};
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> coarse_flux{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::surroundingNodes(coarse_domain, 0), NVAR),
      amrex::FArrayBox(amrex::surroundingNodes(coarse_domain, 1), NVAR),
      amrex::FArrayBox(amrex::surroundingNodes(coarse_domain, 2), NVAR))};
  amrex::FArrayBox coarse_lterm(coarse_domain, NVAR);
  amrex::FArrayBox coarse_vol(coarse_domain, 1);
  coarse_vol.setVal(1.0);

  computeSingleLevelLinearDeposition(
    coarse_domain, coarse_domain, 1.0, coarse_q, coarse_flux, coarse_lterm);

  // Fine patch buffers (with 1-cell ghost so interface face can use central avg).
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> fine_q{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::grow(fine_patch, 1), 4),
      amrex::FArrayBox(amrex::grow(fine_patch, 1), 4),
      amrex::FArrayBox(amrex::grow(fine_patch, 1), 4))};
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> fine_flux{
    AMREX_D_DECL(
      amrex::FArrayBox(amrex::surroundingNodes(fine_patch, 0), NVAR),
      amrex::FArrayBox(amrex::surroundingNodes(fine_patch, 1), NVAR),
      amrex::FArrayBox(amrex::surroundingNodes(fine_patch, 2), NVAR))};
  amrex::FArrayBox fine_lterm(fine_patch, NVAR);
  amrex::FArrayBox fine_vol(fine_patch, 1);
  fine_vol.setVal(0.5);

  computeSingleLevelLinearDeposition(
    fine_patch, fine_domain, 0.5, fine_q, fine_flux, fine_lterm);

  // Coarse-fine interface flux continuity at x=4:
  // coarse face i=4 equals fine face i=8 (same physical face).
  const auto cfx = coarse_flux[0].const_array();
  const auto ffx = fine_flux[0].const_array();
  EXPECT_NEAR(cfx(4, 0, 0, UMX), ffx(8, 0, 0, UMX), tol);
  EXPECT_NEAR(cfx(4, 0, 0, UMY), ffx(8, 0, 0, UMY), tol);
  EXPECT_NEAR(cfx(4, 0, 0, UMZ), ffx(8, 0, 0, UMZ), tol);
  EXPECT_NEAR(cfx(4, 0, 0, UEDEN), 0.0, tol);
  EXPECT_NEAR(ffx(8, 0, 0, UEDEN), 0.0, tol);

  // Composite AMR active-state budget:
  // active coarse cells [0..3] + active fine cells [8..15].
  // Internal coarse-fine interface contributions cancel when interface fluxes match.
  const amrex::Box coarse_active(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(3, 0, 0)));
  const amrex::Real sum_umx =
    sumWeightedComponent(coarse_lterm, coarse_vol, coarse_active, UMX) +
    sumWeightedComponent(fine_lterm, fine_vol, fine_patch, UMX);
  const amrex::Real sum_umy =
    sumWeightedComponent(coarse_lterm, coarse_vol, coarse_active, UMY) +
    sumWeightedComponent(fine_lterm, fine_vol, fine_patch, UMY);
  const amrex::Real sum_umz =
    sumWeightedComponent(coarse_lterm, coarse_vol, coarse_active, UMZ) +
    sumWeightedComponent(fine_lterm, fine_vol, fine_patch, UMZ);

  // Composite external boundaries: left coarse face i=0 and right fine face i=16.
  const amrex::Real rhs_umx =
    -(ffx(16, 0, 0, UMX) - cfx(0, 0, 0, UMX));
  const amrex::Real rhs_umy =
    -(ffx(16, 0, 0, UMY) - cfx(0, 0, 0, UMY));
  const amrex::Real rhs_umz =
    -(ffx(16, 0, 0, UMZ) - cfx(0, 0, 0, UMZ));
  EXPECT_NEAR(sum_umx, rhs_umx, tol);
  EXPECT_NEAR(sum_umy, rhs_umy, tol);
  EXPECT_NEAR(sum_umz, rhs_umz, tol);

  // Momentum-only targeting across both active regions.
  const auto cd = coarse_lterm.const_array();
  for (amrex::IntVect iv = coarse_active.smallEnd(); iv <= coarse_active.bigEnd();
       coarse_active.next(iv)) {
    EXPECT_LE(std::abs(cd(iv, UEDEN)), tol);
    EXPECT_LE(std::abs(cd(iv, URHO)), tol);
    for (int n = UFS; n < UFS + NUM_SPECIES; ++n) {
      EXPECT_LE(std::abs(cd(iv, n)), tol);
    }
#if NUM_AUX > 0
    for (int n = UFX; n < UFX + NUM_AUX; ++n) {
      EXPECT_LE(std::abs(cd(iv, n)), tol);
    }
#endif
#if NUM_ADV > 0
    for (int n = UFA; n < UFA + NUM_ADV; ++n) {
      EXPECT_LE(std::abs(cd(iv, n)), tol);
    }
#endif
  }
  const auto fd = fine_lterm.const_array();
  for (amrex::IntVect iv = fine_patch.smallEnd(); iv <= fine_patch.bigEnd();
       fine_patch.next(iv)) {
    EXPECT_LE(std::abs(fd(iv, UEDEN)), tol);
    EXPECT_LE(std::abs(fd(iv, URHO)), tol);
    for (int n = UFS; n < UFS + NUM_SPECIES; ++n) {
      EXPECT_LE(std::abs(fd(iv, n)), tol);
    }
#if NUM_AUX > 0
    for (int n = UFX; n < UFX + NUM_AUX; ++n) {
      EXPECT_LE(std::abs(fd(iv, n)), tol);
    }
#endif
#if NUM_ADV > 0
    for (int n = UFA; n < UFA + NUM_ADV; ++n) {
      EXPECT_LE(std::abs(fd(iv, n)), tol);
    }
#endif
  }
}

} // namespace pelec_tests

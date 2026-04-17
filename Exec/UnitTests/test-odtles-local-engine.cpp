#include "gtest/gtest.h"

#include <cmath>
#include <cstdlib>
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
#include "ODTRuntimeLESBridge.H"
#include "ODTStepper.H"
#include "ODTTripletMap.H"

namespace pelec_tests
{
namespace
{

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
computeRuntimeODTLESTermSignature(amrex::MultiFab& lterm_out)
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

  constexpr int support_ng =
    (pelec::odtles::ODTLineGeometry::MVPNumSupportCells - 1) / 2;
  amrex::MultiFab state_valid(ba, dm, NVAR, 0);
  amrex::MultiFab state_same_level(ba, dm, NVAR, support_ng);
  amrex::MultiFab state_host_filled(ba, dm, NVAR, support_ng);
  state_valid.setVal(0.0);
  state_same_level.setVal(0.0);
  state_host_filled.setVal(0.0);

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
      s(iv, UEDEN) = 2.5 + 0.02 * x + 0.01 * y;
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
  pelec::odtles::ODTParams params{};
  params.enabled = true;
  params.max_local_substeps = 1;
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
  EXPECT_TRUE(rep.tau_zero_for_single_segment_central_interval);
  EXPECT_LE(rep.overlap_weight_sum_error, tol);
  EXPECT_LE(rep.owner_recovery_max_abs_error, tol);
  EXPECT_LE(rep.tau_max_abs, tol);
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

  // For qx=[x,-x,2x,0] and one-sided physical boundaries used in T6 mapping:
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

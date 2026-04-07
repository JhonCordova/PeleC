#include "gtest/gtest.h"

#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>

#include "Diffterm.H"
#include "IndexDefines.H"
#include "ODTDiffusion.H"
#include "ODTMomentExtractor.H"
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
  EXPECT_NEAR(lterm(iv1, UMX), -1.0, tol);
  EXPECT_NEAR(lterm(iv2, UMX), -1.0, tol);
  EXPECT_NEAR(lterm(iv1, UMY), 1.0, tol);
  EXPECT_NEAR(lterm(iv2, UMY), 1.0, tol);
  EXPECT_NEAR(lterm(iv1, UMZ), -2.0, tol);
  EXPECT_NEAR(lterm(iv2, UMZ), -2.0, tol);

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

} // namespace pelec_tests

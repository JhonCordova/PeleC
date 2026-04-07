#include "gtest/gtest.h"

#include <cmath>

#include "ODTDiffusion.H"
#include "ODTMomentExtractor.H"
#include "ODTStepper.H"
#include "ODTTripletMap.H"

namespace pelec_tests
{

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

} // namespace pelec_tests

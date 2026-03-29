#include "gtest/gtest.h"

#include <cmath>

#include "ODTDiffusion.H"
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

} // namespace pelec_tests

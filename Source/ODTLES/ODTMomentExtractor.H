#ifndef PELEC_ODTLES_ODTMOMENTEXTRACTOR_H
#define PELEC_ODTLES_ODTMOMENTEXTRACTOR_H

#include <array>

#include <AMReX_REAL.H>

#include "ODTLineGeometry.H"
#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTMomentExtractor
{
public:
  struct DirectionalColumn
  {
    int dir = -1; // Fixed line direction j for tau_ij.
    bool valid = false;

    // Geometric averaging metadata over owner-cell central interval.
    amrex::Real owner_interval_length = 0.0;
    amrex::Real overlap_measure = 0.0;
    amrex::Real overlap_weight_sum = 0.0;
    int contributing_segments = 0;

    // Favre-consistent moments over central interval:
    // <rho>, <rho u_i>, <rho u_i u_j> for fixed j=dir.
    amrex::Real mean_rho = 0.0;
    std::array<amrex::Real, AMREX_SPACEDIM> mean_rho_ui{};
    std::array<amrex::Real, AMREX_SPACEDIM> mean_rho_ui_uj{};

    // Directional SGS column at fixed j:
    // tau_ij = <rho u_i u_j> - <rho u_i><rho u_j>/<rho>.
    std::array<amrex::Real, AMREX_SPACEDIM> tau_ij{};
  };

  struct ValidationReport
  {
    bool overlap_weights_partition_unity = false;
    bool owner_central_recovered = false;
    bool tau_zero_for_single_segment_central_interval = false;
    amrex::Real overlap_weight_sum_error = 0.0;
    amrex::Real owner_recovery_max_abs_error = 0.0;
    amrex::Real tau_max_abs = 0.0;
  };

  struct DirectionalMomentumContribution
  {
    static constexpr int MomX = 0;
    static constexpr int MomY = 1;
    static constexpr int MomZ = 2;
    static constexpr int Ueden = 3;
    static constexpr int NumComponents = 4;

    // Fixed line direction j carried explicitly for downstream mapping.
    int dir = -1;
    bool valid = false;
    // Layout-compatible momentum/energy-like slots:
    // [UMX, UMY, UMZ, UEDEN], with UEDEN intentionally inactive for now.
    std::array<amrex::Real, NumComponents> q{{0.0, 0.0, 0.0, 0.0}};
  };

  // Extract directional momentum moments from the owner-cell central interval.
  //
  // Discrete averaging uses geometric overlap weights over support segments:
  // Let I0 be the owner central interval and Ik each support segment interval.
  // For each intersecting segment k:
  //   wk = |I0 ∩ Ik| / |I0|
  // and moments are:
  //   <q> = sum_k wk * qk
  // where qk is the segment piecewise-constant conservative state.
  // For fixed line direction j:
  //   <rho>, <rho u_i>, <rho u_i u_j>
  // and
  //   tau_ij = <rho u_i u_j> - <rho u_i><rho u_j>/<rho>.
  static DirectionalColumn extractDirectionalMomentumColumn(
    const ODTLineGeometry& geom,
    const ODTLineState& line_state);

  // Build directional intermediate momentum SGS contribution from tau_ij
  // for fixed line direction j:
  //   q = [-tau_1j, -tau_2j, -tau_3j, 0]
  // This is a cell-centered intermediate representation only; no face mapping,
  // area scaling, divergence assembly, or LESTerm deposition is performed here.
  static DirectionalMomentumContribution buildDirectionalMomentumContribution(
    const DirectionalColumn& column);

  // Lightweight validation hook for overlap-weighted central extraction:
  // checks partition-of-unity weights and owner-central recovery behavior.
  static ValidationReport runMVPValidationHook();
};

} // namespace pelec::odtles

#endif

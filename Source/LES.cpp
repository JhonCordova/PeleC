// Constant and Dynamic Smagorinsky models are implemented
//
// AT THE MOMENT:
// - no multifluid component of the equations are solved
//
// - Strictly, the formulation assumes perfect gas / constant specific heats /
// GammaLaw EoS.
//     This assumption (e + p/rho = h = Cp*T) is used to relate the turbulent
//     transport of internal energy + pressure term to the turbulent temperature
//     transport and thus the temperature gradient:
//
//     with <> = filter, {} = density weighted filter, state = [rho, rho*E,
//     rho*u_j] <rho>{u_j h} - <rho>{u_j}{h}  = Smagorinsky Model = <rho>
//     nu_T/Pr_T * d{h}/dx_j
//                                     -> modeled as ->  <rho> nu_T/Pr_T *
//                                     Cp(<state>) * dT(<state>)/dx_j
//
//     For non-perfect gasses there are two issues with the formulation: h =/=
//     Cp*T, and T is a nonlinear function of the state variables so T(<state>)
//     =/= {T}. For ideal gases where the specific heats are relatively weak
//     functions of temperaure, the error induced from both of these issues is
//     expected to be small relative to the overall error in the Smagorinsky
//     models because the error in the approximation delta_h = Cp*delta_T is ~
//     dCp/dT*delta_T. But the present formulation is unsuited for real gas
//     equations and superctritical fluids where properties such as Cp are very
//     strong nonlinear functions of the state variables. Converting this term
//     to an enthalpy-based Smagorinsky-like model would alleviate the former
//     issue, but the latter would still remain: h(<state>) =/= {h}.
//
// - The formulation also assumes that the pressure term in the velocity
// equation can be solved
//     using p(<state>), rather than <p>, which are not equivalent for
//     non-perfect gasses (even for ideal gasses where <p> = <rho>*R*{T}
//     (assuming R is constant) {T} is a non-linear function of the state, so
//     this makes p(<state>) =/= <p>. Again, for weak nonlinearity this is
//     probably not the leading order error in the modeling approach. An
//     alternative explanation is that any discrepancy is part of what is
//     modeled by the Smagorinsky approach.
//
// - The last major modeling assumption is that molecular transport may be
// closed using, for example,
//     <q> = < -lambda dT/dx > = -lambda(<state>) d{T}/dx
//     This type of closure is standard in the literature, or at least similar
//     closures using <lambda> instead of lambda(<state>) are standard. Although
//     perhaps it has not been verified for real gas EoS. In any event, for
//     coarse LES grids the molecular transport coeffs are much smaller than the
//     turbulent one, so this is standard modeling approach is probably
//     justified.

#include "LES.H"
#include "ODTMomentExtractor.H"
#include "ODTStepper.H"

void
PeleC::construct_old_les_source(
  amrex::Real time, amrex::Real dt, int /*sub_iteration*/, int /*sub_ncycle*/)
{
  // Add grow cells necessary for explicit filtering of source terms
  if (use_explicit_filter) {
    filtered_les_source.define(
      grids, dmap, NVAR, old_sources[les_src]->nGrow(), amrex::MFInfo(),
      Factory());
    old_sources[les_src]->define(
      grids, dmap, NVAR, old_sources[les_src]->nGrow() + nGrowF,
      amrex::MFInfo(), Factory());
  }

  old_sources[les_src]->setVal(0.0);

  amrex::Real reflux_factor_old = 0.5;
  getLESTerm(time, dt, *old_sources[les_src], reflux_factor_old);

  old_sources[les_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_les_source(
  amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle)
{
  // Add grow cells necessary for explicit filtering of source terms
  if (use_explicit_filter) {
    filtered_les_source.define(
      grids, dmap, NVAR, new_sources[les_src]->nGrow(), amrex::MFInfo(),
      Factory());
    new_sources[les_src]->define(
      grids, dmap, NVAR, new_sources[les_src]->nGrow() + nGrowF,
      amrex::MFInfo(), Factory());
  }

  new_sources[les_src]->setVal(0.0);

  amrex::Real reflux_factor_new = sub_iteration == sub_ncycle - 1 ? 0.5 : 0;
  getLESTerm(time, dt, *new_sources[les_src], reflux_factor_new);
}

// Calculate the LES term by calling an SFS model
// Across all conserved state components, compute the LES "source term"
//    = -Div(LESFlux).
void
PeleC::getLESTerm(
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  amrex::Real reflux_factor)
{
  BL_PROFILE("PeleC::getLESTerm()");

  if (!do_les) {
    LESTerm.setVal(0, 0, NVAR, LESTerm.nGrow());
    return;
  }

  if (verbose != 0) {
    amrex::Print() << "... Computing LES term at time " << time << std::endl;
  }

  switch (les_model) {

  case PeleC::LESModelSmagorinsky:
    getSmagorinskyLESTerm(time, dt, LESTerm, reflux_factor);
    break;

  case PeleC::LESModelDynamicSmagorinsky:
    getDynamicSmagorinskyLESTerm(time, dt, LESTerm, reflux_factor);
    break;

  case PeleC::LESModelWALE:
    getWALELESTerm(time, dt, LESTerm, reflux_factor);
    break;

  case PeleC::LESModelVreman:
    getVremanLESTerm(time, dt, LESTerm, reflux_factor);
    break;

  case PeleC::LESModelODT:
    getODTLESTerm(time, dt, LESTerm, reflux_factor);
    break;

  default:
    amrex::Error("Invalid les_model number.");
    break;
  }

  amrex::Print()
    << "WARNING -- Need to implement LES extrap here, see diffusion"
    << std::endl;

  // Filter the SGS source term
  if (use_explicit_filter) {
    les_filter.apply_filter(LESTerm, filtered_les_source);
    LESTerm.define(
      grids, dmap, NVAR, filtered_les_source.nGrow(), amrex::MFInfo(),
      Factory());
    amrex::MultiFab::Copy(
      LESTerm, filtered_les_source, 0, 0, NVAR, filtered_les_source.nGrow());
  }
}

void
PeleC::getODTLESTerm(
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  amrex::Real reflux_factor)
{
  // ODT runtime hookup:
  // 1) prepare/reconcile local lines,
  // 2) advance each line over host-owned dt_LES,
  // 3) extract directional SGS momentum moments and map to directional face
  //    fluxes, then use standard PeleC conservative divergence assembly.
  odt_manager.initializeLevel(level, grids, dmap);

  constexpr int support_ng =
    (pelec::odtles::ODTLineGeometry::MVPNumSupportCells - 1) / 2;

  // T4 provenance policy:
  // 1) Build valid level state at requested time.
  // 2) Build same-level ghost extension via FillBoundary only.
  // This allows same-level inter-rank/periodic support collection while
  // explicitly avoiding AMR coarse-fine ghost provenance for support means.
  amrex::MultiFab state_valid(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
  FillPatch(*this, state_valid, 0, time, State_Type, 0, NVAR);

  amrex::MultiFab state_same_level(
    grids, dmap, NVAR, support_ng, amrex::MFInfo(), Factory());
  amrex::MultiFab::Copy(state_same_level, state_valid, 0, 0, NVAR, 0);
  state_same_level.FillBoundary(geom.periodicity());

  // Only for provenance surfacing: this contains host-filled ghosts that may
  // include AMR/coarse-fine fill mechanisms. It is not the support source used
  // by the T4 math path unless policy later accepts that provenance.
  amrex::MultiFab state_host_filled(
    grids, dmap, NVAR, support_ng, amrex::MFInfo(), Factory());
  FillPatch(*this, state_host_filled, support_ng, time, State_Type, 0, NVAR);

  // Cell-centered directional intermediate from ODT for each fixed line
  // direction j. Component layout:
  //   [0]=UMX contribution = -tau_1j
  //   [1]=UMY contribution = -tau_2j
  //   [2]=UMZ contribution = -tau_3j
  //   [3]=UEDEN placeholder = 0 (inactive in this phase)
  constexpr int odt_flux_ncomp =
    pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::NumComponents;
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> odt_dir_cc;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    odt_dir_cc[dir].define(
      grids, dmap, odt_flux_ncomp, 1, amrex::MFInfo(), Factory());
    odt_dir_cc[dir].setVal(0.0, 0, odt_flux_ncomp, odt_dir_cc[dir].nGrow());
  }

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

  for (amrex::MFIter mfi(state_valid, false); mfi.isValid(); ++mfi) {
    const amrex::Box& vbx = mfi.validbox();
    for (amrex::IntVect iv = vbx.smallEnd(); iv <= vbx.bigEnd(); vbx.next(iv))
    {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        pelec::odtles::ODTLineGeometry line_geom(geom, level, iv, dir);
        const auto support_data = odt_manager.collectSupportData(
          line_geom, geom, state_valid, state_same_level, state_host_filled);

        const int n_boundary = support_data.count(
          pelec::odtles::ODTManager::SupportProvenance::PhysicalBoundaryGhost);
        const int n_amr = support_data.count(
          pelec::odtles::ODTManager::SupportProvenance::AMRCoarseFineFilled);
        const int n_invalid = support_data.count(
          pelec::odtles::ODTManager::SupportProvenance::InvalidUnsupported);
        const int n_reject_types =
          (n_boundary > 0 ? 1 : 0) + (n_amr > 0 ? 1 : 0) + (n_invalid > 0 ? 1 : 0);

        // Explicit T4 runtime provenance policy:
        // allow: SameLevelValidCell, SameLevelFilledGhost
        // reject: PhysicalBoundaryGhost, AMRCoarseFineFilled, InvalidUnsupported
        if (n_reject_types > 1) {
          ++rejected_mixed_entries;
          continue;
        }
        if (n_boundary > 0) {
          ++rejected_boundary_entries;
          continue;
        }
        if (n_amr > 0) {
          ++rejected_amr_entries;
          continue;
        }
        if (n_invalid > 0) {
          ++rejected_invalid_entries;
          continue;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          support_data.allSameLevelAccepted(),
          "ODTLES support provenance policy mismatch in runtime hook");

        auto* entry = odt_manager.findLineEntry(level, iv, dir);
        if (entry != nullptr && entry->state.valid()) {
          odt_manager.reconcileLineStateFromSupportData(
            level, iv, dir, geom, support_data);
        } else {
          odt_manager.initializeLineStateFromSupportData(
            level, iv, dir, geom, support_data);
        }

        auto* stepped_entry = odt_manager.findLineEntry(level, iv, dir);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          stepped_entry != nullptr && stepped_entry->state.valid(),
          "ODTLES runtime expected a valid persistent line before local stepping");

        pelec::odtles::ODTStepper::Controls step_ctrl{};
        step_ctrl.max_internal_iterations =
          odt_manager.params().max_local_substeps > 0
            ? odt_manager.params().max_local_substeps
            : 1;
        step_ctrl.sampler_controls.deterministic = true;
        step_ctrl.sampler_controls.seed = static_cast<std::uint64_t>(
          (static_cast<unsigned>(level) + 1U) * 73856093U ^
          (static_cast<unsigned>(dir) + 1U) * 19349663U ^
          (static_cast<unsigned>(iv[0]) + 1U) * 83492791U ^
          (static_cast<unsigned>(iv[1]) + 1U) * 2654435761U ^
          (static_cast<unsigned>(iv[2]) + 1U) * 97531U);

        const auto step_rep = pelec::odtles::ODTStepper::advanceOneLESTimestep(
          stepped_entry->geometry, stepped_entry->state, dt, step_ctrl);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          step_rep.reached_dt_les,
          "ODTLES local stepper failed to close dt_LES exactly");

        const auto moment_column =
          pelec::odtles::ODTMomentExtractor::extractDirectionalMomentumColumn(
            stepped_entry->geometry, stepped_entry->state);
        const auto momentum_sgs_intermediate =
          pelec::odtles::ODTMomentExtractor::buildDirectionalMomentumContribution(
            moment_column);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          momentum_sgs_intermediate.valid &&
            momentum_sgs_intermediate.dir == dir,
          "ODTLES directional SGS intermediate metadata mismatch");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          momentum_sgs_intermediate.q[pelec::odtles::ODTMomentExtractor::
                                         DirectionalMomentumContribution::Ueden] == 0.0,
          "ODTLES directional SGS intermediate must keep UEDEN inactive");

        auto const& dir_cc = odt_dir_cc[dir].array(mfi);
        for (int n = 0; n < odt_flux_ncomp; ++n) {
          dir_cc(iv, n) = momentum_sgs_intermediate.q[static_cast<std::size_t>(n)];
        }
        ++moment_columns_built;

        ++stepped_entries;
        stepper_attempted_events += step_rep.attempted_events;
        stepper_applied_events += step_rep.applied_events;
        stepper_rejected_events += step_rep.rejected_events;
        if (step_rep.closed_by_diffusion_only_catchup) {
          ++stepper_diffusion_only_catchup;
        }
      }
    }
  }

  const long rejected_total = rejected_boundary_entries + rejected_amr_entries +
                              rejected_invalid_entries + rejected_mixed_entries;
  if (verbose != 0 && rejected_total > 0) {
    amrex::Print()
      << "ODTLES: rejected " << rejected_total
      << " owner-direction entries by explicit T4 support-provenance policy. "
      << "[boundary_ghost=" << rejected_boundary_entries
      << ", amr_coarse_fine=" << rejected_amr_entries
      << ", invalid=" << rejected_invalid_entries
      << ", mixed=" << rejected_mixed_entries << "]" << std::endl;
  }
  if (verbose != 0 && stepped_entries > 0) {
    amrex::Print() << "ODTLES: stepped " << stepped_entries
                   << " local lines over dt_LES with neutral SGS return. "
                   << "[attempted_events=" << stepper_attempted_events
                   << ", applied_events=" << stepper_applied_events
                   << ", rejected_events=" << stepper_rejected_events
                    << ", diffusion_only_catchup_entries="
                   << stepper_diffusion_only_catchup
                   << ", directional_momentum_columns="
                   << moment_columns_built << "]" << std::endl;
  }

  // Fill same-level/periodic ghosts for centered directional intermediates.
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    odt_dir_cc[dir].FillBoundary(geom.periodicity());
  }

  // Map centered directional ODT SGS contributions to face fluxes in the
  // matching face-normal direction and then deposit conservatively.
  // For fixed direction j and momentum component i:
  //   centered q_i^(j) = -tau_ij
  //   face flux F_i^(j) = A_j * avg_face(q_i^(j))
  // with avg_face as arithmetic average between adjacent centered values.
  // UEDEN is explicitly kept inactive: F_UEDEN^(j)=0.
  const auto domain = geom.Domain();
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(LESTerm.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(LESTerm, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      const amrex::Box cbox = amrex::grow(vbox, 0);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(cbox);
      if (typ != amrex::FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
      if (typ == amrex::FabType::covered) {
        continue;
      }

      const amrex::Array<const amrex::Box, AMREX_SPACEDIM> eboxes = {
        AMREX_D_DECL(
          amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
          amrex::surroundingNodes(cbox, 2))};

      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec;
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      resizeAndSetFlux(eboxes, flx, flux_ec);

      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        a{{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};

      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto qcc = odt_dir_cc[dir].const_array(mfi);
        amrex::ParallelFor(
          eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const amrex::IntVect ivm =
              iv - amrex::IntVect::TheDimensionVector(dir);
            const bool has_hi = domain.contains(iv);
            const bool has_lo = domain.contains(ivm);

            const int il = i - (dir == 0 ? 1 : 0);
            const int jl = j - (dir == 1 ? 1 : 0);
            const int kl = k - (dir == 2 ? 1 : 0);

            const auto centered_face_avg = [&](int comp) -> amrex::Real {
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

            const amrex::Real area_j = a[dir](i, j, k);
            flx[dir](i, j, k, UMX) = area_j * centered_face_avg(
              pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomX);
            flx[dir](i, j, k, UMY) = area_j * centered_face_avg(
              pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomY);
            flx[dir](i, j, k, UMZ) = area_j * centered_face_avg(
              pelec::odtles::ODTMomentExtractor::DirectionalMomentumContribution::MomZ);
            flx[dir](i, j, k, UEDEN) = 0.0;
          });
      }

      computeFluxDiv(LESTerm, mfi, vbox, flx, volume);
      updateFluxRegistersLES(
        reflux_factor, LESTerm, dt, mfi, typ,
        {AMREX_D_DECL(flux_ec.data(), &flux_ec[1], &flux_ec[2])});
    }
  }
}

void
computeTangentialVelDerivs(
  const amrex::Array<const amrex::Box, AMREX_SPACEDIM>& eboxes,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& tander_ec,
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& tanders,
  const amrex::Array4<amrex::Real>& q_ar,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
  BL_PROFILE("PeleC::pc_compute_tangential_vel_derivs()");
#if AMREX_SPACEDIM == 3
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    tander_ec[dir].resize(
      eboxes[dir], GradUtils::nCompTan, amrex::The_Async_Arena());
    tanders[dir] = tander_ec[dir].array();
    setV(eboxes[dir], GradUtils::nCompTan, tanders[dir], 0);
    const amrex::Real d1 = dir == 0 ? dx[1] : dx[0];
    const amrex::Real d2 = dir == 2 ? dx[1] : dx[2];
    amrex::ParallelFor(
      eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_compute_tangential_vel_derivs(
          i, j, k, q_ar, dir, d1, d2, tanders[dir]);
      });
  }
#else
  amrex::ignore_unused(eboxes, tander_ec, tanders, q_ar, dx);
  amrex::Abort("computeTangentialVelDerivs: only supported in 3D");
#endif
}

void
resizeAndSetFlux(
  const amrex::Array<const amrex::Box, AMREX_SPACEDIM>& eboxes,
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& flux_ec)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    flux_ec[dir].resize(eboxes[dir], NVAR, amrex::The_Async_Arena());
    flx[dir] = flux_ec[dir].array();
    setV(eboxes[dir], NVAR, flx[dir], 0);
  }
}

void
computeFluxDiv(
  amrex::MultiFab& LESTerm,
  const amrex::MFIter& mfi,
  const amrex::Box& vbox,
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::MultiFab& volume)
{
  auto const& Lterm = LESTerm.array(mfi);
  setV(vbox, NVAR, Lterm, 0.0);
  {
    BL_PROFILE("PeleC::pc_flux_div()");
    auto const& vol = volume.array(mfi);
    amrex::ParallelFor(
      vbox, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        pc_flux_div(
          i, j, k, n, AMREX_D_DECL(flx[0], flx[1], flx[2]), vol, Lterm);
      });
  }
}

void
PeleC::updateFluxRegistersLES(
  amrex::Real reflux_factor,
  const amrex::MultiFab& LESTerm,
  amrex::Real dt,
  const amrex::MFIter& mfi,
  amrex::FabType typ,
  const amrex::Array<amrex::FArrayBox const*, AMREX_SPACEDIM>& flux_ec)
{
  if (do_reflux && reflux_factor != 0) {
    amrex::FArrayBox dm_as_fine(
      amrex::Box::TheUnitBox(), LESTerm.nComp(), amrex::The_Async_Arena());
    update_flux_registers(reflux_factor * dt, mfi, typ, flux_ec, dm_as_fine);
  }
}

// Calculate the LES term using the Smagorinsky SFS model
void
PeleC::getSmagorinskyLESTerm(
#if AMREX_SPACEDIM < 3
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  amrex::MultiFab& /*LESTerm*/,
  amrex::Real /*reflux_factor*/)
{
  amrex::Abort("LES only implemented in 3D for now");
#else
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  amrex::Real reflux_factor)
{
  // Only use this functionality for 3D
  const int ngrow = 1;
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

  amrex::MultiFab S(grids, dmap, NVAR, ngrow, amrex::MFInfo(), Factory());
  FillPatch(*this, S, ngrow, time, State_Type, 0, NVAR);

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      const amrex::Box gbox = amrex::grow(vbox, ngrow);
      const amrex::Box cbox = amrex::grow(vbox, ngrow - 1);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(cbox);
      if (typ != amrex::FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
      if (typ == amrex::FabType::covered) {
        continue;
      }

      auto const& s = S.array(mfi);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      amrex::FArrayBox q(gbox, QVAR, amrex::The_Async_Arena());
      amrex::FArrayBox qaux(gbox, nqaux, amrex::The_Async_Arena());
      auto const& q_ar = q.array();
      auto const& qauxar = qaux.array();

      // Get primitives, Q, including (Y, T, p, rho) from conserved state
      // required for L term
      {
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, s, q_ar, qauxar);
          });
      }

      // Get the tangential derivatives
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> tander_ec;
      const amrex::Array<const amrex::Box, AMREX_SPACEDIM> eboxes = {
        AMREX_D_DECL(
          amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
          amrex::surroundingNodes(cbox, 2))};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> tanders;

      computeTangentialVelDerivs(eboxes, tander_ec, tanders, q_ar, dx);

      // Compute extensive LES fluxes, F.A
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec;
      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        a{{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      resizeAndSetFlux(eboxes, flx, flux_ec);
      {
        BL_PROFILE("PeleC::pc_smagorinsky_sfs_term()");
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          amrex::Real Cs_local = PeleC::Cs;
          amrex::Real CI_local = PeleC::CI;
          amrex::Real PrT_local = PeleC::PrT;
          amrex::ParallelFor(
            eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_smagorinsky_sfs_term(
                i, j, k, q_ar, tanders[dir], a[dir], dx, dir, Cs_local,
                CI_local, PrT_local, flx[dir]);
            });
        }
      }

      // Compute flux divergence (1/Vol).Div(F.A)
      computeFluxDiv(LESTerm, mfi, vbox, flx, volume);

      // Refluxing
      updateFluxRegistersLES(
        reflux_factor, LESTerm, dt, mfi, typ,
        {AMREX_D_DECL(flux_ec.data(), &flux_ec[1], &flux_ec[2])});
    } // End of MFIter scope
  } // End of OMP scope
#endif
}

// Calculate the LES term using the dynamic Smagorinsky SFS model
void
PeleC::getDynamicSmagorinskyLESTerm(
#if AMREX_SPACEDIM < 3
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  amrex::MultiFab& /*LESTerm*/,
  amrex::Real /*reflux_factor*/)
{
  amrex::Abort("LES only implemented in 3D for now");
#else
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  amrex::Real reflux_factor)
{
  // clang-format off
  /*
    Note on the grow cells:
    N + 0                             (cbox) |----->         LESTerm
    N + 0**                           (ebox) |----->         flux_ec, coeff_ec, alphaij_ec, alpha_ec, flux_T_ec
    N + 1                            (g4box) |------>        filtered_coeff_cc [= LES_Coeffs]
    N + 1 + nGrowC                   (g3box) |-------->      coeff_cc, filtered_(K, RUT, alphaij, alpha, flux_T)
    N + 1 + nGrowC + nGrowD          (g2box) |---------->    filtered_(S, Q, Qaux)
    N + 1 + nGrowC + nGrowT          (g1box) |----------->   K, RUT, alphaij, alpha, flux_T
    N + 1 + nGrowC + nGrowT + nGrowD (g0box) |-------------> S, Q, Qaux
       |----------------------------|
       This is the number of grow cells on each side

    where
    nGrowD = number of grow cells necessary for the diffusion operator
    nGrowC = number of grow cells necessary for filtering the Smagorinsky coefficients
    nGrowT = number of grow cells necessary for filtering the derived quantities (test level)

    ** Everything is calculated at cell centers (cc), then at the very end the coefficients and the stress terms (e.g. alpha_ij)
       are moved to edge/faces centers (ec) to calculate the fluxes. ec quantities have length N+1 in the face-normal
       direction and length N in the other two directions.

  */
  // clang-format on
  Filter test_filter = Filter(les_test_filter_type, les_test_filter_fgr);
  Filter coeff_filter = Filter(box, 6);

  const int nGrowD = 1;
  const int nGrowC = coeff_filter.get_filter_ngrow();
  const int nGrowT = test_filter.get_filter_ngrow();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

  // 1. Get state variable data
  amrex::MultiFab S(
    grids, dmap, NVAR, nGrowD + nGrowC + nGrowT + 1, amrex::MFInfo(),
    Factory());
  FillPatch(*this, S, nGrowD + nGrowC + nGrowT + 1, time, State_Type, 0, NVAR);
  LES_Coeffs.setVal(0.0);

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      const amrex::Box g0box = amrex::grow(vbox, nGrowD + nGrowC + nGrowT + 1);
      const amrex::Box g1box = amrex::grow(vbox, nGrowC + nGrowT + 1);
      const amrex::Box g2box = amrex::grow(vbox, nGrowD + nGrowC + 1);
      const amrex::Box g3box = amrex::grow(vbox, nGrowC + 1);
      const amrex::Box g4box = amrex::grow(vbox, 1);
      const amrex::Box cbox = amrex::grow(vbox, 0);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(cbox);
      if (typ != amrex::FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
      if (typ == amrex::FabType::covered) {
        continue;
      }

      auto const& s = S.array(mfi);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      amrex::FArrayBox q(g0box, QVAR, amrex::The_Async_Arena());
      amrex::FArrayBox qaux(g0box, nqaux, amrex::The_Async_Arena());
      auto const& q_ar = q.array();
      auto const& qauxar = qaux.array();

      // 1. Get primitives, Q, including (Y, T, p, rho) from conserved state
      // required for L term
      {
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          g0box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, s, q_ar, qauxar);
          });
      }

      // 2. Get dynamic Smagorinsky derived quantities after setting the
      // BC. These quantities need to be stored because we need to filter
      // them at the test filter level. All are located at cell centers.
      const int upper_triangle_n =
        static_cast<int>(0.5 * AMREX_SPACEDIM * (AMREX_SPACEDIM + 1));
      amrex::FArrayBox K(g1box, upper_triangle_n, amrex::The_Async_Arena());
      amrex::FArrayBox RUT(g1box, AMREX_SPACEDIM, amrex::The_Async_Arena());
      amrex::FArrayBox alphaij(
        g1box, AMREX_SPACEDIM * AMREX_SPACEDIM, amrex::The_Async_Arena());
      amrex::FArrayBox alpha(g1box, 1, amrex::The_Async_Arena());
      amrex::FArrayBox flux_T(g1box, AMREX_SPACEDIM, amrex::The_Async_Arena());

      auto const& K_ar = K.array();
      auto const& RUT_ar = RUT.array();
      auto const& alphaij_ar = alphaij.array();
      auto const& alpha_ar = alpha.array();
      auto const& flux_T_ar = flux_T.array();
      {
        const int les_filter_fgr_local = PeleC::les_filter_fgr;
        BL_PROFILE("PeleC::pc_smagorinsky_sfs_term()");
        amrex::ParallelFor(
          g1box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_dynamic_smagorinsky_quantities(
              i, j, k, q_ar, les_filter_fgr_local, dx, K_ar, RUT_ar, alphaij_ar,
              alpha_ar, flux_T_ar);
          });
      }

      // 3. Filter the state variables and the derived quantities at the
      // test filter level - still at cell centers
      amrex::FArrayBox filtered_S(g2box, NVAR, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_Q(g2box, QVAR, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_Qaux(
        g2box, NQAUX > 0 ? NQAUX : 1, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_K(
        g3box, upper_triangle_n, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_RUT(
        g3box, AMREX_SPACEDIM, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_alphaij(
        g3box, AMREX_SPACEDIM * AMREX_SPACEDIM, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_alpha(g3box, 1, amrex::The_Async_Arena());
      amrex::FArrayBox filtered_flux_T(
        g3box, AMREX_SPACEDIM, amrex::The_Async_Arena());

      auto const& filtered_S_ar = filtered_S.array();
      auto const& filtered_Q_ar = filtered_Q.array();
      auto const& filtered_Qaux_ar = filtered_Qaux.array();

      const amrex::FArrayBox& Sfab = S[mfi];
      test_filter.apply_filter(g2box, Sfab, filtered_S);
      {
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          g2box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, filtered_S_ar, filtered_Q_ar, filtered_Qaux_ar);
          });
      }
      test_filter.apply_filter(g3box, K, filtered_K);
      test_filter.apply_filter(g3box, RUT, filtered_RUT);
      test_filter.apply_filter(g3box, alphaij, filtered_alphaij);
      test_filter.apply_filter(g3box, alpha, filtered_alpha);
      test_filter.apply_filter(g3box, flux_T, filtered_flux_T);

      // 4. Calculate the dynamic Smagorinsky coefficients - still at cell
      // centers
      amrex::FArrayBox coeff_cc(g3box, nCompC, amrex::The_Async_Arena());
      auto const& coeff_cc_ar = coeff_cc.array();
      auto const& filtered_K_ar = filtered_K.array();
      auto const& filtered_RUT_ar = filtered_RUT.array();
      auto const& filtered_alphaij_ar = filtered_alphaij.array();
      auto const& filtered_alpha_ar = filtered_alpha.array();
      auto const& filtered_flux_T_ar = filtered_flux_T.array();
      {
        const int les_test_filter_fgr_local = PeleC::les_test_filter_fgr;
        BL_PROFILE("PeleC::pc_dynamic_smagorinsky_coeffs()");
        amrex::ParallelFor(
          g3box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_dynamic_smagorinsky_coeffs(
              i, j, k, filtered_Q_ar, les_test_filter_fgr_local, dx,
              filtered_K_ar, filtered_RUT_ar, filtered_alphaij_ar,
              filtered_alpha_ar, filtered_flux_T_ar, coeff_cc_ar);
          });
      }

      // 5. Filter to smooth the dynamic coefficients - still at cell
      // centers
      auto const& LES_Coeffs_ar = LES_Coeffs[mfi].array();
      coeff_filter.apply_filter(g4box, coeff_cc, LES_Coeffs[mfi]);

      // 6. Get the SFS term

      // First step: move everything needed to compute fluxes to ec (faces)
      const amrex::Array<const amrex::Box, AMREX_SPACEDIM> eboxes = {
        AMREX_D_DECL(
          amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
          amrex::surroundingNodes(cbox, 2))};
      amrex::FArrayBox coeff_ec[AMREX_SPACEDIM];
      amrex::FArrayBox alphaij_ec[AMREX_SPACEDIM];
      amrex::FArrayBox alpha_ec[AMREX_SPACEDIM];
      amrex::FArrayBox flux_T_ec[AMREX_SPACEDIM];
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> coeff_ec_arr;
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
        alphaij_ec_arr;
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> alpha_ec_arr;
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flux_T_ec_arr;

      const bool l_transport_harmonic_mean = transport_harmonic_mean;
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        coeff_ec[dir].resize(eboxes[dir], nCompC, amrex::The_Async_Arena());
        alphaij_ec[dir].resize(
          eboxes[dir], AMREX_SPACEDIM, amrex::The_Async_Arena());
        alpha_ec[dir].resize(eboxes[dir], 1, amrex::The_Async_Arena());
        flux_T_ec[dir].resize(eboxes[dir], 1, amrex::The_Async_Arena());
        coeff_ec_arr[dir] = coeff_ec[dir].array();
        alphaij_ec_arr[dir] = alphaij_ec[dir].array();
        alpha_ec_arr[dir] = alpha_ec[dir].array();
        flux_T_ec_arr[dir] = flux_T_ec[dir].array();
        amrex::ParallelFor(
          eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::Real c[nCompC] = {0.0};
            for (int n = 0; n < nCompC; n++) {
              pc_move_transcoefs_to_ec(
                AMREX_D_DECL(i, j, k), n, LES_Coeffs_ar, c, dir,
                l_transport_harmonic_mean);
              coeff_ec_arr[dir](i, j, k, n) = c[n];
            }

            amrex::Real aijc[AMREX_SPACEDIM * AMREX_SPACEDIM] = {0.0};
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              pc_move_transcoefs_to_ec(
                AMREX_D_DECL(i, j, k), dir * AMREX_SPACEDIM + n, alphaij_ar,
                aijc, dir, l_transport_harmonic_mean);
              alphaij_ec_arr[dir](i, j, k, n) = aijc[dir * AMREX_SPACEDIM + n];
            }

            amrex::Real ac[1] = {0.0};
            pc_move_transcoefs_to_ec(
              AMREX_D_DECL(i, j, k), 0, alpha_ar, ac, dir,
              l_transport_harmonic_mean);
            alpha_ec_arr[dir](i, j, k, 0) = ac[0];

            amrex::Real tc[1] = {0.0};
            pc_move_transcoefs_to_ec(
              AMREX_D_DECL(i, j, k), 0, flux_T_ar, tc, dir,
              l_transport_harmonic_mean);
            flux_T_ec_arr[dir](i, j, k, 0) = tc[0];
          });
      }

      // Compute the fluxes at the faces: all values passed are at faces
      // except for Q, V, and Lterm
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec;
      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        a{{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      resizeAndSetFlux(eboxes, flx, flux_ec);
      {
        BL_PROFILE("PeleC::pc_dynamic_smagorinsky_sfs_term()");
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          amrex::ParallelFor(
            eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_dynamic_smagorinsky_sfs_term(
                i, j, k, q_ar, alphaij_ec_arr[dir], alpha_ec_arr[dir],
                flux_T_ec_arr[dir], coeff_ec_arr[dir], a[dir], dir, flx[dir]);
            });
        }
      }

      // Compute flux divergence (1/Vol).Div(F.A)
      computeFluxDiv(LESTerm, mfi, vbox, flx, volume);

      // Refluxing
      updateFluxRegistersLES(
        reflux_factor, LESTerm, dt, mfi, typ,
        {AMREX_D_DECL(flux_ec.data(), &flux_ec[1], &flux_ec[2])});
    }
  }
#endif
}

// Calculate the LES term using the WALE SFS model
void
PeleC::getWALELESTerm(
#if AMREX_SPACEDIM < 3
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  amrex::MultiFab& /*LESTerm*/,
  amrex::Real /*reflux_factor*/)
{
  amrex::Abort("LES only implemented in 3D for now");
#else
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  amrex::Real reflux_factor)
{
  // Only use this functionality for 3D
  const int ngrow = 1;
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

  amrex::MultiFab S(grids, dmap, NVAR, ngrow, amrex::MFInfo(), Factory());
  FillPatch(*this, S, ngrow, time, State_Type, 0, NVAR);

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      const amrex::Box gbox = amrex::grow(vbox, ngrow);
      const amrex::Box cbox = amrex::grow(vbox, ngrow - 1);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(cbox);
      if (typ != amrex::FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
      if (typ == amrex::FabType::covered) {
        continue;
      }

      auto const& s = S.array(mfi);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      amrex::FArrayBox q(gbox, QVAR, amrex::The_Async_Arena());
      amrex::FArrayBox qaux(gbox, nqaux, amrex::The_Async_Arena());
      auto const& q_ar = q.array();
      auto const& qauxar = qaux.array();

      // Get primitives, Q, including (Y, T, p, rho) from conserved state
      // required for L term
      {
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, s, q_ar, qauxar);
          });
      }

      // Get the tangential derivatives
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> tander_ec;
      const amrex::Array<const amrex::Box, AMREX_SPACEDIM> eboxes = {
        AMREX_D_DECL(
          amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
          amrex::surroundingNodes(cbox, 2))};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> tanders;
      computeTangentialVelDerivs(eboxes, tander_ec, tanders, q_ar, dx);

      // Compute extensive LES fluxes, F.A
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec;

      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        a{{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      resizeAndSetFlux(eboxes, flx, flux_ec);
      {
        BL_PROFILE("PeleC::pc_wale_sfs_term()");
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          amrex::Real Cw_local = PeleC::Cw;
          amrex::Real CI_local = PeleC::CI;
          amrex::Real PrT_local = PeleC::PrT;
          amrex::ParallelFor(
            eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_wale_sfs_term(
                i, j, k, q_ar, tanders[dir], a[dir], dx, dir, Cw_local,
                CI_local, PrT_local, flx[dir]);
            });
        }
      }

      // Compute flux divergence (1/Vol).Div(F.A)
      computeFluxDiv(LESTerm, mfi, vbox, flx, volume);

      // Refluxing
      updateFluxRegistersLES(
        reflux_factor, LESTerm, dt, mfi, typ,
        {AMREX_D_DECL(flux_ec.data(), &flux_ec[1], &flux_ec[2])});
    } // End of MFIter scope
  } // End of OMP scope
#endif
}

// Calculate the LES term using the Vreman SFS model
void
PeleC::getVremanLESTerm(
#if AMREX_SPACEDIM < 3
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  amrex::MultiFab& /*LESTerm*/,
  amrex::Real /*reflux_factor*/)
{
  amrex::Abort("LES only implemented in 3D for now");
#else
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& LESTerm,
  amrex::Real reflux_factor)
{
  // Only use this functionality for 3D
  const int ngrow = 1;
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

  amrex::MultiFab S(grids, dmap, NVAR, ngrow, amrex::MFInfo(), Factory());
  FillPatch(*this, S, ngrow, time, State_Type, 0, NVAR);

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      const amrex::Box gbox = amrex::grow(vbox, ngrow);
      const amrex::Box cbox = amrex::grow(vbox, ngrow - 1);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(cbox);
      if (typ != amrex::FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
      if (typ == amrex::FabType::covered) {
        continue;
      }

      auto const& s = S.array(mfi);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      amrex::FArrayBox q(gbox, QVAR, amrex::The_Async_Arena());
      amrex::FArrayBox qaux(gbox, nqaux, amrex::The_Async_Arena());
      auto const& q_ar = q.array();
      auto const& qauxar = qaux.array();

      // Get primitives, Q, including (Y, T, p, rho) from conserved state
      // required for L term
      {
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, s, q_ar, qauxar);
          });
      }

      // Get the tangential derivatives
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> tander_ec;
      const amrex::Array<const amrex::Box, AMREX_SPACEDIM> eboxes = {
        AMREX_D_DECL(
          amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
          amrex::surroundingNodes(cbox, 2))};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> tanders;

      computeTangentialVelDerivs(eboxes, tander_ec, tanders, q_ar, dx);

      // Compute extensive LES fluxes, F.A
      amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM> flux_ec;
      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        a{{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      resizeAndSetFlux(eboxes, flx, flux_ec);
      {
        BL_PROFILE("PeleC::pc_vreman_sfs_term()");
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          amrex::Real Cs_local = PeleC::Cs;
          amrex::Real CI_local = PeleC::CI;
          amrex::Real PrT_local = PeleC::PrT;
          amrex::ParallelFor(
            eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_vreman_sfs_term(
                i, j, k, q_ar, tanders[dir], a[dir], dx, dir, Cs_local,
                CI_local, PrT_local, flx[dir]);
            });
        }
      }

      // Compute flux divergence (1/Vol).Div(F.A)
      computeFluxDiv(LESTerm, mfi, vbox, flx, volume);

      // Refluxing
      updateFluxRegistersLES(
        reflux_factor, LESTerm, dt, mfi, typ,
        {AMREX_D_DECL(flux_ec.data(), &flux_ec[1], &flux_ec[2])});
    } // End of MFIter scope
  } // End of OMP scope
#endif
}

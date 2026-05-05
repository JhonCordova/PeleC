#include "ODTDiffusion.H"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

#include "ODTThermoBridge.H"
#include "PeleC.H"
#include "PelePhysics.H"

namespace pelec::odtles
{

namespace
{

using ConservativeCell = ODTLineState::ConservativeCell;
constexpr int NComp = static_cast<int>(ODTDiffusion::Component::NumComponents);

amrex::Real
getComponent(const ConservativeCell& c, int comp)
{
  AMREX_ALWAYS_ASSERT(comp >= 0 && comp < NComp);
  if (comp == 0) {
    return c.rho;
  }
  if (comp == 1) {
    return c.rhou;
  }
  if (comp == 2) {
    return c.rhov;
  }
  if (comp == 3) {
    return c.rhow;
  }
  return c.rhoE;
}

void
setComponent(ConservativeCell& c, int comp, amrex::Real val)
{
  AMREX_ALWAYS_ASSERT(comp >= 0 && comp < NComp);
  if (comp == 0) {
    c.rho = val;
  } else if (comp == 1) {
    c.rhou = val;
  } else if (comp == 2) {
    c.rhov = val;
  } else if (comp == 3) {
    c.rhow = val;
  } else {
    c.rhoE = val;
  }
}

void
solveTridiagonalInPlace(
  std::vector<amrex::Real>& a,
  std::vector<amrex::Real>& b,
  std::vector<amrex::Real>& c,
  std::vector<amrex::Real>& rhs,
  std::vector<amrex::Real>& x)
{
  const int n = static_cast<int>(rhs.size());
  AMREX_ALWAYS_ASSERT(n > 0);
  AMREX_ALWAYS_ASSERT(
    static_cast<int>(a.size()) == n && static_cast<int>(b.size()) == n &&
    static_cast<int>(c.size()) == n && static_cast<int>(x.size()) == n);

  for (int i = 1; i < n; ++i) {
    AMREX_ALWAYS_ASSERT(b[i - 1] != 0.0);
    const amrex::Real m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    rhs[i] -= m * rhs[i - 1];
  }

  AMREX_ALWAYS_ASSERT(b[n - 1] != 0.0);
  x[n - 1] = rhs[n - 1] / b[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    AMREX_ALWAYS_ASSERT(b[i] != 0.0);
    x[i] = (rhs[i] - c[i] * x[i + 1]) / b[i];
  }
}

amrex::Real
maxAbsStateChange(const ODTLineState& a, const ODTLineState& b)
{
  AMREX_ALWAYS_ASSERT(a.numCells() == b.numCells());
  amrex::Real max_abs = 0.0;
  for (int i = 0; i < a.numCells(); ++i) {
    for (int comp = 0; comp < NComp; ++comp) {
      const amrex::Real da = getComponent(a.cell(i), comp);
      const amrex::Real db = getComponent(b.cell(i), comp);
      max_abs = std::max(max_abs, std::abs(db - da));
    }
  }
  return max_abs;
}

bool
isMomentumComponent(int comp)
{
  return (
    comp == static_cast<int>(ODTDiffusion::Component::RhoU) ||
    comp == static_cast<int>(ODTDiffusion::Component::RhoV) ||
    comp == static_cast<int>(ODTDiffusion::Component::RhoW));
}

amrex::Real
computeDynamicViscosity(
  const ODTThermoBridge::RecoveredState& rec)
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
    get_xi, get_mu, get_lam, get_Ddiag, get_chi, rec.temperature, rec.rho, Y.data(),
    nullptr, nullptr, mu, xi_dummy, lam_dummy, tparm);
  return mu;
}

amrex::Real
lineSegmentSpacing(const ODTLineGeometry& geom)
{
  const auto& segs = geom.supportCells();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    !segs.empty(), "ODTDiffusion requires non-empty line support");

  amrex::Real sum_w = 0.0;
  amrex::Real min_w = std::numeric_limits<amrex::Real>::max();
  amrex::Real max_w = 0.0;
  for (const auto& seg : segs) {
    const amrex::Real w = seg.s_interval.hi - seg.s_interval.lo;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      w > 0.0,
      "ODTDiffusion requires positive physical widths on support segments");
    sum_w += w;
    min_w = std::min(min_w, w);
    max_w = std::max(max_w, w);
  }

  const amrex::Real dx = sum_w / static_cast<amrex::Real>(segs.size());
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    dx > 0.0, "ODTDiffusion requires positive average segment spacing");

  const amrex::Real tol = 1.0e-12 * (std::abs(dx) + 1.0);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    (max_w - min_w) <= tol,
    "ODTDiffusion expects near-uniform support-segment spacing for the "
    "current implicit tridiagonal operator");

  return dx;
}

} // namespace

ODTDiffusion::MomentumViscosityProfile
ODTDiffusion::buildMomentumViscosityProfile(const ODTLineState& line_state)
{
  MomentumViscosityProfile profile{};
  profile.dynamic_viscosity_mu.assign(
    static_cast<std::size_t>(line_state.numCells()), 0.0);

  for (int i = 0; i < line_state.numCells(); ++i) {
    const auto rec = ODTThermoBridge::recoverCellThermoState(line_state, i);
    if (!rec.success) {
      profile.success = false;
      profile.failure_cell = i;
      return profile;
    }
    const amrex::Real mu = computeDynamicViscosity(rec);
    const amrex::Real mu_nonneg = std::max<amrex::Real>(mu, 0.0);
    profile.dynamic_viscosity_mu[static_cast<std::size_t>(i)] = mu_nonneg;
  }

  profile.success = true;
  profile.failure_cell = -1;
  return profile;
}

void
ODTDiffusion::applyImplicitUniform(
  const ODTLineGeometry& geom,
  ODTLineState& line_state,
  amrex::Real dt,
  const Controls& controls)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTDiffusion requires defined line geometry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTDiffusion requires an initialized and valid line state");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.numCells() == geom.supportCellCount(),
    "ODTDiffusion requires state/geometry support-size consistency");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    dt >= 0.0, "ODTDiffusion requires non-negative dt");

  const int n = line_state.numCells();
  if (n <= 1 || dt == 0.0) {
    return;
  }

  // Use the physical support-segment spacing (host-cell width / nsub) so
  // subsegment refinement is represented consistently in the local operator.
  const amrex::Real dx = lineSegmentSpacing(geom);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    dx > 0.0, "ODTDiffusion requires positive support-segment spacing");
  const amrex::Real inv_dx2 = 1.0 / (dx * dx);

  std::vector<amrex::Real> q_old(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> a(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> b(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> c(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> rhs(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> q_new(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> coeff_face(static_cast<std::size_t>(n - 1), 0.0);

  MomentumViscosityProfile mu_profile{};
  bool built_mu_profile = false;

  for (int comp = 0; comp < NComp; ++comp) {
    if (!controls.diffuse_component[static_cast<std::size_t>(comp)]) {
      continue;
    }

    if (
      controls.use_molecular_viscosity_for_momentum &&
      isMomentumComponent(comp)) {
      if (!built_mu_profile) {
        mu_profile = buildMomentumViscosityProfile(line_state);
        if (
          !mu_profile.success &&
          controls.fail_on_molecular_viscosity_recovery_failure) {
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            false,
            "ODTDiffusion molecular-viscosity bridge failed thermochemical "
            "recovery on line cell");
        }
        built_mu_profile = true;
      }
      if (mu_profile.success) {
        const amrex::Real mu_floor =
          std::max<amrex::Real>(controls.molecular_viscosity_floor, 0.0);
        for (int f = 0; f < n - 1; ++f) {
          const amrex::Real mu_l = std::max<amrex::Real>(
            mu_profile.dynamic_viscosity_mu[static_cast<std::size_t>(f)], mu_floor);
          const amrex::Real mu_r = std::max<amrex::Real>(
            mu_profile.dynamic_viscosity_mu[static_cast<std::size_t>(f + 1)],
            mu_floor);
          const amrex::Real mu_face = std::max<amrex::Real>(0.0, 0.5 * (mu_l + mu_r));
          if (
            controls.momentum_molecular_form ==
            MomentumMolecularForm::ConservativeVariableNuOnRhoU) {
            const amrex::Real rho_l = line_state.cell(f).rho;
            const amrex::Real rho_r = line_state.cell(f + 1).rho;
            const amrex::Real nu_l = mu_l / rho_l;
            const amrex::Real nu_r = mu_r / rho_r;
            coeff_face[static_cast<std::size_t>(f)] =
              std::max<amrex::Real>(0.0, 0.5 * (nu_l + nu_r));
          } else {
            coeff_face[static_cast<std::size_t>(f)] = mu_face;
          }
        }
      } else {
        const amrex::Real coeff_const =
          controls.diffusivity[static_cast<std::size_t>(comp)];
        if (coeff_const <= 0.0) {
          continue;
        }
        for (int f = 0; f < n - 1; ++f) {
          coeff_face[static_cast<std::size_t>(f)] = coeff_const;
        }
      }
    } else {
      const amrex::Real coeff_const =
        controls.diffusivity[static_cast<std::size_t>(comp)];
      if (coeff_const <= 0.0) {
        continue;
      }
      for (int f = 0; f < n - 1; ++f) {
        coeff_face[static_cast<std::size_t>(f)] = coeff_const;
      }
    }

    for (int i = 0; i < n; ++i) {
      q_old[static_cast<std::size_t>(i)] = getComponent(line_state.cell(i), comp);
      rhs[static_cast<std::size_t>(i)] = q_old[static_cast<std::size_t>(i)];
    }

    if (
      controls.use_molecular_viscosity_for_momentum && isMomentumComponent(comp) &&
      mu_profile.success &&
      controls.momentum_molecular_form ==
        MomentumMolecularForm::VelocityGradientMuFlux) {
      // Semi-implicit linearization of:
      //   d(rho*u)/dt = d/ds(mu * d(u)/ds), with u = (rho*u)/rho.
      // Here rho and mu are frozen from current line state.
      for (int i = 0; i < n; ++i) {
        const amrex::Real rho_i = line_state.cell(i).rho;
        const amrex::Real mu_l =
          (i > 0) ? coeff_face[static_cast<std::size_t>(i - 1)] : 0.0;
        const amrex::Real mu_r =
          (i < n - 1) ? coeff_face[static_cast<std::size_t>(i)] : 0.0;
        const amrex::Real k = dt * inv_dx2;

        a[static_cast<std::size_t>(i)] =
          (i > 0) ? -k * (mu_l / line_state.cell(i - 1).rho) : 0.0;
        b[static_cast<std::size_t>(i)] = 1.0 + k * ((mu_l + mu_r) / rho_i);
        c[static_cast<std::size_t>(i)] =
          (i < n - 1) ? -k * (mu_r / line_state.cell(i + 1).rho) : 0.0;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        const amrex::Real coeff_l =
          (i > 0) ? coeff_face[static_cast<std::size_t>(i - 1)] : 0.0;
        const amrex::Real coeff_r =
          (i < n - 1) ? coeff_face[static_cast<std::size_t>(i)] : 0.0;
        const amrex::Real r_l = dt * coeff_l * inv_dx2;
        const amrex::Real r_r = dt * coeff_r * inv_dx2;

        // No-flux (Neumann) closure at line ends through zero boundary-face
        // coefficient and interior face diffusion only.
        a[static_cast<std::size_t>(i)] = -r_l;
        b[static_cast<std::size_t>(i)] = 1.0 + r_l + r_r;
        c[static_cast<std::size_t>(i)] = -r_r;
      }
    }

    solveTridiagonalInPlace(a, b, c, rhs, q_new);

    for (int i = 0; i < n; ++i) {
      auto cell = line_state.cell(i);
      setComponent(cell, comp, q_new[static_cast<std::size_t>(i)]);
      line_state.setCell(i, cell);
    }
  }
}

ODTDiffusion::ValidationReport
ODTDiffusion::runMVPValidationHook()
{
  ValidationReport report{};

  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(2, 0, 0)));
  amrex::RealBox rb(
    AMREX_D_DECL(0.0, 0.0, 0.0), AMREX_D_DECL(3.0, 1.0, 1.0));
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(1, 0, 0));

  ODTLineGeometry line_geom(geom, 0, owner, 0);

  Controls ctrl{};
  ctrl.diffuse_component.fill(false);
  ctrl.diffuse_component[static_cast<int>(Component::RhoU)] = true;
  ctrl.use_molecular_viscosity_for_momentum = false;
  ctrl.diffusivity.fill(0.0);
  ctrl.diffusivity[static_cast<int>(Component::RhoU)] = 0.1;

  constexpr amrex::Real dt = 1.0e-2;
  constexpr amrex::Real tol = 1.0e-12;

  // Uniform-profile invariance check.
  ODTLineState uniform_state;
  uniform_state.initialize(line_geom);
  for (int i = 0; i < uniform_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0;
    c.rhou = 2.5;
    c.rhov = -0.4;
    c.rhow = 0.1;
    c.rhoE = 10.0;
    uniform_state.setCell(i, c);
  }
  uniform_state.setValid(true);
  ODTLineState uniform_before = uniform_state;
  applyImplicitUniform(line_geom, uniform_state, dt, ctrl);
  report.uniform_max_abs_change = maxAbsStateChange(uniform_before, uniform_state);
  report.uniform_invariant = (report.uniform_max_abs_change <= tol);

  // Nonuniform-profile evolution check.
  ODTLineState nonuniform_state;
  nonuniform_state.initialize(line_geom);
  for (int i = 0; i < nonuniform_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0;
    c.rhou = static_cast<amrex::Real>(i + 1);
    c.rhov = 0.0;
    c.rhow = 0.0;
    c.rhoE = 5.0;
    nonuniform_state.setCell(i, c);
  }
  nonuniform_state.setValid(true);
  ODTLineState nonuniform_before = nonuniform_state;
  applyImplicitUniform(line_geom, nonuniform_state, dt, ctrl);
  report.nonuniform_max_abs_change =
    maxAbsStateChange(nonuniform_before, nonuniform_state);
  report.nonuniform_evolved = (report.nonuniform_max_abs_change > tol);

  return report;
}

} // namespace pelec::odtles

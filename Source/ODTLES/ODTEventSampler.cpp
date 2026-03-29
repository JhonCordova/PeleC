#include "ODTEventSampler.H"

#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX.H>

namespace pelec::odtles
{

ODTEventSampler::ODTEventSampler(const Controls& controls)
{
  setControls(controls);
}

ODTEventSampler::ODTEventSampler() : ODTEventSampler(Controls{}) {}

void
ODTEventSampler::setControls(const Controls& controls)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    controls.event_rate >= 0.0,
    "ODTEventSampler requires non-negative event_rate");
  m_ctrl = controls;

  if (m_ctrl.deterministic) {
    m_rng.seed(static_cast<std::mt19937_64::result_type>(m_ctrl.seed));
  } else {
    std::random_device rd;
    const auto s0 = static_cast<std::uint64_t>(rd());
    const auto s1 = static_cast<std::uint64_t>(rd());
    m_rng.seed(static_cast<std::mt19937_64::result_type>((s1 << 32U) ^ s0));
  }
}

void
ODTEventSampler::reseed(std::uint64_t seed)
{
  m_ctrl.seed = seed;
  m_rng.seed(static_cast<std::mt19937_64::result_type>(seed));
}

ODTEventSampler::EventSample
ODTEventSampler::sample(int line_cells)
{
  EventSample out{};

  if (line_cells <= 0 || m_ctrl.event_rate <= 0.0) {
    out.wait_time = std::numeric_limits<amrex::Real>::infinity();
    out.valid = false;
    return out;
  }

  std::exponential_distribution<amrex::Real> wait_dist(m_ctrl.event_rate);
  out.wait_time = wait_dist(m_rng);

  const int min_n = std::max(1, m_ctrl.min_interval_size);
  const int max_n = line_cells;
  if (min_n > max_n) {
    out.valid = false;
    return out;
  }

  std::uniform_int_distribution<int> len_dist(min_n, max_n);
  const int n = len_dist(m_rng);

  std::uniform_int_distribution<int> lo_dist(0, line_cells - n);
  const int i_lo = lo_dist(m_rng);

  out.interval.i_lo = i_lo;
  out.interval.i_hi = i_lo + n - 1;
  out.valid = true;
  return out;
}

} // namespace pelec::odtles

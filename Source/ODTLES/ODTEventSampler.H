#ifndef PELEC_ODTLES_ODTEVENTSAMPLER_H
#define PELEC_ODTLES_ODTEVENTSAMPLER_H

#include <limits>
#include <cstdint>
#include <random>

#include <AMReX_REAL.H>

#include "ODTTripletMap.H"

namespace pelec::odtles
{

class ODTEventSampler
{
public:
  struct Controls
  {
    // Poisson candidate sampling rate for internal waiting-time draws.
    amrex::Real event_rate = 0.0;
    bool deterministic = true;
    std::uint64_t seed = 1ULL;

    // MVP size/position sampling window in local support cells.
    int min_interval_size = 2;
  };

  struct EventSample
  {
    amrex::Real wait_time = std::numeric_limits<amrex::Real>::infinity();
    ODTTripletMap::EddyInterval interval{};
    bool valid = false;
  };

  ODTEventSampler();
  explicit ODTEventSampler(const Controls& controls);

  void setControls(const Controls& controls);
  const Controls& controls() const noexcept { return m_ctrl; }
  void reseed(std::uint64_t seed);
  amrex::Real sampleUnitUniform();
  amrex::Real sampleExponentialWait(amrex::Real rate);

  // Sample waiting time + interval over the local line support.
  EventSample sample(int line_cells);

private:
  Controls m_ctrl{};
  std::mt19937_64 m_rng;
};

} // namespace pelec::odtles

#endif

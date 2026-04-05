#pragma once

#include <tuple>
#include <utility>
#include <valarray>

#include "global_def.hh"
#include "ppm/ppm_constant.hh"
#include "vector_wrapper.hh"

namespace ppm
{

constexpr double xwig = 0.0;

constexpr double omega1 = 0.75;
constexpr double omega2 = 5.0;
constexpr double epsilon = 0.33;

constexpr index_type N_STATEV = 3;
constexpr index_type NZ_GHOST = 6;

// constexpr size_t GRID_NX = 100;
// constexpr size_t GRID_NY = 100;

template <index_type N>
using ppmlr_vf_t = vector_frame_t<N, NZ_GHOST>;

template <index_type N>
using ppmlr_vector_t = vector_wrapper_t<N, NZ_GHOST>;

template <index_type N>
struct ppmlr_base_t {
  constexpr static index_type vlen = N;

  using vector_type = ppmlr_vector_t<N>;
  using bp = ppmlr_vf_t<N>;
};


}  // namespace ppm

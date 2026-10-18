#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include "global_def.hpp"
#include "ppmlr/ppm_constant.hpp"

namespace ppm
{

constexpr double xwig = 0.0;

constexpr double omega1 = 0.75;
constexpr double omega2 = 5.0;
constexpr double epsilon = 0.33;

constexpr index_type N_STATEV = 3;
constexpr index_type NZ_GHOST = 6;

// Frame concept for ghost cell indexing
template <index_type N, index_type B>
struct vector_frame_t {
  constexpr static index_type i0 = 0;
  constexpr static index_type iN = N + 2 * B;
  constexpr static index_type j0 = B;
  constexpr static index_type jM = N + B;
};

template <index_type N>
using ppmlr_vf_t = vector_frame_t<N, NZ_GHOST>;

template <index_type N>
struct ppmlr_base_t {
  constexpr static index_type vlen = N;

  using vector_type = std::vector<double>;
  using bp = ppmlr_vf_t<N>;

  // Helper to create properly sized vector with ghost cells
  static vector_type make_vector() {
    return vector_type(bp::iN, 0.0);
  }
};


}  // namespace ppm

#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>

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

  // Kokkos view for device offloading
  // Uses default execution space memory (device for CUDA/HIP, host for Serial/OpenMP)
  using vector_type = Kokkos::View<double*>;
  using bp = ppmlr_vf_t<N>;

  // Helper to create properly sized vector with ghost cells
  // Label string aids debugging with Kokkos profiling tools
  static vector_type make_vector() {
    return vector_type("vector", bp::iN);
  }

  // Helper to create named vector (for better profiling visibility)
  static vector_type make_vector(const char* label) {
    return vector_type(label, bp::iN);
  }
};


}  // namespace ppm

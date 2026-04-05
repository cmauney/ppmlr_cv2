#pragma once

#include <cstddef>
#include <iostream>

//#include "hydro_physics.hh"
#include "geometry.hh"
#include "ppm_constant.hh"
#include "ppm_def.hh"
#include "ppm_vector.hh"

namespace ppm
{

template <typename GeometorType, index_type N>
struct ppm_grid_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using geo_type = GeometorType;

  // grid values
  std::array<vector_type, 2> xe, xc, dx, dvol;

  // storage space for wiggle
  vector_type __xe, __xc, __dx, __dvol;

  vector_type amid, xnolap, xndiff;

  //  sweep_vector_t<NZ_SWEEP> r, u, v, w, p, e, f;
  //  sweep_vector_t<NZ_SWEEP> q;

  // geometry variable
  geo_type geometer;
  double radius;

  ppm_grid_t() : radius(1.0) {}

  template <index_type I>
  constexpr inline void build_volume() {
    for (index_type i = bp::i0; i < bp::iN; ++i) {
      dvol[I][i] = geometer.volume(xe[I][i], dx[I][i]);
    }
  }

  constexpr inline void build_amid() {
    for (index_type i = bp::i0; i < bp::iN; ++i) {
      amid[i] = geometer.area_mid(xe[0][i], xe[1][i]);
    }
  }

  constexpr inline void build_overlap() {
    for (index_type i = bp::i0; i < bp::iN; ++i) {
      xnolap[i] = geometer.overlap_volume(xe[0][i], xe[1][i]);
    }

    for (index_type i = bp::j0; i < bp::jM + 1; ++i) {
      xndiff[i] = xe[1][i] - xe[0][i];
    }
  }

  inline void store_coords() {
    for (index_type i = bp::i0; i < bp::iN; ++i) {
      __xe[i] = xe[0][i];
      __xc[i] = xc[0][i];
      __dx[i] = dx[0][i];
      __dvol[i] = dvol[0][i];
    }
    //  __xe[vector_type::iN] = xe[0][vector_type::iN];
  }

  inline void retieve_coords() {
    for (index_type i = bp::i0; i < bp::iN; ++i) {
      xe[0][i] = __xe[i];
      xc[0][i] = __xc[i];
      dx[0][i] = __dx[i];
      dvol[0][i] = __dvol[i];
    }
    //    xe[0][vector_type::iN] = __xe[vector_type::iN];
  }
};

} // namespace ppm

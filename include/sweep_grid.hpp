#pragma once

#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>

#include "hydro_physics.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "sweep_vector.hpp"

namespace ppm
{
namespace geo {
struct planar_geo_t {
  constexpr planar_geo_t() {}
  constexpr inline double volume(const double x, const double dx) const {
    return dx;
  }
  constexpr inline double area_mid(const double x0, const double xn) const {
    return 1.0;
  }
  constexpr inline double overlap_volume(const double x0,
                                         const double xn) const {
    return xn - x0;
  }
};

struct cylin_geo_t {
  constexpr cylin_geo_t() {}
  constexpr inline double volume(const double x, const double dx) const {
    return dx * (x + half * dx);
  }
  constexpr inline double area_mid(const double x0, const double xn) const {
    return half * (x0 + xn);
  }
  constexpr inline double overlap_volume(const double x0,
                                         const double xn) const {
    auto dx = xn - x0;

    return dx * (x0 + half * dx);
  }
};

}  // namespace geo

template <typename PPMVector, typename GeometorType>
struct ppm_grid_t {
  using vector_type = PPMVector;
  using index_type = typename vector_type::index_type;
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
    for (auto i = 0; i < vector_type::i0; ++i) {
      dvol[I][i] = geometer.volume(xe[I][i], dx[I][i]);
    }
  }

  constexpr inline void build_amid() {
    for (auto i = 0; i < vector_type::i0; ++i) {
      amid[i] = geometer.overlap_area(xe[0][i], xe[1][i]);
    }
  }

  constexpr inline void build_overlap() {
    for (auto i = 0; i < vector_type::i0; ++i) {
      xnolap[i] = geometer.overlap_volume(xe[0][i], xe[0][i]);
    }

    for (auto i = vector_type::j0; i < vector_type::jM + 1; ++i) {
      xndiff[i] = xe[1][i] - xe[0][i];
    }
  }

  inline void store_coords() {
    for (index_type i = vector_type::i0(); i < vector_type::iN(); ++i) {
      __xe[i] = xe[0][i];
      __xc[i] = xc[0][i];
      __dx[i] = dx[0][i];
      __dvol[i] = dvol[0][i];
    }
    //  __xe[vector_type::iN()] = xe[0][vector_type::iN()];
  }

  inline void retieve_coords() {
    for (index_type i = vector_type::i0; i < vector_type::iN; ++i) {
      xe[0][i] = __xe[i];
      xc[0][i] = __xc[i];
      dx[0][i] = __dx[i];
      dvol[0][i] = __dvol[i];
    }
    //    xe[0][vector_type::iN()] = __xe[vector_type::iN()];
  }
};

} // namespace ppm

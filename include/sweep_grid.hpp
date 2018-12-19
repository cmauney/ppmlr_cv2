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

template <typename SweepVector>
struct sweep_grid_t {
  using vector_type = SweepVector;
  using index_type = typename SweepVector::index_type;

  // grid values
  std::array<vector_type, 2> xe, xc, dx, dvol;

  // storage space for wiggle
  vector_type __xe, __xc, __dx, __dvol;

  vector_type amid, xnolap, xndiff;

  //  sweep_vector_t<NZ_SWEEP> r, u, v, w, p, e, f;
  //  sweep_vector_t<NZ_SWEEP> q;

  // geometry variable
  GEO geom;
  double radius;

  sweep_grid_t()
      : geom(GEO::CART),
        radius(1.0) {}

  inline void build_volume(const index_type ti) {
    switch (geom) {
      case GEO::CART:
        for (index_type i = vector_type::i0(); i < vector_type::iN(); ++i)
          dvol[ti][i] = dx[ti][i];
        break;
      default:
        std::cerr << "not a geo i kno\n";
        std::exit(1);
    }
  }

  inline void build_amid() {
    switch (geom) {
      case GEO::CART:
        for (index_type i = vector_type::i0(); i < vector_type::iN(); ++i)
          amid[i] = 1.0;
        break;
      case GEO::CYLIN:
        for (index_type i = vector_type::i0(); i < vector_type::iN(); ++i)
          amid[i] = 0.5 * (xe[1][i] + xe[0][i]);
        break;
      default:
        std::cerr << "not a geo i kno\n";
        std::exit(1);
    }
  }

  inline void build_overlap() {
    switch (geom) {
      case GEO::CART:
        for (index_type i = vector_type::j0(); i < vector_type::jM() + 1; ++i) {
          xnolap[i] = xe[1][i] - xe[0][i];
        }
        break;
      default:
        std::cerr << "not a geo i kno\n";
        std::exit(1);
    }
    for (index_type i = vector_type::j0(); i < vector_type::jM() + 1; ++i) {
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
  __xe[vector_type::iN()] = xe[0][vector_type::iN()];
  }

  inline void retieve_coords() {
    for (index_type i = vector_type::i0(); i < vector_type::iN(); ++i) {
      xe[0][i] = __xe[i];
      xc[0][i] = __xc[i];
      dx[0][i] = __dx[i];
      dvol[0][i] = __dvol[i];
    }
    xe[0][vector_type::iN()] = __xe[vector_type::iN()];
  }
};

} // namespace ppm

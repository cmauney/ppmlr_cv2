#pragma once

#include <cmath>
#include <iostream>
#include <valarray>

#include "global_const.hpp"
#include "ppmlr/ppm_constant.hpp"
#include "vector_wrapper.hpp"

namespace hydro {

constexpr double courant = 0.5;

struct hydro_physics_t {
  const double gamma, gammam1;
  double last_dt, dt, hdt;

  double svel;

  hydro_physics_t(double gam)
      : gamma(gam), gammam1(gam - 1), last_dt(ppm::small), svel(zero) {}

  template <index_type N>
  constexpr inline double ridt(const double svel,
                               const vector_wrapper_t<N>& xvel,
                               const vector_wrapper_t<N>& dx) {
    vector_wrapper_t<N> vm, vabs;

    vabs = xvel.abs();
    for (index_type i = 0; i < N; ++i) {
      vm[i] = vabs[i] / dx[i];
    }

    return std::max(svel, (vm.max()));
  }

  constexpr inline void set_dt(double timestep) {
    last_dt = dt;

    dt = timestep;
    hdt = half * dt;
  }

  template <index_type N>
  constexpr inline void first_dt(const vector_wrapper_t<N>& prs,
                                 const vector_wrapper_t<N>& rho,
                                 const vector_wrapper_t<N>& xvel,
                                 const vector_wrapper_t<N>& dx) {
    vector_wrapper_t<N> ctdx;
    for (index_type i = 0; i < N; ++i)
      ctdx[i] = std::sqrt(gamma * prs[i] / rho[i]) / dx[i];

    auto svel = ctdx.max();

    dt = courant / ridt(svel, xvel, dx);
    hdt = half * dt;
    last_dt = dt;
  }

  template <index_type N>
  constexpr inline void dt_courant(double cdx,
                                   const vector_wrapper_t<N>& xvel,
                                   const vector_wrapper_t<N>& dx) {
    auto dt_n = courant / ridt(cdx, xvel, dx);

    auto dtu = last_dt * 1.1;
    set_dt(std::min(dtu, dt_n));
  }
};

}  // namespace hydro

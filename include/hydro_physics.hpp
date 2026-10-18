#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "global_const.hpp"
#include "global_def.hpp"
#include "ppmlr/ppm_constant.hpp"

namespace hydro {

constexpr double courant = 0.5;

struct hydro_physics_t {
  const double gamma, gammam1;
  double last_dt, dt, hdt;

  double svel;

  hydro_physics_t(double gam)
      : gamma(gam), gammam1(gam - 1), last_dt(ppm::small), svel(zero) {}

  inline double ridt(const double svel,
                     const std::vector<double>& xvel,
                     const std::vector<double>& dx) {
    // Modern C++20: use std::transform_reduce for max calculation
    const size_t N = xvel.size();
    double max_vm = 0.0;
    for (size_t i = 0; i < N; ++i) {
      max_vm = std::max(max_vm, std::abs(xvel[i]) / dx[i]);
    }
    return std::max(svel, max_vm);
  }

  constexpr inline void set_dt(double timestep) {
    last_dt = dt;

    dt = timestep;
    hdt = half * dt;
  }

  inline void first_dt(const std::vector<double>& prs,
                       const std::vector<double>& rho,
                       const std::vector<double>& xvel,
                       const std::vector<double>& dx) {
    // Modern C++20: compute max sound speed
    const size_t N = prs.size();
    double svel = 0.0;
    for (size_t i = 0; i < N; ++i) {
      svel = std::max(svel, std::sqrt(gamma * prs[i] / rho[i]) / dx[i]);
    }

    dt = courant / ridt(svel, xvel, dx);
    hdt = half * dt;
    last_dt = dt;
  }

  inline void dt_courant(double cdx,
                         const std::vector<double>& xvel,
                         const std::vector<double>& dx) {
    auto dt_n = courant / ridt(cdx, xvel, dx);

    auto dtu = last_dt * 1.1;
    set_dt(std::min(dtu, dt_n));
  }
};

}  // namespace hydro

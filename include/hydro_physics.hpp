#pragma once

#include <cmath>
#include <iostream>
#include <valarray>

#include "ppm_constant.hpp"
#include "ppm_def.hpp"

namespace ppm
{

    struct hydro_physics_t
    {
      const double gamma, gammam1;
      double last_dt, dt, hdt;

      double svel;

      hydro_physics_t(double gam)
          : gamma(gam), gammam1(gam - 1), last_dt(small), svel(zero) {}

      template <typename V>
      constexpr inline double ridt(const double svel,
                                   const V& xvel,
                                   const V& dx) {
        return std::max(svel, (std::abs(xvel) / dx).max());
      }

      constexpr inline void set_dt(double timestep) {
        last_dt = dt;

        dt = timestep;
        hdt = 0.5 * dt;
      }

      template <typename V>
      constexpr inline void first_dt(const V& prs,
                                     const V& rho,
                                     const V& xvel,
                                     const V& dx) {
        auto svel = (std::sqrt(gamma * prs / rho) / dx).max();
        set_dt(courant / ridt(svel, xvel, dx));
      }

      template <typename V>
      constexpr inline void dt_courant(double cdx, const V& xvel, const V& dx) {
        auto dt_n = courant / ridt(cdx, xvel, dx);

        auto dtu = last_dt * 1.1;
        set_dt(std::min(dtu, dt_n));
      }
    };

} // namespace ppm

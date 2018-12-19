#pragma once

#include <cmath>
#include <iostream>
#include <valarray>

#include "ppm_def.hpp"

namespace ppm
{

    struct hydro_physics_t
    {
      const double gamma, gammam1;
      double last_dt, dt, hdt;

      double svel;

      hydro_physics_t(double gam) : gamma(gam), gammam1(gam - 1), svel(0.0) {}

      constexpr inline void set_dt(double timestep) {
        auto dtu = last_dt * 1.1;
        if (dtu == 0.0)
          dt = timestep;
        else
          dt = std::min(dtu, timestep);

        hdt = 0.5 * dt;
        last_dt = dt;
      }

      template <typename V>
      constexpr inline void dt_courant(double cdx, const V& xvel, const V& dx) {
        double ridt = 0.0;

        ridt = (std::abs(xvel) / dx).max();
        ridt = std::max(cdx, ridt);
        auto dt_n = courant / ridt;
        set_dt(dt_n);
      }
    };

} // namespace ppm

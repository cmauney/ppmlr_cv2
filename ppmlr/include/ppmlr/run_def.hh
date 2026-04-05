#pragma once

#include "ppm/ppm_grid.hh"

#include "geometry.hh"
#include "global_const.hh"
#include "global_def.hh"

namespace hydro {
constexpr size_t NDIM = 1;

// X dim definitions
using XGEO = geo::planar_geo_t;
constexpr size_t NX = 250;

}  // namespace hydro

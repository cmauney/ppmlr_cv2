#pragma once

#include "ppmlr/ppm_grid.hpp"

#include "geometry.hpp"
#include "global_const.hpp"
#include "global_def.hpp"

namespace hydro
{
constexpr size_t NDIM = 1;

// X dim definitions
using XGEO = geo::planar_geo_t;
constexpr size_t NX = 250;

} // namespace hydro

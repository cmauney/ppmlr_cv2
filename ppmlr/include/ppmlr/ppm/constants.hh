#pragma once

#include <cstddef>

namespace ppmlr::ppm {

constexpr std::size_t N_STATE_VAR = 3;
constexpr std::size_t N_GHOST = 6;

constexpr double small = 1.0E-15;
constexpr double smallp = 1.0E-12;
constexpr double smallr = 1.0E-12;
constexpr double small_riemann = 1.0E-25;

}  // namespace ppmlr::ppm

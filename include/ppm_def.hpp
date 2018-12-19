#pragma once

#include <tuple>
#include <utility>
#include <valarray>

namespace ppm
{

enum GEO { CART = 0, CYLIN, SPHERI, CYLIN_THETA, SPHERI_THETA, SPHERI_PHI };

enum FLUID_VI { RHO, PRS, XVL, YVL, ZVL, ENT, EIN, FLT };

constexpr std::array<FLUID_VI, 3> PRIMATIVE_VARS{RHO, PRS, XVL};

constexpr std::array<FLUID_VI, 6> CONSERVATIVE_VARS{RHO, XVL, YVL,
                                                    ZVL, EIN, ENT};

constexpr size_t N_FLUIDV = 8;
constexpr size_t N_STATEV = 3;

constexpr size_t NZ_GHOST = 6;

constexpr size_t GRID_NX = 100;
constexpr size_t GRID_NY = 100;

constexpr double xwig = 0.0;
constexpr double courant = 0.5;

constexpr double omega1 = 0.75;
constexpr double omega2 = 5.0;
constexpr double epsilon = 0.33;

}  // namespace ppm

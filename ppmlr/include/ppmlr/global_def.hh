#pragma once

#include <array>
#include <utility>

// index type
using index_type = size_t;

// vector storage type
// template <index_type N>
// using fixed_storage_type = std::array<double, N>;

// fluid variable IDs
enum FLUID_VI { RHO, PRS, XVL, YVL, ZVL, ENT, EIN, FLT };

// IDs of the primative variables (technically prs is not primative but easier
// to use. these define the states of the riemann problem during ppm
constexpr std::array<FLUID_VI, 3> PRIMATIVE_VARS{RHO, PRS, XVL};

// conserved quantities. note these are 'specific', i.e. per-unit-mass.
constexpr std::array<FLUID_VI, 6> CONSERVATIVE_VARS{RHO, XVL, YVL,
                                                    ZVL, EIN, ENT};

// total number of advected quantities
constexpr index_type N_FLUIDV = 8;

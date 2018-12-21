#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "gsl/gsl"

#include "hydro_physics.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "sweep_grid.hpp"
#include "sweep_vector.hpp"

namespace ppm {

template <typename PPMVector, typename Grid>
struct boundary_conditions_t {
  using vector_type = PPMVector;
  using index_type = typename vector_type::index_type;
  using grid_type = Grid;

  boundary_conditions_t() = default;
  ~boundary_conditions_t() = default;

  void reflect_bc(vector_type& v, double symm = 1.0) {
    for (index_type i = vector_type::i0(); i < vector_type::j0(); ++i) {
      v[vector_type::j0() - i - 1] = symm * v[i + vector_type::j0()];
    }

    for (index_type i = vector_type::i0(); i < vector_type::j0(); ++i) {
      v[vector_type::jM() + i] = symm * v[vector_type::jM() - i - 1];
    }
  }

  inline void apply(grid_type& g, FluidVars<PPMVector>& vars) {
    reflect_bc(g.dx[0]);
    reflect_bc(vars[RHO]);
    reflect_bc(vars[PRS]);
    reflect_bc(vars[XVL], -1.0);
    reflect_bc(vars[YVL]);
    reflect_bc(vars[ZVL]);
    reflect_bc(vars[ENT]);
    reflect_bc(vars[FLT]);

    for (size_t i = 0; i < vector_type::j0(); ++i) {
      g.xe[0][vector_type::j0() - i - 1] =
          g.xe[0][vector_type::j0() - i] - g.dx[0][vector_type::j0() - i - 1];

      auto ri = vector_type::jM();
      g.xe[0][ri + i] = g.xe[0][ri + i - 1] + g.dx[0][ri + i - 1];
    }

    //    g.xe[0][vector_type::iN()] =
    //        g.xe[0][vector_type::iN() - 1] + g.dx[0][vector_type::iN() - 1];
    /*

        for (size_t i = 0; i < NZ_GHOST; ++i) {
    //      auto idL = BCi_left_reflect(i);

          g.dx[0][idL.first] = g.dx[0][idL.second];

          for (size_t f = 0; f < N_FLUIDV; ++f) {
            vars[f][idL.first] = vars[f][idL.second];
          }
          flt.f[idL.first] = flt.f[idL.second];

          g.xe[0][idL.first] = g.xe[0][idL.first + 1] - g.dx[0][idL.first];
          g.xc[0][idL.first] = g.xe[0][idL.first] - 0.5 * g.dx[0][idL.first];

          auto idR = BCi_right_reflect(nx, i);
          g.dx[0][idR.first] = g.dx[0][idR.second];

          for (size_t f = 0; f < N_FLUIDV; ++f) {
            vars[f][idR.first] = vars[f][idR.second];
          }
          flt.f[idR.first] = flt.f[idR.second];

          g.xe[0][idR.first] = g.xe[0][idR.first - 1] + g.dx[0][idR.first - 1];
          g.xc[0][idR.first] = g.xe[0][idR.first] + 0.5 * g.dx[0][idR.first -
    1];

          // hack, make real later
          vars[XVL][idL.first] = -vars[XVL][idL.second];
          vars[XVL][idR.first] = -vars[XVL][idR.second];
        }*/
  }
};

}  // namespace ppm

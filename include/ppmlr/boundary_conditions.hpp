#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <utility>

//#include "gsl/gsl"

//#include "hydro_physics.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"
#include "ppm_vector.hpp"

namespace ppm {

template <typename Grid, index_type N>
struct boundary_conditions_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using grid_type = Grid;

  boundary_conditions_t() = default;
  ~boundary_conditions_t() = default;

  void reflect_bc(vector_type& v, double symm = 1.0) {
    for (index_type i = bp::i0; i < bp::j0; ++i) {
      v[bp::j0 - i - 1] = symm * v[i + bp::j0];
    }

    for (index_type i = bp::i0; i < bp::j0; ++i) {
      v[bp::jM + i] = symm * v[bp::jM - i - 1];
    }
  }

  inline void apply(grid_type& g, FluidVars<vector_type>& vars) {
    reflect_bc(g.dx[0]);
    reflect_bc(vars[RHO]);
    reflect_bc(vars[PRS]);
    reflect_bc(vars[XVL], -1.0);
    reflect_bc(vars[YVL]);
    reflect_bc(vars[ZVL]);
    reflect_bc(vars[ENT]);
    reflect_bc(vars[FLT]);

    for (index_type i = bp::i0; i < bp::j0; ++i) {
      g.xe[0][bp::j0 - i - 1] = g.xe[0][bp::j0 - i] - g.dx[0][bp::j0 - i - 1];

      g.xe[0][bp::jM + i] = g.xe[0][bp::jM + i - 1] + g.dx[0][bp::jM + i - 1];
    }
  }
};

}  // namespace ppm

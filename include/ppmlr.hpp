#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "boundary_conditions.hpp"
#include "hydro_physics.hpp"
#include "parabola.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"
#include "ppm_util.hpp"
#include "ppm_vector.hpp"
#include "riemann.hpp"
#include "state.hpp"

namespace ppm
{

template <typename PPMVector, typename GeoType, size_t NX>
struct ppmlr_t {
  using vector_type = PPMVector;
  using index_type = typename vector_type::index_type;

  using grid_type = ppm_grid_t<vector_type, GeoType>;
  using bc_type = boundary_conditions_t<vector_type, grid_type>;
  using parabola_type = parabola_t<vector_type>;
  using state_domain_type = state_domain_t<vector_type>;
  using state_type = state_t<parabola_type, state_domain_type>;

  using parabola_set = parabola_set_t<parabola_type>;
  using state_set = state_set_t<state_type, N_STATEV>;
  using riemann = riemann_t<vector_type, state_type>;

  const hydro_physics_t& hydro;  // reference to 'globals'

  // coordinate vars
  grid_type grid;
  // fluid vars
  FluidVars<vector_type> fluid_v;

  // boundary do-er
  bc_type bc;

  // flatten stuff
  //  flat_type flat;
  // parabolas
  parabola_set parabolas;
  // states
  state_set states;
  // riemann-er
  riemann riemann_solve;

  vector_type flat, z;
  // kinetic E store, cspeed
  vector_type ekin, cdx;
  // evolve stores
  vector_type dm, dm0, dtdm, pmid, umid, upmid, uold;
  // remap stores
  std::array<vector_type, N_FLUIDV> fluxes;
  // max svel in sweep
  double max_cdx;

  ppmlr_t(const hydro_physics_t& _hydro)
      : hydro(_hydro), max_cdx(0.0), states(PRIMATIVE_VARS) {}

  ~ppmlr_t() = default;

  template <typename V>
  inline void put_vars(const V& xe,
                       const V& xc,
                       const V& dx,
                       const V& rho,
                       const V& prs,
                       const V& xvel,
                       const V& yvel,
                       const V& zvel,
                       const V& f) {
    grid.xe[0].assign(xe);
    grid.xc[0].assign(xc);
    grid.dx[0].assign(dx);

    fluid_v[RHO].assign(rho);
    fluid_v[PRS].assign(prs);
    fluid_v[XVL].assign(xvel);
    fluid_v[YVL].assign(yvel);
    fluid_v[ZVL].assign(zvel);
    fluid_v[FLT].assign(f);

    build_ekin();
    for (index_type i = vector_type::j0; i < vector_type::jM; ++i)
      fluid_v[ENT][i] =
          fluid_v[PRS][i] / (fluid_v[RHO][i] * hydro.gammam1) + ekin[i];
  }

  template <typename V>
  inline void get_vars(V& rho, V& prs, V& xvel, V& yvel, V& zvel, V& f) {
    for (index_type i = vector_type::j0; i < vector_type::jM; ++i) {
      index_type n = i - vector_type::j0;
      rho[n] = fluid_v[RHO][i];
      prs[n] = fluid_v[PRS][i];
      xvel[n] = fluid_v[XVL][i];
      yvel[n] = fluid_v[YVL][i];
      zvel[n] = fluid_v[ZVL][i];
      f[n] = fluid_v[FLT][i];
    }
  }

  inline void build_ekin() {
    for (index_type i = vector_type::i0; i < vector_type::iN; ++i) {
      ekin[i] = fluid_v[XVL][i] * fluid_v[XVL][i] +
                fluid_v[YVL][i] * fluid_v[YVL][i] +
                fluid_v[ZVL][i] * fluid_v[ZVL][i];
      ekin[i] *= 0.5;
    }
  }

  inline void build_cdx() {
    max_cdx = 0.0;
    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      cdx[i] = sqrt(hydro.gamma * fluid_v[PRS][i] / fluid_v[RHO][i]) /
               (grid.dx[0][i] * grid.radius);
      max_cdx = std::max(max_cdx, cdx[i]);
    }
  }

  inline void flatten() {
    double delp1, delp2, zsten;
    double shock, old_flat;

    for (index_type i = vector_type::j0 - 4; i < vector_type::jM + 4; ++i) {
      delp1 = fluid_v[PRS][i + 1] - fluid_v[PRS][i - 1];
      delp2 = fluid_v[PRS][i + 2] - fluid_v[PRS][i - 2];
      if (std::abs(delp2) < small)
        delp2 = small;
      shock =
          std::abs(delp1) / std::min(fluid_v[PRS][i + 1], fluid_v[PRS][i - 1]) -
          epsilon;
      shock = std::max(0.0, shock);
      if (shock > 0.0)
        shock = 1.0;
      if (fluid_v[XVL][i - 1] < fluid_v[XVL][i + 1])
        shock = 0.0;

      z[i] = shock * std::max(0.0, (delp1 / delp2 - omega1) * omega2);
    }

    z[1] = z[2];
    z[vector_type::iN - 2] = z[vector_type::iN - 3];

    flat.fill(0.0);
    //  for (const auto& [i, im2, im1, ip1, ip2] : r1) {
    for (index_type i = vector_type::j0 - 4; i < vector_type::jM + 4; ++i) {
      zsten = std::max(z[i - 1], std::max(z[i], z[i + 1]));
      flat[i] = std::max(0.0, std::min(0.5, zsten));
    }

    for (index_type i = vector_type::j0 - 3; i < vector_type::jM + 3; ++i) {
      old_flat = fluid_v[FLT][i] - static_cast<int>(fluid_v[FLT][i]);
      if (flat[i] > 0.0) {
        flat[i] = std::max(flat[i], old_flat);
        fluid_v[FLT][i] = std::max(flat[i] - 1.0, 0.0);
      } else {
        fluid_v[FLT][i] = std::max(fluid_v[FLT][i] - 1.0, 0.0);
        flat[i] = old_flat;
      }
    }
  }

  inline void evolve() {
    grid.template build_volume<0>();

    for (index_type i = vector_type::i0 + 3; i < vector_type::iN - 2; ++i) {
      dm[i] = fluid_v[RHO][i] * grid.dvol[0][i];
      dtdm[i] = hydro.dt / dm[i];
      grid.xe[1][i] = grid.xe[0][i] + hydro.dt * umid[i] / grid.radius;
      upmid[i] = umid[i] * pmid[i];
    }

    for (index_type i = 0; i < 2; ++i) {
      grid.xe[1][vector_type::i0 + 2 - i] = grid.xe[0][vector_type::i0 + 2 - i];
      grid.xe[1][vector_type::iN - 2 + i] = grid.xe[0][vector_type::iN - 2 + i];
    }

    //    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 1; ++i)
    //    {
    for (index_type i = vector_type::i0; i < vector_type::iN - 1; ++i) {
      grid.dx[1][i] = grid.xe[1][i + 1] - grid.xe[1][i];
      grid.xc[1][i] = grid.xe[1][i] + 0.5 * grid.dx[1][i];
    }

    grid.template build_volume<1>();
    grid.build_amid();

    // dump_vec_and_keep_calm(vector_type::i0, vector_type::iN, grid.xe[0],
    //                       grid.xe[1]);

    for (index_type i = vector_type::i0 + 3; i < vector_type::iN - 3; ++i) {
      fluid_v[RHO][i] *= (grid.dvol[0][i] / grid.dvol[1][i]);

      uold[i] = fluid_v[XVL][i];

      fluid_v[XVL][i] -= dtdm[i] * (pmid[i + 1] - pmid[i]) * 0.5 *
                         (grid.amid[i + 1] + grid.amid[i]);

      fluid_v[ENT][i] -=
          dtdm[i] * (grid.amid[i + 1] * upmid[i + 1] - grid.amid[i] * upmid[i]);

      build_ekin();
      fluid_v[EIN][i] = fluid_v[ENT][i] - ekin[i];

      fluid_v[EIN][i] =
          std::max(fluid_v[EIN][i], smallp / (hydro.gammam1 * fluid_v[RHO][i]));
    }
  }

  inline void remap() {
    double fn1, fn2;

    parabolas.prepare(grid.dx[1]);

    parabolas.build(fluid_v, flat, CONSERVATIVE_VARS);

    /*    for (auto i = vector_type::j0; i < vector_type::jM + 1; ++i) {
          std::cout << i << " " << i + 1 << " " << fluid_v[RHO][i] << " "
                    << parabolas[RHO].da[i] << " " << parabolas[RHO].a6[i] << "
       "
                    << parabolas[RHO].al[i] << " " << flat[i];
          std::cout << std::endl;
        }
        std::exit(1);
    */
    grid.build_overlap();

    for (index_type f = 0; f < N_FLUIDV; ++f)
      fluxes[f].fill(0.0);

    for (index_type i = vector_type::j0; i < vector_type::jM + 1; ++i) {
      if (grid.xndiff[i] >= 0.0) {
        fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i - 1];
        fn2 = 1.0 - fourthird * fn1;

        fluxes[RHO][i] =
            parabolas[RHO].rc_left(i - 1, fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < N_FLUIDV; ++f) {
          fluxes[f][i] = parabolas[f].rc_left(i - 1, fn1, fn2) * fluxes[RHO][i];
        }
      } else {
        // tripple check
        fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i];
        fn2 = 1.0 + fourthird * fn1;

        fluxes[RHO][i] = parabolas[RHO].rc_right(i, -fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < N_FLUIDV; ++f) {
          fluxes[f][i] = parabolas[f].rc_right(i, -fn1, fn2) * fluxes[RHO][i];
        }
      }
    }

    /*    for (auto i = vector_type::j0; i < vector_type::jM + 1; ++i) {
          fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i - 1];
          fn2 = 1.0 - fourthird * fn1;

          std::cout << i << " " << i + 1 << " " << parabolas[RHO].al[i - 1] << "
       "
                    << parabolas[RHO].da[i - 1] << " " << parabolas[RHO].a6[i -
       1]
                    << " " << flat[i - 1];

          std::cout << std::endl;
        }
        std::exit(1);
    */
    //    dump_vec_and_keep_calm(vector_type::j0, vector_type::jM + 1,
    //    fluxes[RHO]);

    for (index_type i = vector_type::j0 - 1; i < vector_type::jM + 1; ++i) {
      dm[i] = fluid_v[RHO][i] * grid.dvol[1][i];
      dm0[i] = (dm[i] + fluxes[RHO][i] - fluxes[RHO][i + 1]);

      fluid_v[RHO][i] = dm0[i] / grid.dvol[0][i];
      fluid_v[RHO][i] = std::max(smallr, fluid_v[RHO][i]);

      dm0[i] = 1. / (fluid_v[RHO][i] * grid.dvol[0][i]);

      for (index_type f = 1; f < N_FLUIDV; ++f) {
        fluid_v[f][i] =
            (fluid_v[f][i] * dm[i] + fluxes[f][i] - fluxes[f][i + 1]) * dm0[i];
      }
    }

    // supersonic flow
    build_ekin();
    for (index_type i = vector_type::j0; i < vector_type::jM; ++i) {
      if (ekin[i] / fluid_v[EIN][i] < 100.0)
        fluid_v[EIN][i] = fluid_v[ENT][i] - ekin[i];
      fluid_v[PRS][i] = hydro.gammam1 * fluid_v[RHO][i] * fluid_v[EIN][i];
      fluid_v[PRS][i] = std::max(smallp, fluid_v[PRS][i]);
    }

  }  // remap

  inline void wiggle() {
    if (xwig * flat.sum() != 0.0) {
      //      std::cout << "DO WHIG\n";
      grid.store_coords();

      //      for (size_t i = NZ_GHOST + 1; i < nx + NZ_GHOST; ++i) {
      for (index_type i = vector_type::j0 + 1; i < vector_type::jM; ++i) {
        if (std::max(flat[i - 1], flat[i]) > 0.0) {
          grid.xe[0][i] += xwig * grid.dx[0][i];
          grid.dx[0][i - 1] = grid.xe[0][i] - grid.xe[0][i - 1];
        }
      }

      grid.template build_volume<0>();
      remap();

      grid.retieve_coords();
    }
  }

  inline void operator()() {
    bc.apply(grid, fluid_v);

    build_cdx();

    flatten();

    parabolas.prepare(grid.dx[0]);

    parabolas.build(fluid_v, flat, PRIMATIVE_VARS);

    states.prepare(hydro, cdx);

    states.build(parabolas.arr());

    riemann_solve(hydro, states, pmid, umid);

    evolve();

    wiggle();

    remap();
  }
};

}  // namespace ppm

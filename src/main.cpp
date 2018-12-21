#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <valarray>

#include "gsl/gsl"

#include "hydro_physics.hpp"
#include "ppm_grid.hpp"
#include "ppm_vector.hpp"
#include "ppmlr.hpp"

#define _GNU_SOURCE
#include <fenv.h>

constexpr size_t NCYC = 50;

constexpr double pr = 0.1;
constexpr double dr = 0.125;
constexpr double pl = 1.0;
constexpr double dl = 1.0;

constexpr double my_gamma = 1.4;

using VecType = std::valarray<double>;
using SweepVec = ppm::ppmlr_vector_t<VecType, ppm::NZ_GHOST, ppm::GRID_NX>;

using SWEEP_XGEO = ppm::geo::planar_geo_t;

template <typename V>
auto gridder(size_t n, double xmin, double xmax) {
  V xe(n + 1), xc(n), dx(n);

  dx = (xmax - xmin) / static_cast<double>(n);

  for (size_t i = 0; i < n; ++i) {
    xe[i] = xmin + static_cast<double>(i) * dx[i];
    xc[i] = xe[i] + 0.5 * dx[i];
  }
  xe[n] = xe[n-1] + dx[n-1];

  return std::make_tuple(xe, xc, dx);
}

template <typename V>
auto initial_conditions(size_t n, std::valarray<double>& xc) {
  V rho(n), prs(n), ene(n), fla(n);
  std::array<V, 3> vels;

  fla = 0.0;

  for (size_t i = 0; i < n; ++i) {
    if (xc[i] <= 0.5) {
      rho[i] = dl;
      prs[i] = pl;
    } else {
      rho[i] = dr;
      prs[i] = pr;
    }
  }

  for (auto& v : vels) {
    v.resize(n);
    v = 0.0;
  }
  /*  for (size_t i = 0; i < 2; ++i) {
      vels[0][n - i - 1] = -0.05;
    }
  */

  ene = prs / (rho * (my_gamma - 1)) +
        0.5 * (vels[0] * vels[0] + vels[1] * vels[1] + vels[2] * vels[2]);

  return std::make_tuple(rho, prs, vels, ene, fla);
}

template <typename... Vs>
void printout(size_t N, Vs&&... vs) {
  std::string s;
  auto commd = [](double a) { return std::to_string(a) + ", "; };
  for (size_t i = 0; i < N; ++i) {
    s = "";
    s = (commd(vs[i]) + ...);
    std ::cout << i << ": " << s << std::endl;
  }
}

template <typename... Vs>
void dump(size_t N, double t, Vs&&... vs) {
  std::ofstream of("run00.dat", std::ofstream::out | std::ofstream::app);

  std::string s;
  auto commd = [](double a) { return std::to_string(a) + ", "; };
  of << t << std::endl;
  for (size_t i = 0; i < N; ++i) {
    s = "";
    s = (commd(vs[i]) + ...);
    of << i << " " << s << std::endl;
  }
  of << "\n\n";
  of.close();
}

/*template <typename V>
void courant(ppm::hydro_physics_t& h,
             const V& vx,
             const V& xdx) {
  double ridt = 0.0;

  //  for (size_t i = 0; i < nx; ++i) {
  //    ridt = std::max(std::abs(vx[i]) / xdx[i], ridt);
  //  }

  ridt = (std::abs(vx) / xdx).max();
  ridt = std::max(h.svel, ridt);

  h.set_dt(ppm::courant / ridt);

}
*/

struct sweep_x {
  ppm::ppmlr_t<SweepVec, SWEEP_XGEO, ppm::GRID_NX> ppmlr_drv;

  sweep_x(const ppm::hydro_physics_t& h) : ppmlr_drv(h) {}
};
int main()
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  ppm::hydro_physics_t hydro(my_gamma);
  sweep_x swpx(hydro);

  auto [xe, xc, dx] = gridder<VecType>(ppm::GRID_NX, 0.0, 1.0);
  auto [r, p, u, e, f] = initial_conditions<VecType>(ppm::GRID_NX, xc);

  /*  auto svel = (std::sqrt(hydro.gamma * p / r) / dx).max();*/
  hydro.first_dt(p, r, u[0], dx);
  double time = 0.0;
  auto start_hydro = std::chrono::high_resolution_clock::now();

  size_t ncycle = 0;
  for (size_t it = 0; it < NCYC; ++it) {
    time += hydro.dt;
    ++ncycle;

    e = p / (r * hydro.gammam1) +
        0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

    dump(ppm::GRID_NX, time, xe, r, p, u[0], e, f);

    swpx.ppmlr_drv.put_vars(xe, xc, dx, r, p, u[0], u[1], u[2], f);

    //    sweepx();
    swpx.ppmlr_drv();

    swpx.ppmlr_drv.get_vars(r, p, u[0], u[1], u[2], f);

    hydro.dt_courant(swpx.ppmlr_drv.max_cdx, u[0], dx);
  }
    auto end_hydro = std::chrono::high_resolution_clock::now();

    std::cout << time << std::endl;
    auto hydro_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                              end_hydro - start_hydro)
                              .count();

    std::cout << "ttf = " << hydro_duration << " ms\n";

  /* for (size_t i = 0; i < NX; ++i) */
  /*   std::cout << i << " : " << swp.fluid_v[ppm::RHO][i + ppm::NZ_GHOST] */
  /*             << std::endl; */

  return 0;
}


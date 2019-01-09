#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <valarray>

#include "hydro_physics.hpp"
#include "ppmlr/ppm_grid.hpp"
#include "ppmlr/ppm_vector.hpp"
#include "ppmlr/ppmlr.hpp"
#include "run_def.hpp"
#include "vector_wrapper.hpp"

//#define _GNU_SOURCE
//#include <fenv.h>

constexpr size_t NCYC = 1500;

constexpr double pr = 0.1;
constexpr double dr = 0.125;
constexpr double pl = 1.0;
constexpr double dl = 1.0;

constexpr double my_gamma = 1.4;

template <index_type N>
using VecType = vector_wrapper_t<N>;

template <typename V>
constexpr inline auto etot(size_t n,
                           const V& p,
                           const V& r,
                           const std::array<V, 3>& vs,
                           V& e) {
  double gm1 = (my_gamma - 1.0);

  for (size_t i = 0; i < n; ++i) {
    e[i] =
        p[i] / (r[i] * gm1) +
        0.5 * (vs[0][i] * vs[0][i] + vs[1][i] * vs[1][i] + vs[2][i] * vs[2][i]);
  }
}

auto gridder(double xmin, double xmax) {
  VecType<hydro::NX> xe;
  VecType<hydro::NX> xc, dx;

  auto ddx = (xmax - xmin) / static_cast<double>(hydro::NX);

  for (size_t i = 0; i < hydro::NX; ++i) {
    dx[i] = ddx;
    xe[i] = xmin + static_cast<double>(i) * dx[i];
    xc[i] = xe[i] + 0.5 * dx[i];
  }
  //  xe[hydro::NX] = xe[hydro::NX - 1] + dx[hydro::NX - 1];

  return std::make_tuple(xe, xc, dx);
}

template <template <index_type> class V, index_type N>
auto initial_conditions(const V<N>& xc) {
  V<N> rho, prs, ene, fla;
  std::array<V<N>, 3> vels;

  for (size_t i = 0; i < N; ++i) {
    fla[i] = 0.0;
    if (xc[i] <= 0.5) {
      rho[i] = dl;
      prs[i] = pl;
    } else {
      rho[i] = dr;
      prs[i] = pr;
    }
  }

  for (auto& v : vels) {
    for (size_t i = 0; i < N; ++i) {
      v[i] = 0.0;
    }
  }

  etot(N, prs, rho, vels, ene);

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
  ppm::ppmlr_t<hydro::XGEO, hydro::NX> ppmlr_drv;

  sweep_x(const hydro::hydro_physics_t& h) : ppmlr_drv(h) {}
};
int main()
{
  //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  hydro::hydro_physics_t hydro(my_gamma);
  sweep_x swpx(hydro);

  auto [xe, xc, dx] = gridder(0.0, 1.0);
  auto [r, p, u, e, f] = initial_conditions<VecType, hydro::NX>(xc);

  std::chrono::duration<double> rt_s(0);

  hydro.first_dt(p, r, u[0], dx);
  double time = 0.0;

  size_t ncycle = 0;
  for (size_t it = 0; it < NCYC; ++it) {
    dump(hydro::NX, time, xe, r, p, u[0], e, f);

    auto start_hydro = std::chrono::high_resolution_clock::now();

    time += hydro.dt;
    ++ncycle;

    etot(hydro::NX, p, r, u, e);

    auto end_hydro = std::chrono::high_resolution_clock::now();
    rt_s += (end_hydro - start_hydro);

    swpx.ppmlr_drv.put_vars(xe, xc, dx, r, p, u[0], u[1], u[2], f);

    start_hydro = std::chrono::high_resolution_clock::now();
    //    sweepx();
    swpx.ppmlr_drv();

    end_hydro = std::chrono::high_resolution_clock::now();
    rt_s += (end_hydro - start_hydro);

    swpx.ppmlr_drv.get_vars(r, p, u[0], u[1], u[2], f);

    start_hydro = std::chrono::high_resolution_clock::now();

    hydro.dt_courant(swpx.ppmlr_drv.max_cdx, u[0], dx);

    end_hydro = std::chrono::high_resolution_clock::now();
    rt_s += (end_hydro - start_hydro);
  }

    std::cout << time << std::endl;

    //    std::cout << "ttf = " << hydro_duration.count() << " us\n";
    std::cout << "ttf = " << rt_s.count() * 1.0E-6 << " us\n";

    /* for (size_t i = 0; i < NX; ++i) */
    /*   std::cout << i << " : " << swp.fluid_v[ppm::RHO][i + ppm::NZ_GHOST] */
    /*             << std::endl; */

    return 0;
}


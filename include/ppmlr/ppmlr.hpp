#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>

#include "boundary_conditions.hpp"
#include "compute_kernels.hpp"
#include "fluid_data.hpp"
#include "hydro_physics.hpp"
#include "parabola.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"
#include "ppm_vector.hpp"
#include "riemann_solve.hpp"
#include "state.hpp"

namespace ppm
{

template <typename GeoType, index_type N>
struct ppmlr_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using grid_type = ppm_grid_t<GeoType, N>;
  using bc_type = boundary_conditions_t<grid_type, N>;
  using parabola_type = parabola_t<N>;
  using state_domain_type = state_domain_t<N>;
  using state_type = state_t<parabola_type, N>;

  using parabola_set = parabola_set_t<parabola_type, N>;
  using state_set = state_set_t<state_type, N>;
  using riemann = riemann_t<state_type, N>;

  const hydro::hydro_physics_t& hydro;  // reference to 'globals'

  // coordinate vars
  grid_type grid;

  // fluid vars - individual vectors for each variable
  vector_type rho, prs, xvl, yvl, zvl, ent, ein, flt;

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

  // remap stores - individual flux vectors for each variable
  vector_type flux_rho, flux_prs, flux_xvl, flux_yvl, flux_zvl;
  vector_type flux_ent, flux_ein, flux_flt;

  // max svel in sweep
  double max_cdx;

  ppmlr_t(const hydro::hydro_physics_t& _hydro)
      : hydro(_hydro), max_cdx(zero),
        rho(base::make_vector()), prs(base::make_vector()),
        xvl(base::make_vector()), yvl(base::make_vector()),
        zvl(base::make_vector()), ent(base::make_vector()),
        ein(base::make_vector()), flt(base::make_vector()),
        flat(base::make_vector()), z(base::make_vector()),
        ekin(base::make_vector()), cdx(base::make_vector()),
        dm(base::make_vector()), dm0(base::make_vector()),
        dtdm(base::make_vector()), pmid(base::make_vector()),
        umid(base::make_vector()), upmid(base::make_vector()),
        uold(base::make_vector()),
        flux_rho(base::make_vector()), flux_prs(base::make_vector()),
        flux_xvl(base::make_vector()), flux_yvl(base::make_vector()),
        flux_zvl(base::make_vector()), flux_ent(base::make_vector()),
        flux_ein(base::make_vector()), flux_flt(base::make_vector()) {}

  ~ppmlr_t() = default;

  template <typename V>
  inline void put_vars(const V& xe,
                       const V& xc,
                       const V& dx,
                       const V& rho_in,
                       const V& prs_in,
                       const V& xvel,
                       const V& yvel,
                       const V& zvel,
                       const V& f) {
    // Explicit memory management: Host → Device transfer
    // Create unmanaged views wrapping std::vector data (no copy)
    const auto offset = bp::j0;
    const auto count = rho_in.size();

    // Wrap input std::vectors as unmanaged host views
    auto h_xe = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(xe.data(), count);
    auto h_xc = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(xc.data(), count);
    auto h_dx = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(dx.data(), count);
    auto h_rho = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(rho_in.data(), count);
    auto h_prs = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(prs_in.data(), count);
    auto h_xvel = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(xvel.data(), count);
    auto h_yvel = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(yvel.data(), count);
    auto h_zvel = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(zvel.data(), count);
    auto h_f = Kokkos::View<const double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(f.data(), count);

    // Create subviews of device arrays (active zone only, excluding ghost cells)
    auto d_xe = Kokkos::subview(grid.xe[0], Kokkos::make_pair(offset, offset + count));
    auto d_xc = Kokkos::subview(grid.xc[0], Kokkos::make_pair(offset, offset + count));
    auto d_dx = Kokkos::subview(grid.dx[0], Kokkos::make_pair(offset, offset + count));
    auto d_rho = Kokkos::subview(rho, Kokkos::make_pair(offset, offset + count));
    auto d_prs = Kokkos::subview(prs, Kokkos::make_pair(offset, offset + count));
    auto d_xvl = Kokkos::subview(xvl, Kokkos::make_pair(offset, offset + count));
    auto d_yvl = Kokkos::subview(yvl, Kokkos::make_pair(offset, offset + count));
    auto d_zvl = Kokkos::subview(zvl, Kokkos::make_pair(offset, offset + count));
    auto d_flt = Kokkos::subview(flt, Kokkos::make_pair(offset, offset + count));

    // Explicit deep copy: Host → Device
    Kokkos::deep_copy(d_xe, h_xe);
    Kokkos::deep_copy(d_xc, h_xc);
    Kokkos::deep_copy(d_dx, h_dx);
    Kokkos::deep_copy(d_rho, h_rho);
    Kokkos::deep_copy(d_prs, h_prs);
    Kokkos::deep_copy(d_xvl, h_xvel);
    Kokkos::deep_copy(d_yvl, h_yvel);
    Kokkos::deep_copy(d_zvl, h_zvel);
    Kokkos::deep_copy(d_flt, h_f);

    build_ekin();

    // DEVICE OFFLOAD READY: Entropy computation
    kernels::compute_entropy(
        prs.data(), rho.data(), ekin.data(), ent.data(),
        hydro.gammam1, bp::j0, bp::jM);
  }

  template <typename V>
  inline void get_vars(V& rho_out, V& prs_out, V& xvel, V& yvel, V& zvel, V& f) {
    // Explicit memory management: Device → Host transfer
    const auto count = rho_out.size();

    // Create subviews of device arrays (active zone only)
    auto d_rho = Kokkos::subview(rho, Kokkos::make_pair(bp::j0, bp::jM));
    auto d_prs = Kokkos::subview(prs, Kokkos::make_pair(bp::j0, bp::jM));
    auto d_xvl = Kokkos::subview(xvl, Kokkos::make_pair(bp::j0, bp::jM));
    auto d_yvl = Kokkos::subview(yvl, Kokkos::make_pair(bp::j0, bp::jM));
    auto d_zvl = Kokkos::subview(zvl, Kokkos::make_pair(bp::j0, bp::jM));
    auto d_flt = Kokkos::subview(flt, Kokkos::make_pair(bp::j0, bp::jM));

    // Wrap output std::vectors as unmanaged host views
    auto h_rho = Kokkos::View<double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(rho_out.data(), count);
    auto h_prs = Kokkos::View<double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(prs_out.data(), count);
    auto h_xvel = Kokkos::View<double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(xvel.data(), count);
    auto h_yvel = Kokkos::View<double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(yvel.data(), count);
    auto h_zvel = Kokkos::View<double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(zvel.data(), count);
    auto h_f = Kokkos::View<double*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>(f.data(), count);

    // Explicit deep copy: Device → Host
    Kokkos::deep_copy(h_rho, d_rho);
    Kokkos::deep_copy(h_prs, d_prs);
    Kokkos::deep_copy(h_xvel, d_xvl);
    Kokkos::deep_copy(h_yvel, d_yvl);
    Kokkos::deep_copy(h_zvel, d_zvl);
    Kokkos::deep_copy(h_f, d_flt);
  }

  inline void build_ekin() {
    // DEVICE OFFLOAD READY: Use compute kernel for kinetic energy
    // This is embarrassingly parallel and perfect for GPU offload
    kernels::compute_kinetic_energy(
        xvl.data(), yvl.data(), zvl.data(), ekin.data(),
        bp::i0, bp::iN);
  }

  inline void build_cdx() {
    // DEVICE OFFLOAD READY: Sound speed with parallel reduction
    max_cdx = kernels::compute_sound_speed(
        prs.data(), rho.data(), grid.dx[0].data(), cdx.data(),
        hydro.gamma, grid.radius,
        bp::i0 + 2, bp::iN - 2);
  }

  inline void flatten() {
    double delp1, delp2, zsten;
    double shock, old_flat;

    for (index_type i = bp::j0 - 4; i < bp::jM + 4; ++i) {
      delp1 = prs[i + 1] - prs[i - 1];
      delp2 = prs[i + 2] - prs[i - 2];
      if (std::abs(delp2) < small)
        delp2 = small;
      shock =
          std::abs(delp1) / std::min(prs[i + 1], prs[i - 1]) -
          epsilon;
      shock = std::max(0.0, shock);
      if (shock > 0.0)
        shock = 1.0;
      if (xvl[i - 1] < xvl[i + 1])
        shock = 0.0;

      z[i] = shock * std::max(0.0, (delp1 / delp2 - omega1) * omega2);
    }

    z[1] = z[2];
    z[bp::iN - 2] = z[bp::iN - 3];

    // Zero-fill flat vector using Kokkos
    Kokkos::deep_copy(flat, 0.0);
    //  for (const auto& [i, im2, im1, ip1, ip2] : r1) {
    for (index_type i = bp::j0 - 4; i < bp::jM + 4; ++i) {
      zsten = std::max(z[i - 1], std::max(z[i], z[i + 1]));
      flat[i] = std::max(0.0, std::min(0.5, zsten));
    }

    for (index_type i = bp::j0 - 3; i < bp::jM + 3; ++i) {
      old_flat = flt[i] - static_cast<int>(flt[i]);
      if (flat[i] > 0.0) {
        flat[i] = std::max(flat[i], old_flat);
        flt[i] = std::max(flat[i] - 1.0, 0.0);
      } else {
        flt[i] = std::max(flt[i] - 1.0, 0.0);
        flat[i] = old_flat;
      }
    }
  }

  inline void evolve() {
    grid.template build_volume<0>();

    // DEVICE OFFLOAD READY: Grid mass and timestep calculation
    kernels::compute_grid_mass(
        dm.data(), dtdm.data(), rho.data(), grid.dvol[0].data(),
        hydro.dt, bp::i0 + 3, bp::iN - 2);

    // Grid update (sequential, stays on host)
    for (index_type i = bp::i0 + 3; i < bp::iN - 2; ++i) {
      grid.xe[1][i] = grid.xe[0][i] + hydro.dt * umid[i] / grid.radius;
      upmid[i] = umid[i] * pmid[i];
    }

    for (index_type i = 0; i < 2; ++i) {
      grid.xe[1][bp::i0 + 2 - i] = grid.xe[0][bp::i0 + 2 - i];
      grid.xe[1][bp::iN - 2 + i] = grid.xe[0][bp::iN - 2 + i];
    }

    for (index_type i = bp::i0; i < bp::iN - 1; ++i) {
      grid.dx[1][i] = grid.xe[1][i + 1] - grid.xe[1][i];
      grid.xc[1][i] = grid.xe[1][i] + 0.5 * grid.dx[1][i];
    }

    grid.template build_volume<1>();
    grid.build_amid();

    // DEVICE OFFLOAD READY: Density evolution (Lagrangian step)
    kernels::evolve_density(
        rho.data(), grid.dvol[0].data(), grid.dvol[1].data(),
        bp::i0 + 3, bp::iN - 3);

    // DEVICE OFFLOAD READY: Velocity evolution (momentum update)
    kernels::evolve_velocity(
        xvl.data(), uold.data(), pmid.data(), grid.amid.data(), dtdm.data(),
        bp::i0 + 3, bp::iN - 3);

    // DEVICE OFFLOAD READY: Energy evolution
    kernels::evolve_energy(
        ent.data(), grid.amid.data(), upmid.data(), dtdm.data(),
        bp::i0 + 3, bp::iN - 3);

    build_ekin();

    // DEVICE OFFLOAD READY: Internal energy computation
    kernels::compute_internal_energy(
        ein.data(), ent.data(), ekin.data(), rho.data(),
        hydro.gammam1, bp::i0 + 3, bp::iN - 3);
  }

  inline void remap() {
    double fn1, fn2;

    parabolas.prepare(grid.dx[1]);

    parabolas.build_conservative(rho, xvl, yvl, zvl, ein, ent, flat);

    grid.build_overlap();

    // Clear all flux vectors using Kokkos deep_copy
    Kokkos::deep_copy(flux_rho, 0.0);
    Kokkos::deep_copy(flux_prs, 0.0);
    Kokkos::deep_copy(flux_xvl, 0.0);
    Kokkos::deep_copy(flux_yvl, 0.0);
    Kokkos::deep_copy(flux_zvl, 0.0);
    Kokkos::deep_copy(flux_ent, 0.0);
    Kokkos::deep_copy(flux_ein, 0.0);
    Kokkos::deep_copy(flux_flt, 0.0);

    for (index_type i = bp::j0; i < bp::jM + 1; ++i) {
      if (grid.xndiff[i] >= 0.0) {
        fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i - 1];
        fn2 = 1.0 - fourthird * fn1;

        flux_rho[i] = parabolas.para_rho.rc_left(i - 1, fn1, fn2) * grid.xnolap[i];

        // Compute fluxes for other variables
        flux_prs[i] = parabolas.para_prs.rc_left(i - 1, fn1, fn2) * flux_rho[i];
        flux_xvl[i] = parabolas.para_xvl.rc_left(i - 1, fn1, fn2) * flux_rho[i];
        flux_yvl[i] = parabolas.para_yvl.rc_left(i - 1, fn1, fn2) * flux_rho[i];
        flux_zvl[i] = parabolas.para_zvl.rc_left(i - 1, fn1, fn2) * flux_rho[i];
        flux_ent[i] = parabolas.para_ent.rc_left(i - 1, fn1, fn2) * flux_rho[i];
        flux_ein[i] = parabolas.para_ein.rc_left(i - 1, fn1, fn2) * flux_rho[i];
        flux_flt[i] = parabolas.para_flt.rc_left(i - 1, fn1, fn2) * flux_rho[i];
      } else {
        // tripple check
        fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i];
        fn2 = 1.0 + fourthird * fn1;

        flux_rho[i] = parabolas.para_rho.rc_right(i, -fn1, fn2) * grid.xnolap[i];

        // Compute fluxes for other variables
        flux_prs[i] = parabolas.para_prs.rc_right(i, -fn1, fn2) * flux_rho[i];
        flux_xvl[i] = parabolas.para_xvl.rc_right(i, -fn1, fn2) * flux_rho[i];
        flux_yvl[i] = parabolas.para_yvl.rc_right(i, -fn1, fn2) * flux_rho[i];
        flux_zvl[i] = parabolas.para_zvl.rc_right(i, -fn1, fn2) * flux_rho[i];
        flux_ent[i] = parabolas.para_ent.rc_right(i, -fn1, fn2) * flux_rho[i];
        flux_ein[i] = parabolas.para_ein.rc_right(i, -fn1, fn2) * flux_rho[i];
        flux_flt[i] = parabolas.para_flt.rc_right(i, -fn1, fn2) * flux_rho[i];
      }
    }

    /*    for (auto i = bp::j0; i < bp::jM + 1; ++i) {
          fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i - 1];
          fn2 = 1.0 - fourthird * fn1;

          std::cout << i << " " << i + 1 << " " << parabolas.para_rho.al[i - 1] << "
       "
                    << parabolas.para_rho.da[i - 1] << " " << parabolas.para_rho.a6[i -
       1]
                    << " " << flat[i - 1];

          std::cout << std::endl;
        }
        std::exit(1);
    */

    // DEVICE OFFLOAD READY: Density remap
    kernels::remap_density(
        rho.data(), flux_rho.data(),
        grid.dvol[0].data(), grid.dvol[1].data(),
        dm.data(), dm0.data(),
        bp::j0 - 1, bp::jM + 1);

    // DEVICE OFFLOAD READY: Velocity and energy remap
    kernels::remap_scalar(xvl.data(), flux_xvl.data(), dm.data(), dm0.data(),
                          bp::j0 - 1, bp::jM + 1);
    kernels::remap_scalar(yvl.data(), flux_yvl.data(), dm.data(), dm0.data(),
                          bp::j0 - 1, bp::jM + 1);
    kernels::remap_scalar(zvl.data(), flux_zvl.data(), dm.data(), dm0.data(),
                          bp::j0 - 1, bp::jM + 1);
    kernels::remap_scalar(ein.data(), flux_ein.data(), dm.data(), dm0.data(),
                          bp::j0 - 1, bp::jM + 1);
    kernels::remap_scalar(ent.data(), flux_ent.data(), dm.data(), dm0.data(),
                          bp::j0 - 1, bp::jM + 1);

    // DEVICE OFFLOAD READY: Supersonic flow and pressure update
    build_ekin();
    kernels::compute_pressure(
        prs.data(), rho.data(), ein.data(), ent.data(), ekin.data(),
        hydro.gammam1, bp::j0, bp::jM);

  }  // remap

  inline void wiggle() {
    // Sum flat array on host (small operation, wiggle is kept on CPU)
    // Create host mirror and copy
    auto h_flat = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, flat);
    double flat_sum = 0.0;
    for (index_type i = 0; i < bp::iN; ++i) {
      flat_sum += h_flat(i);
    }

    if (xwig * flat_sum != 0.0) {
      //      std::cout << "DO WHIG\n";
      grid.store_coords();

      //      for (size_t i = NZ_GHOST + 1; i < nx + NZ_GHOST; ++i) {
      for (index_type i = bp::j0 + 1; i < bp::jM; ++i) {
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
    bc.apply(grid, rho, prs, xvl, yvl, zvl, ent, flt);

    build_cdx();

    flatten();

    parabolas.prepare(grid.dx[0]);

    parabolas.build_primitive(rho, prs, xvl, flat);

    states.prepare(hydro, cdx);

    states.build(parabolas);

    riemann_solve(hydro, states, pmid, umid);

    evolve();

    wiggle();

    remap();
  }
};

}  // namespace ppm

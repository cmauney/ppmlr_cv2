#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

#include <Kokkos_Core.hpp>

#include "hydro_physics.hpp"
#include "ppm_constant.hpp"
#include "ppm_vector.hpp"
#include "state.hpp"

namespace ppm
{

    constexpr double NewRapTol = 1.0E-5;
    constexpr size_t NewRapItr = 12;

    template <typename StateType, index_type N>
    struct riemann_t : ppmlr_base_t<N> {
      using base = ppmlr_base_t<N>;
      using bp = typename base::bp;
      using vector_type = typename base::vector_type;

      using state_set = state_set_t<StateType, N, N_STATEV>;

      double gamfac1, gamfac2;

      vector_type cl, cr;
      vector_type pli, pri, rli, rri;
      vector_type wl, wr, zl, zr;
      vector_type umidl, umidr, pmold;

      riemann_t() : cl(base::make_vector()), cr(base::make_vector()),
                    pli(base::make_vector()), pri(base::make_vector()),
                    rli(base::make_vector()), rri(base::make_vector()),
                    wl(base::make_vector()), wr(base::make_vector()),
                    zl(base::make_vector()), zr(base::make_vector()),
                    umidl(base::make_vector()), umidr(base::make_vector()),
                    pmold(base::make_vector()) {}

      inline void operator()(const hydro::hydro_physics_t& hydro,
                             const state_set& states,
                             vector_type& pmid,
                             vector_type& umid) {
        const double gamma = hydro.gamma;
        const double gamfac1 = gamma + one;
        const double gamfac2 = half * gamfac1 / gamma;

        // Phase 1: Compute left state sound speeds and reciprocals (GPU-accelerated)
        Kokkos::parallel_for("riemann_left_setup",
            Kokkos::RangePolicy<>(bp::i0 + 3, bp::iN - 2),
            KOKKOS_LAMBDA(const index_type i) {
              cl(i) = Kokkos::sqrt(gamma * states[PRS].l(i) * states[RHO].l(i));
              pli(i) = one / states[PRS].l(i);
              rli(i) = one / states[RHO].l(i);
            });

        // Phase 2: Compute right state sound speeds and reciprocals (GPU-accelerated)
        Kokkos::parallel_for("riemann_right_setup",
            Kokkos::RangePolicy<>(bp::i0 + 3, bp::iN - 2),
            KOKKOS_LAMBDA(const index_type i) {
              cr(i) = Kokkos::sqrt(gamma * states[PRS].r(i) * states[RHO].r(i));
              pri(i) = one / states[PRS].r(i);
              rri(i) = one / states[RHO].r(i);
            });

        // Phase 3: Initial pressure guess from acoustic approximation (GPU-accelerated)
        Kokkos::parallel_for("riemann_initial_guess",
            Kokkos::RangePolicy<>(bp::i0 + 3, bp::iN - 2),
            KOKKOS_LAMBDA(const index_type i) {
              double pm = states[PRS].r(i) - states[PRS].l(i) -
                          cr(i) * (states[XVL].r(i) - states[XVL].l(i));
              pm = states[PRS].l(i) + pm * cl(i) / (cl(i) + cr(i));
              pmid(i) = Kokkos::fmax(small_riemann, pm);
            });

        // Phase 4: Newton-Raphson iteration (GPU-accelerated)
        // Fixed iteration count for GPU efficiency (no early exit)
        Kokkos::parallel_for("riemann_newton_raphson",
            Kokkos::RangePolicy<>(bp::i0 + 3, bp::iN - 2),
            KOKKOS_LAMBDA(const index_type i) {
              // Newton-Raphson: 12 iterations (sufficient for convergence in practice)
              // Early exit removed for GPU warp efficiency
              for (index_type n = 0; n < NewRapItr; ++n) {
                const double pm = pmid(i);

                // Compute wave speeds
                double wl_val = one + gamfac2 * (pm - states[PRS].l(i)) * pli(i);
                double wr_val = one + gamfac2 * (pm - states[PRS].r(i)) * pri(i);

                wl_val = cl(i) * Kokkos::sqrt(wl_val);
                wr_val = cr(i) * Kokkos::sqrt(wr_val);

                // Compute impedances
                double zl_val = 4.0 * rli(i) * wl_val * wl_val;
                double zr_val = 4.0 * rri(i) * wr_val * wr_val;

                zl_val = -zl_val * wl_val / (zl_val - gamfac1 * (pm - states[PRS].l(i)));
                zr_val = zr_val * wr_val / (zr_val - gamfac1 * (pm - states[PRS].r(i)));

                // Store for final velocity computation
                wl(i) = wl_val;
                wr(i) = wr_val;

                // Compute velocities
                const double ul = states[XVL].l(i) - (pm - states[PRS].l(i)) / wl_val;
                const double ur = states[XVL].r(i) + (pm - states[PRS].r(i)) / wr_val;

                // Newton-Raphson update
                double pm_new = pm + (ur - ul) * (zl_val * zr_val) / (zr_val - zl_val);
                pmid(i) = Kokkos::fmax(small_riemann, pm_new);
              }
            });

        // Phase 5: Final velocity computation (GPU-accelerated)
        Kokkos::parallel_for("riemann_final_velocity",
            Kokkos::RangePolicy<>(bp::i0 + 3, bp::iN - 2),
            KOKKOS_LAMBDA(const index_type i) {
              const double ul = states[XVL].l(i) - (pmid(i) - states[PRS].l(i)) / wl(i);
              const double ur = states[XVL].r(i) + (pmid(i) - states[PRS].r(i)) / wr(i);
              umid(i) = half * (ul + ur);
            });
      }
    };

  }  // namespace ppm

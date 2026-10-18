#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

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
        double gamfac1 = hydro.gamma + one;
        double gamfac2 = half * gamfac1 / hydro.gamma;

        //        auto lo = NZ_GHOST - 3, hi = nx + 2 * NZ_GHOST - 2;
        for (index_type i = bp::i0 + 3; i < bp::iN - 2; ++i) {
          cl[i] = std::sqrt(hydro.gamma * states[PRS].l[i] * states[RHO].l[i]);
          pli[i] = one / states[PRS].l[i];
          rli[i] = one / states[RHO].l[i];
        }

        for (index_type i = bp::i0 + 3; i < bp::iN - 2; ++i) {
          cr[i] = std::sqrt(hydro.gamma * states[PRS].r[i] * states[RHO].r[i]);
          pri[i] = one / states[PRS].r[i];
          rri[i] = one / states[RHO].r[i];
        }

        for (index_type i = bp::i0 + 3; i < bp::iN - 2; ++i) {
          pmid[i] = states[PRS].r[i] - states[PRS].l[i] -
                    cr[i] * (states[XVL].r[i] - states[XVL].l[i]);
          pmid[i] = states[PRS].l[i] + pmid[i] * cl[i] / (cl[i] + cr[i]);
          pmid[i] = std::max(small_riemann, pmid[i]);
        }

        for (index_type i = bp::i0 + 3; i < bp::iN - 2; ++i) {
          for (index_type n = 0; n < NewRapItr; ++n) {
            pmold[i] = pmid[i];
            wl[i] = one + gamfac2 * (pmid[i] - states[PRS].l[i]) * pli[i];
            wr[i] = one + gamfac2 * (pmid[i] - states[PRS].r[i]) * pri[i];

            wl[i] = cl[i] * std::sqrt(wl[i]);
            wr[i] = cr[i] * std::sqrt(wr[i]);

            zl[i] = 4.0 * rli[i] * wl[i] * wl[i];
            zr[i] = 4.0 * rri[i] * wr[i] * wr[i];

            zl[i] = -zl[i] * wl[i] /
                    (zl[i] - gamfac1 * (pmid[i] - states[PRS].l[i]));
            zr[i] = zr[i] * wr[i] /
                    (zr[i] - gamfac1 * (pmid[i] - states[PRS].r[i]));

            umidl[i] = states[XVL].l[i] - (pmid[i] - states[PRS].l[i]) / wl[i];
            umidr[i] = states[XVL].r[i] + (pmid[i] - states[PRS].r[i]) / wr[i];

            pmid[i] +=
                (umidr[i] - umidl[i]) * (zl[i] * zr[i]) / (zr[i] - zl[i]);
            pmid[i] = std::max(small_riemann, pmid[i]);
            if (std::abs(pmid[i] - pmold[i]) / pmid[i] < NewRapTol)
              break;
          }
        }

        for (index_type i = bp::i0 + 3; i < bp::iN - 2; ++i) {
          umidl[i] = states[XVL].l[i] - (pmid[i] - states[PRS].l[i]) / wl[i];
          umidr[i] = states[XVL].r[i] + (pmid[i] - states[PRS].r[i]) / wr[i];
          umid[i] = half * (umidl[i] + umidr[i]);
        }
      }
    };

  }  // namespace ppm

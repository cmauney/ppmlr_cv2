#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <valarray>

#include "hydro_physics.hpp"
#include "ppm_constant.hpp"
#include "ppm_util.hpp"
#include "state.hpp"
#include "sweep_vector.hpp"

namespace ppm
{

    constexpr double NewRapTol = 1.0E-5;
    constexpr size_t NewRapItr = 12;

    template <typename SweepVector>
    struct riemann_t {
      using vector_type = SweepVector;
      using index_type = typename SweepVector::index_type;
      using state_set = state_set_t<vector_type>;

      double gamfac1, gamfac2;

      // solutions
      //  sweep_vector_t<NZ_SWEEP> pmid, umid;
      // newton-raphson solver vecs
      vector_type cl, cr;
      vector_type pli, pri, rli, rri;
      vector_type wl, wr, zl, zr;
      vector_type umidl, umidr, pmold;

      riemann_t() = default;

      inline void operator()(const hydro_physics_t& hydro,
                             const state_set& states,
                             vector_type& pmid,
                             vector_type& umid) {
        double gamfac1 = hydro.gamma + 1.0;
        double gamfac2 = 0.5 * gamfac1 / hydro.gamma;

        //        auto lo = NZ_GHOST - 3, hi = nx + 2 * NZ_GHOST - 2;
        for (index_type i = vector_type::i0() + 3; i < vector_type::iN() - 2;
             ++i) {
          //          cl[i] = std::sqrt(hydro.gamma * states[PRS].l[i] *
          //          states[RHO].l[i]);
          cl[i] = std::sqrt(hydro.gamma * states[PRS].l[i] * states[RHO].l[i]);
          pli[i] = 1.0 / states[PRS].l[i];
          //          pri[i] = 1.0 / states[PRS].r[i];
          rli[i] = 1.0 / states[RHO].l[i];
          //          rri[i] = 1.0 / states[RHO].r[i];

          //          pmid[i] = states[PRS].r[i] - states[PRS].l[i] -
          //                    cr[i] * (states[XVL].r[i] - states[XVL].l[i]);
          //          pmid[i] = states[PRS].l[i] + pmid[i] * cl[i] / (cl[i] +
          //          cr[i]); pmid[i] = std::max(small_riemann, pmid[i]);
        }

        for (index_type i = vector_type::i0() + 3; i < vector_type::iN() - 2;
             ++i) {
          cr[i] = std::sqrt(hydro.gamma * states[PRS].r[i] * states[RHO].r[i]);
          pri[i] = 1.0 / states[PRS].r[i];
          rri[i] = 1.0 / states[RHO].r[i];
        }

        for (index_type i = vector_type::i0() + 3; i < vector_type::iN() - 2;
             ++i) {
          pmid[i] = states[PRS].r[i] - states[PRS].l[i] -
                    cr[i] * (states[XVL].r[i] - states[XVL].l[i]);
        pmid[i] = states[PRS].l[i] + pmid[i] * cl[i] / (cl[i] + cr[i]);
        pmid[i] = std::max(small_riemann, pmid[i]);
        }

        /* auto dbmsg = "riemann solve " + std::to_string(__LINE__); */
        /* chknan(lo, hi, rli, "rli", dbmsg); */
        /* chknan(lo, hi, rri, "rri", dbmsg, true); */

        for (index_type i = vector_type::i0() + 3; i < vector_type::iN() - 2;
             ++i) {
          for (size_t n = 0; n < NewRapItr; ++n) {
            pmold[i] = pmid[i];
            wl[i] = 1.0 + gamfac2 * (pmid[i] - states[PRS].l[i]) * pli[i];
            wr[i] = 1.0 + gamfac2 * (pmid[i] - states[PRS].r[i]) * pri[i];

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

        for (index_type i = vector_type::i0() + 3; i < vector_type::iN() - 2;
             ++i) {
          umidl[i] = states[XVL].l[i] - (pmid[i] - states[PRS].l[i]) / wl[i];
          umidr[i] = states[XVL].r[i] + (pmid[i] - states[PRS].r[i]) / wr[i];
          umid[i] = 0.5 * (umidl[i] + umidr[i]);
        }

        /* chkzero(lo, hi, wl, "wl", "riemann solve"); */
        /* chknan(lo, hi, wl, "wl", "riemann solve"); */
        /* chknan(lo, hi, wr, "wr", "riemann solve"); */
        /* chknan(lo, hi, umid, "umid", "riemann solve"); */
      }
    };

  }  // namespace ppm

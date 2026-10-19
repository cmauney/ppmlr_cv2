#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include <Kokkos_Core.hpp>

#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"
#include "ppm_vector.hpp"

namespace ppm
{
template <index_type N>
struct parabola_basis_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using basis_type = Kokkos::View<double*>;

  // internal variables
  vector_type a, b, c, d, ai, bi, ci;

  // basis
  basis_type basis;

  // spans/views
  parabola_basis_t() : a(base::make_vector()), b(base::make_vector()),
                       c(base::make_vector()), d(base::make_vector()),
                       ai(base::make_vector()), bi(base::make_vector()),
                       ci(base::make_vector()), basis("basis", bp::iN * 5) {}
  ~parabola_basis_t() = default;

  constexpr inline void clear() {
    // Zero-fill using Kokkos deep_copy
    Kokkos::deep_copy(a, 0.0);
    Kokkos::deep_copy(b, 0.0);
    Kokkos::deep_copy(c, 0.0);
    Kokkos::deep_copy(d, 0.0);
    Kokkos::deep_copy(ai, 0.0);
    Kokkos::deep_copy(bi, 0.0);
    Kokkos::deep_copy(ci, 0.0);
    Kokkos::deep_copy(basis, 0.0);
  }

  constexpr inline void prepare(const vector_type& dx) {
    // Phase 1: Compute geometric factors (GPU-accelerated)
    Kokkos::parallel_for("basis_prep_phase1",
        Kokkos::RangePolicy<>(bp::i0, bp::iN - 1),
        KOKKOS_LAMBDA(const index_type i) {
          a(i) = dx(i) + dx(i + 1);
          b(i) = a(i) + dx(i);
          c(i) = a(i) + dx(i + 1);
          ai(i) = 1.0 / a(i);
          bi(i) = 1.0 / b(i);
          ci(i) = 1.0 / c(i);
        });

    // Phase 2: Compute basis coefficients (GPU-accelerated)
    Kokkos::parallel_for("basis_prep_phase2",
        Kokkos::RangePolicy<>(bp::i0 + 1, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          const size_t b_row = i * 5;
          d(i) = 1.0 / (a(i - 1) + a(i + 1));
          basis(b_row + 0) = dx(i) * ai(i) + 2.0 * dx(i + 1) * dx(i) * d(i) * ai(i) *
                             (a(i - 1) * bi(i) - a(i + 1) * ci(i));

          basis(b_row + 1) = -d(i) * dx(i) * a(i - 1) * bi(i);
          basis(b_row + 2) = d(i) * dx(i + 1) * a(i + 1) * ci(i);

          // Second basis set
          d(i) = dx(i) / (a(i - 1) + dx(i + 1));
          basis(b_row + 3) = d(i) * b(i - 1) * ai(i);
          basis(b_row + 4) = d(i) * c(i) * ai(i - 1);
        });
  }
};
template <index_type N>
struct parabola_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using basis_type = typename parabola_basis_t<N>::basis_type;

  // parabola coeffs
  vector_type da, a6, al, ar;

  // internals
  vector_type diff_an, del_a;
  vector_type two_diff_an_abs, two_diff_an_min;
  vector_type tmp1, tmp2, tmp3;

  parabola_t() : da(base::make_vector()), a6(base::make_vector()),
                 al(base::make_vector()), ar(base::make_vector()),
                 diff_an(base::make_vector()), del_a(base::make_vector()),
                 two_diff_an_abs(base::make_vector()), two_diff_an_min(base::make_vector()),
                 tmp1(base::make_vector()), tmp2(base::make_vector()),
                 tmp3(base::make_vector()) {}
  ~parabola_t() = default;

  constexpr inline void clear() {
    // Zero-fill using Kokkos deep_copy
    Kokkos::deep_copy(diff_an, 0.0);
    Kokkos::deep_copy(del_a, 0.0);
    Kokkos::deep_copy(tmp1, 0.0);
    Kokkos::deep_copy(tmp2, 0.0);
    Kokkos::deep_copy(tmp3, 0.0);
    Kokkos::deep_copy(da, 0.0);
    Kokkos::deep_copy(a6, 0.0);
    Kokkos::deep_copy(al, 0.0);
    Kokkos::deep_copy(ar, 0.0);
  }

  constexpr inline void construct(const vector_type& a,
                                  const vector_type& flat,
                                  const basis_type& pbase) {
    // Phase 1: Compute differences and absolute values (GPU-accelerated)
    Kokkos::parallel_for("para_diff",
        Kokkos::RangePolicy<>(bp::i0, bp::iN - 1),
        KOKKOS_LAMBDA(const index_type i) {
          diff_an(i) = a(i + 1) - a(i);
          two_diff_an_abs(i) = 2.0 * Kokkos::fabs(diff_an(i));
        });

    // Phase 2: Compute minimums (GPU-accelerated)
    Kokkos::parallel_for("para_min",
        Kokkos::RangePolicy<>(bp::i0 + 1, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          two_diff_an_min(i) = Kokkos::fmin(two_diff_an_abs(i - 1), two_diff_an_abs(i));
        });

    // Phase 3: Monotonicity constraint and delta computation (GPU-accelerated)
    Kokkos::parallel_for("para_delta",
        Kokkos::RangePolicy<>(bp::i0 + 1, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          const index_type b_row = i * 5;
          const double delta = pbase(b_row + 3) * diff_an(i) + pbase(b_row + 4) * diff_an(i - 1);

          // Monotonicity constraint with copysign
          double del = Kokkos::copysign(Kokkos::fmin(Kokkos::fabs(delta), two_diff_an_min(i)), delta);

          // Zero-crossing check (monotonicity) - branch acceptable on modern GPUs
          if (diff_an(i - 1) * diff_an(i) < 0.0)
            del = 0.0;

          del_a(i) = del;
        });

    // Phase 4: Edge state computation (GPU-accelerated)
    Kokkos::parallel_for("para_edges",
        Kokkos::RangePolicy<>(bp::i0 + 1, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          const index_type b_row = i * 5;
          const double ar_val = a(i) + pbase(b_row + 0) * diff_an(i) +
                                pbase(b_row + 1) * del_a(i + 1) + pbase(b_row + 2) * del_a(i);
          ar(i) = ar_val;
          al(i + 1) = ar_val;
        });

    // Phase 5: Flattening application (GPU-accelerated)
    Kokkos::parallel_for("para_flatten",
        Kokkos::RangePolicy<>(bp::i0 + 2, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          const double onemfl = 1.0 - flat(i);
          ar(i) = flat(i) * a(i) + onemfl * ar(i);
          al(i) = flat(i) * a(i) + onemfl * al(i);
        });

    // Phase 6: Compute parabola coefficients and temps (GPU-accelerated)
    Kokkos::parallel_for("para_coeffs",
        Kokkos::RangePolicy<>(bp::i0 + 2, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          const double da_val = ar(i) - al(i);
          const double a6_val = 6.0 * (a(i) - half * (al(i) + ar(i)));

          tmp1(i) = (ar(i) - a(i)) * (a(i) - al(i));
          tmp2(i) = da_val * da_val;
          tmp3(i) = da_val * a6_val;

          da(i) = da_val;
          a6(i) = a6_val;
        });

    // Phase 7: Apply monotonicity constraints (GPU-accelerated)
    // Branch divergence acceptable - ~10-20% overhead on modern GPUs
    Kokkos::parallel_for("para_monotonicity",
        Kokkos::RangePolicy<>(bp::i0 + 2, bp::iN - 2),
        KOKKOS_LAMBDA(const index_type i) {
          bool redo = false;
          double ar_new = ar(i);
          double al_new = al(i);

          if (tmp1(i) <= zero) {
            ar_new = a(i);
            al_new = a(i);
            redo = true;
          }
          if (tmp2(i) < tmp3(i)) {
            al_new = three * a(i) - two * ar_new;
            redo = true;
          }
          if (tmp2(i) < -tmp3(i)) {
            ar_new = three * a(i) - two * al_new;
            redo = true;
          }

          if (redo) {
            ar(i) = ar_new;
            al(i) = al_new;
            da(i) = ar_new - al_new;
            a6(i) = 6.0 * (a(i) - half * (al_new + ar_new));
          }
        });
  }

  constexpr inline double rc_left(const size_t i,
                                  const double f1,
                                  const double f2) const {
    return ar[i] - f1 * (da[i] - f2 * a6[i]);
  }

  constexpr inline double rc_right(const size_t i,
                                   const double f1,
                                   const double f2) const {
    return al[i] + f1 * (da[i] + f2 * a6[i]);
  }
};

template <typename ParabolaType, index_type N>
struct parabola_set_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using parabola_type = ParabolaType;
  using parabola_basis = parabola_basis_t<N>;

  parabola_basis basis;

  // Individual parabola instances for each fluid variable
  parabola_type para_rho;
  parabola_type para_prs;
  parabola_type para_xvl;
  parabola_type para_yvl;
  parabola_type para_zvl;
  parabola_type para_ent;
  parabola_type para_ein;
  parabola_type para_flt;

  parabola_set_t() = default;
  ~parabola_set_t() = default;

  inline void prepare(const vector_type& dx) { basis.prepare(dx); }

  template <size_t FN>
  constexpr inline void build(const vector_type& rho, const vector_type& prs,
                              const vector_type& xvl, const vector_type& yvl,
                              const vector_type& zvl, const vector_type& ent,
                              const vector_type& ein, const vector_type& flt,
                              const vector_type& flat,
                              const std::array<FLUID_VI, FN>& fidx) {
    for (FLUID_VI fi : fidx) {
      switch(fi) {
        case RHO: para_rho.construct(rho, flat, basis.basis); break;
        case PRS: para_prs.construct(prs, flat, basis.basis); break;
        case XVL: para_xvl.construct(xvl, flat, basis.basis); break;
        case YVL: para_yvl.construct(yvl, flat, basis.basis); break;
        case ZVL: para_zvl.construct(zvl, flat, basis.basis); break;
        case ENT: para_ent.construct(ent, flat, basis.basis); break;
        case EIN: para_ein.construct(ein, flat, basis.basis); break;
        case FLT: para_flt.construct(flt, flat, basis.basis); break;
      }
    }
  }

  // Helper to get parabola by enum (for compatibility, but prefer direct access)
  constexpr inline const parabola_type& get(FLUID_VI idx) const {
    switch(idx) {
      case RHO: return para_rho;
      case PRS: return para_prs;
      case XVL: return para_xvl;
      case YVL: return para_yvl;
      case ZVL: return para_zvl;
      case ENT: return para_ent;
      case EIN: return para_ein;
      case FLT: return para_flt;
    }
    return para_rho; // fallback
  }

  constexpr inline parabola_type& get(FLUID_VI idx) {
    switch(idx) {
      case RHO: return para_rho;
      case PRS: return para_prs;
      case XVL: return para_xvl;
      case YVL: return para_yvl;
      case ZVL: return para_zvl;
      case ENT: return para_ent;
      case EIN: return para_ein;
      case FLT: return para_flt;
    }
    return para_rho; // fallback
  }
};
}// namespace ppm

#pragma once

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "ppm_constant.hh"
#include "ppm_def.hh"
#include "ppm_grid.hh"
#include "ppm_vector.hh"

namespace ppm
{
template <index_type N>
struct parabola_basis_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  // using basis_type = fixed_storage_type<bp::iN * 5>;
  using basis_type = vector_wrapper_t<bp::iN * 5>;

  // internal variables
  vector_type a, b, c, d, ai, bi, ci;

  // basis:
  //  std::array<double, bp::iN * 5> basis;
  basis_type basis;

  // spans/views
  parabola_basis_t() = default;
  ~parabola_basis_t() = default;

  constexpr inline void clear() {
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    ai.clear();
    bi.clear();
    ci.clear();
    basis.clear();
  }

  constexpr inline void prepare(const vector_type& dx) {
    for (index_type i = bp::i0; i < bp::iN - 1; ++i) {
      a[i] = dx[i] + dx[i + 1];
      b[i] = a[i] + dx[i];
      c[i] = a[i] + dx[i + 1];
      ai[i] = 1. / a[i];
      bi[i] = 1. / b[i];
      ci[i] = 1. / c[i];
    }

    for (index_type i = bp::i0 + 1; i < bp::iN - 2; ++i) {
      size_t b_row = i * 5;
      d[i] = 1. / (a[i - 1] + a[i + 1]);
      basis[b_row + 0] =
          dx[i] * ai[i] + 2.0 * dx[i + 1] * dx[i] * d[i] * ai[i] *
                              (a[i - 1] * bi[i] - a[i + 1] * ci[i]);

      basis[b_row + 1] = -d[i] * dx[i] * a[i - 1] * bi[i];
      basis[b_row + 2] = d[i] * dx[i + 1] * a[i + 1] * ci[i];

      // check this range - upd: should be fine, one extra term worth not
      // restarting loop;
      d[i] = dx[i] / (a[i - 1] + dx[i + 1]);
      basis[b_row + 3] = d[i] * b[i - 1] * ai[i];
      basis[b_row + 4] = d[i] * c[i] * ai[i - 1];
    }
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

  parabola_t() = default;
  ~parabola_t() = default;

  constexpr inline void clear() {
    diff_an.clear();
    del_a.clear();

    tmp1.clear();
    tmp2.clear();
    tmp3.clear();

    da.clear();
    a6.clear();
    al.clear();
    ar.clear();
  }

  constexpr inline void construct(const vector_type& a,
                                  const vector_type& flat,
                                  const basis_type& pbase) {
    // clear();

    for (index_type i = bp::i0; i < bp::iN - 1; ++i) {
      diff_an[i] = a[i + 1] - a[i];
      two_diff_an_abs[i] = 2.0 * std::abs(diff_an[i]);
    }

    for (index_type i = bp::i0 + 1; i < bp::iN - 2; ++i) {
      two_diff_an_min[i] = std::min(two_diff_an_abs[i - 1], two_diff_an_abs[i]);
    }

    for (index_type i = bp::i0 + 1; i < bp::iN - 2; ++i) {
      index_type b_row = i * 5;
      del_a[i] =
          pbase[b_row + 3] * diff_an[i] + pbase[b_row + 4] * diff_an[i - 1];

      del_a[i] = std::copysign(
          //          std::min(std::abs(del_a[i]), std::min(2.0 *
          //          std::abs(diff_an[i - 1]),
          //                                                2.0 *
          //                                                std::abs(diff_an[i]))),
          std::min(std::abs(del_a[i]), two_diff_an_min[i]),

          del_a[i]);

      if (diff_an[i - 1] * diff_an[i] < 0.0)
        del_a[i] = 0.0;
    }

    for (index_type i = bp::i0 + 1; i < bp::iN - 2; ++i) {
      index_type b_row = i * 5;
      ar[i] = a[i] + pbase[b_row + 0] * diff_an[i] +
              pbase[b_row + 1] * del_a[i + 1] + pbase[b_row + 2] * del_a[i];
      al[i + 1] = ar[i];
    }

    for (index_type i = bp::i0 + 2; i < bp::iN - 2; ++i) {
      double onemfl = 1.0 - flat[i];
      ar[i] = flat[i] * a[i] + onemfl * ar[i];
      al[i] = flat[i] * a[i] + onemfl * al[i];
    }

    for (index_type i = bp::i0 + 2; i < bp::iN - 2; ++i) {
      da[i] = ar[i] - al[i];
      a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));

      tmp1[i] = (ar[i] - a[i]) * (a[i] - al[i]);
      tmp2[i] = da[i] * da[i];
      tmp3[i] = da[i] * a6[i];
    }

    for (index_type i = bp::i0 + 2; i < bp::iN - 2; ++i) {
      bool redo = false;
      if (tmp1[i] <= zero) {
        ar[i] = a[i];
        al[i] = a[i];
        redo = true;
      }
      if (tmp2[i] < +tmp3[i]) {
        al[i] = three * a[i] - two * ar[i];
        redo = true;
      }
      if (tmp2[i] < -tmp3[i]) {
        ar[i] = three * a[i] - two * al[i];
        redo = true;
      }
      if (redo) {
        da[i] = ar[i] - al[i];
        a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));
      }
    }

    /*
       for (index_type i = bp::i0 + 2; i < bp::iN - 2; ++i) {
          da[i] = ar[i] - al[i];
          a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));
        }
    */
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
  using parabola_array = FluidVars<parabola_type>;

  parabola_basis basis;
  parabola_array parabolas;

  parabola_set_t() = default;
  ~parabola_set_t() = default;

  inline void prepare(const vector_type& dx) { basis.prepare(dx); }

  template <size_t FN>
  constexpr inline void build(const FluidVars<vector_type>& fluid_v,
                              const vector_type& flat,
                              const std::array<FLUID_VI, FN>& fidx) {
    for (FLUID_VI fi : fidx) {
      parabolas[fi].construct(fluid_v[fi], flat, basis.basis);
    }
  }

  constexpr inline const parabola_type& operator[](const index_type idx) const {
    return parabolas[idx];
  }

  constexpr inline const parabola_array& arr() const { return parabolas; }
};
}// namespace ppm

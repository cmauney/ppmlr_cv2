#pragma once

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_util.hpp"
#include "sweep_grid.hpp"
#include "sweep_vector.hpp"

namespace ppm
{
template <typename PPMVector>
struct parabola_basis_t {
  using vector_type = PPMVector;
  using index_type = typename vector_type::index_type;
  using basis_type = std::array<double, vector_type::iN * 5>;

  // internal variables
  vector_type a, b, c, d, ai, bi, ci;

  // basis
  std::array<double, vector_type::iN() * 5> basis;

  // spans/views
  parabola_basis_t() = default;
  ~parabola_basis_t() = default;

  inline void clear() {
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    ai.clear();
    bi.clear();
    ci.clear();
    basis.clear();
  }

  inline void prepare(const vector_type& dx) {
    for (index_type i = vector_type::i0; i < vector_type::iN - 1; ++i) {
      a[i] = dx[i] + dx[i + 1];
      b[i] = a[i] + dx[i];
      c[i] = a[i] + dx[i + 1];
      ai[i] = 1. / a[i];
      bi[i] = 1. / b[i];
      ci[i] = 1. / c[i];
    }

    for (index_type i = vector_type::i0 + 1; i < vector_type::iN - 2; ++i) {
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
template <typename PPMVector>
struct parabola_t {
  using vector_type = PPMVector;
  using index_type = typename vector_type::index_type;
  using basis_type = typename parabola_basis_t<vector_type>::basis_type;

  // parabola coeffs
  vector_type da, a6, al, ar;

  // internals
  vector_type diff_an, del_a;
  vector_type tmp1, tmp2, tmp3;

  parabola_t() = default;
  ~parabola_t() = default;

  inline void clear() {
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

  inline void construct(const vector_type& a,
                        const vector_type& flat,
                        const basis_type& pbase) {
    clear();

    for (index_type i = vector_type::i0; i < vector_type::iN - 1; ++i) {
      diff_an[i] = a[i + 1] - a[i];
    }

    for (index_type i = vector_type::i0 + 1; i < vector_type::iN - 2; ++i) {
      size_t b_row = i * 5;
      del_a[i] =
          pbase[b_row + 3] * diff_an[i] + pbase[b_row + 4] * diff_an[i - 1];

      del_a[i] = std::copysign(
          std::min(std::abs(del_a[i]), std::min(2.0 * std::abs(diff_an[i - 1]),
                                                2.0 * std::abs(diff_an[i]))),
          del_a[i]);

      if (diff_an[i - 1] * diff_an[i] < 0.0)
        del_a[i] = 0.0;
    }

    for (index_type i = vector_type::i0 + 1; i < vector_type::iN - 2; ++i) {
      size_t b_row = i * 5;
      ar[i] = a[i] + pbase[b_row + 0] * diff_an[i] +
              pbase[b_row + 1] * del_a[i + 1] + pbase[b_row + 2] * del_a[i];
      al[i + 1] = ar[i];
    }

    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      double onemfl = 1.0 - flat[i];
      ar[i] = flat[i] * a[i] + onemfl * ar[i];
      al[i] = flat[i] * a[i] + onemfl * al[i];
    }

    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      da[i] = ar[i] - al[i];
      a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));

      tmp1[i] = (ar[i] - a[i]) * (a[i] - al[i]);
      tmp2[i] = da[i] * da[i];
      tmp3[i] = da[i] * a6[i];
    }
    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      if (tmp1[i] <= zero) {
        ar[i] = a[i];
        al[i] = a[i];
      }
      if (tmp2[i] < +tmp3[i])
        al[i] = three * a[i] - two * ar[i];
      if (tmp2[i] < -tmp3[i])
        ar[i] = three * a[i] - two * al[i];
    }

    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      da[i] = ar[i] - al[i];
      a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));
    }
  }

  inline double rc_left(const size_t i,
                        const double f1,
                        const double f2) const {
    return ar[i] - f1 * (da[i] - f2 * a6[i]);
  }

  inline double rc_right(const size_t i,
                         const double f1,
                         const double f2) const {
    return al[i] + f1 * (da[i] + f2 * a6[i]);
  }
};
template <typename ParabolaType>
struct parabola_set_t {
  using vector_type = typename ParabolaType::vector_type;
  using index_type = typename vector_type::index_type;

  using parabola = ParabolaType;
  using parabola_basis = typename parabola::basis_type;
  using parabola_array = FluidVars<parabola>;

  parabola_basis base;
  parabola_array parabolas;

  parabola_set_t() = default;
  ~parabola_set_t() = default;

  inline void prepare(const vector_type& dx) { base.prepare(dx); }

  template <index_type FN>
  constexpr inline void build(const FluidVars<vector_type>& fluid_v,
                              const vector_type& flat,
                              const std::array<FLUID_VI, FN>& fidx) {
    for (FLUID_VI fi : fidx) {
      parabolas[fi].construct(fluid_v[fi], flat, base.basis);
    }
  }

  constexpr inline const parabola& operator[](const index_type idx) const {
    return parabolas[idx];
  }
};
}// namespace ppm

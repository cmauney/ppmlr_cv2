#pragma once

#include <array>
#include <experimental/type_traits>
#include <iostream>
#include <utility>
#include <vector>

#include "gsl/gsl"

#include "ppm_def.hpp"

namespace ppm {

namespace detail {
template <class V>
using sum_t = decltype(std::declval<V&>().sum());

template <class V>
using max_t = decltype(std::declval<V&>().max());

template <typename V, size_t N>
inline constexpr bool needs_alloc =
    !(std::is_array_v<V> && std::is_same_v<V, std::array<double, N>>);
template <typename V>
inline constexpr bool has_sum = std::experimental::is_detected_v<sum_t, V>;

template <typename V>
inline constexpr bool has_max = std::experimental::is_detected_v<max_t, V>;
}  // namespace detail

template <typename VectorType, size_t NG, size_t NX>
struct __attribute__((aligned(16))) ppmlr_vector_t {
  using storage_type = VectorType;
  using index_type = std::ptrdiff_t;

  constexpr static index_type i0 = 0;
  constexpr static index_type iN = NX + 2 * NG;

  constexpr static index_type j0 = NG;
  constexpr static index_type jM = NG + NX;

  VectorType storage_;

  ppmlr_vector_t() {
    if constexpr (detail::needs_alloc<VectorType, iN>) {
      storage_.resize(iN);
    }
    clear();
  };

  ~ppmlr_vector_t() = default;

  constexpr inline double& operator[](const index_type idx) {
    return storage_[idx];
  };
  constexpr inline const double& operator[](const index_type idx) const {
    return storage_[idx];
  };


  // constexpr inline double inz(size_t idx) const;
  constexpr inline void clear() { fill(0.0); };
  constexpr inline void fill(const double value) {
    for (auto i = i0; i < iN; ++i)
      storage_[i] = value;
  };

  template <typename OtherVec>
  constexpr inline void assign(const OtherVec& v) {
    for (auto i = j0; i < jM; ++i)
      storage_[i] = v[i - j0];
  }

  constexpr inline double sum() const {
    if constexpr (detail::has_sum<VectorType>) {
      return storage_.sum();
    } else {
      double sum = 0.0;
      for (auto i = 0; i < iN; ++i) {
        sum += storage_[i];
      }
      return sum;
    }
  }
  constexpr inline double max() const {
    if constexpr (detail::has_max<VectorType>)
      return storage_.max();
    else {
      double maxv = std::numeric_limits<double>::lowest();
      for (auto i = 0; i < iN; ++i) {
        maxv = std::max(maxv, storage_[i]);
      }
      return maxv;
    }
  }
};

}  // namespace ppm


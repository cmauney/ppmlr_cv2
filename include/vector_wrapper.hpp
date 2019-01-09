#pragma once

#include <array>
#include <experimental/type_traits>
#include <iostream>
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include "global_def.hpp"

template <typename V, class Enable = void>
constexpr inline auto vsum(const V& v) {
  double sum = 0.0;
  for (index_type i = 0; i < v.size(); ++i) {
    sum += v[i];
  }
  return sum;
}

template <typename V, class Enable = void>
constexpr inline auto vmax(const V& v) {
  double maxv = std::numeric_limits<double>::lowest();
  for (index_type i = 0; i < v.size(); ++i) {
    maxv = std::max(maxv, v[i]);
  }
  return maxv;
}

template <typename V, typename std::enable_if_t<has_sum<V>>>
constexpr inline auto vsum(const V& v) {
  return v.sum();
}

template <typename V, typename std::enable_if_t<has_max<V>>>
constexpr inline auto vmax(const V& v) {
  return v.max();
}

template <index_type N, index_type B>
struct vector_frame_t {
  constexpr static index_type i0 = 0;
  constexpr static index_type iN = N + 2 * B;
  constexpr static index_type j0 = B;
  constexpr static index_type jM = N + B;
};

template <index_type N, index_type B = 0>
struct __attribute__((aligned(32))) vector_wrapper_t {
  using vf = vector_frame_t<N, B>;

  //  using vtype = std::array<double, vf::iN>;
  using vtype = std::valarray<double>;
  //  using vtype = std::vector<double>;

  //  fixed_storage_type<vf::iN> storage_;
  //  std::array<double, vf::iN> storage_;
  // std::vector<double> storage_;
  // std::valarray<double> storage_;
  vtype storage_;

  vector_wrapper_t() {
    //    storage_.resize(vf::iN);
    if constexpr (has_resize<vtype>)
      storage_.resize(vf::iN);

    clear();
  };

  vector_wrapper_t(decltype(storage_)&& v) : storage_(std::move(v)) {}

  ~vector_wrapper_t() = default;

  constexpr inline double& operator[](const index_type idx) {
    return storage_[idx];
  };
  constexpr inline const double& operator[](const index_type idx) const {
    return storage_[idx];
  };

  // constexpr inline double inz(size_t idx) const;
  constexpr inline void fill(const double value) {
    // storage_ = value;
    for (index_type i = vf::i0; i < vf::iN; ++i)
      storage_[i] = value;
  };

  constexpr inline void clear() { fill(0.0); };

  template <typename OtherVec>
  constexpr inline void assign(const OtherVec& v) {
    for (index_type i = vf::j0; i < vf::jM; ++i)
      storage_[i] = v[i - vf::j0];
  }

  constexpr inline double sum() const { return vsum(storage_); }

  constexpr inline double max() const { return vmax(storage_); }

  constexpr inline auto abs() const {
    // return vector_wrapper_t<N>(std::abs(storage_));
    vector_wrapper_t<N, B> rv;
    for (index_type i = vf::i0; i < vf::iN; ++i) {
      rv[i] = std::abs(storage_[i]);
        }
        return rv;
  }
};


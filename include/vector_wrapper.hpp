#pragma once

#include <array>
#include <experimental/type_traits>
#include <iostream>
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include "eigen3/Eigen/Dense"

#include "global_def.hpp"
#include "vector_adapt.hpp"

template <index_type N, index_type B>
struct vector_frame_t {
  constexpr static index_type i0 = 0;
  constexpr static index_type iN = N + 2 * B;
  constexpr static index_type j0 = B;
  constexpr static index_type jM = N + B;
};

template <index_type N, index_type B = 0>
// struct __attribute__((aligned(32))) vector_wrapper_t {
struct vector_wrapper_t {
  using vf = vector_frame_t<N, B>;

  // using vtype = double[vf::iN];
  //  using vtype = std::array<double, vf::iN>;
  using vtype = std::valarray<double>;
  // using vtype = std::vector<double>;
  // using vtype = Eigen::Array<double, vf::iN, 1>;

  //  fixed_storage_type<vf::iN> storage_;
  //  std::array<double, vf::iN> storage_;
  // std::vector<double> storage_;

  vtype storage_;

  vector_wrapper_t() {
    //    storage_.resize(vf::iN);
    if constexpr (vec_adt::has_resize<vtype>)
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

  constexpr inline double sum() const {
    return vec_adt::vsum(storage_, vf::i0, vf::iN);
  }

  constexpr inline double max() const {
    return vec_adt::vmax(storage_, vf::i0, vf::iN);
  }

  constexpr inline auto abs() const {
    return vec_adt::vabs(storage_, vf::i0, vf::iN);
  }
};


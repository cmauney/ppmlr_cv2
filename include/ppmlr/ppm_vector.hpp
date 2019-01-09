#pragma once

#include <array>
#include <experimental/type_traits>
#include <iostream>
#include <utility>
#include <vector>

#include "ppm_def.hpp"

namespace ppm {
/*
template <index_type N>
struct __attribute__((aligned(16))) ppmlr_vector_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  fixed_storage_type<bp::iN> storage_;

  ppmlr_vector_t() { clear(); };

  ~ppmlr_vector_t() = default;

  constexpr inline double& operator[](const index_type idx) {
    return storage_[idx];
  };
  constexpr inline const double& operator[](const index_type idx) const {
    return storage_[idx];
  };


  // constexpr inline double inz(size_t idx) const;
  constexpr inline void fill(const double value) {
    for (index_type i = bp::i0; i < bp::iN; ++i)
      storage_[i] = value;
  };

  constexpr inline void clear() { fill(0.0); };

  template <typename OtherVec>
  constexpr inline void assign(const OtherVec& v) {
    for (index_type i = bp::j0; i < bp::jM; ++i)
      storage_[i] = v[i - bp::j0];
  }

  constexpr inline double sum() const {
    double sum = 0.0;
    for (auto i = 0; i < bp::iN; ++i) {
      sum += storage_[i];
    }
    return sum;
  }
  constexpr inline double max() const {
      double maxv = std::numeric_limits<double>::lowest();
      for (index_type i = bp::i0; i < bp::iN; ++i) {
        maxv = std::max(maxv, storage_[i]);
      }
      return maxv;
  }
};
*/
template <typename Container>
using FluidVars = std::array<Container, N_FLUIDV>;


}  // namespace ppm


#pragma once

#include <array>
#include <experimental/type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include "eigen3/Eigen/Dense"

#include "global_def.hpp"

namespace vec_adt {

#define OP_CHECK(op, op_sym)                                                   \
  template <typename V>                                                        \
  using op##_t = decltype(std::declval<V&>() op_sym std::declval<const V&>()); \
  template <typename V>                                                        \
  inline constexpr bool has_##op =                                             \
      std::experimental::is_detected_v<op##_t, V>;

#define METHOD_CHK(fn, ...)                                             \
  template <typename V, typename... Vs>                                 \
  using fn##_t = decltype(std::declval<V>().fn(std::declval<Vs>()...)); \
  template <typename V>                                                 \
  inline constexpr bool has_##fn =                                      \
      std::experimental::is_detected_v<fn##_t, V, ##__VA_ARGS__>;

#define ARGUMENT_CHK(fn, ...)                                            \
  template <typename V, typename... Vs>                                  \
  using fn##_t = decltype(fn(std::declval<V>(), std::declval<Vs>()...)); \
  template <typename V>                                                  \
  inline constexpr bool has_##fn =                                       \
      std::experimental::is_detected_v<fn##_t, V, ##__VA_ARGS__>;

METHOD_CHK(sum)
METHOD_CHK(max)
METHOD_CHK(maxCoeff)  // for Eigen
METHOD_CHK(resize, int)

ARGUMENT_CHK(abs)

OP_CHECK(plus, +)
OP_CHECK(minus, -)
OP_CHECK(times, *)
OP_CHECK(divide, /)

/*
 * compile-time assert some test cases
 * can be disregarded if things running smoothly
 */
static_assert(!has_resize<std::array<double, 1>>,
              "array should not have resize");
static_assert(has_resize<std::vector<double>>, "vector should have resize");
static_assert(has_resize<std::valarray<double>>, "valarray should have resize");

static_assert(!has_sum<std::array<double, 1>>, "array should not have sum");
static_assert(!has_sum<std::vector<double>>, "vector should not have sum");
static_assert(has_sum<std::valarray<double>>, "valarray should have sum");

static_assert(!has_max<std::array<double, 1>>, "array should not have max");
static_assert(!has_max<std::vector<double>>, "vector should not have max");
static_assert(has_max<std::valarray<double>>, "valarray should have max");

static_assert(has_plus<std::valarray<double>>,
              "valarray should have operator+");
static_assert(has_divide<std::valarray<double>>,
              "valarray should have operator/");
static_assert(!has_plus<std::array<double, 1>>, "array has no operator+");
static_assert(!has_divide<std::array<double, 1>>, "array has no operator/");

static_assert(!has_abs<std::vector<double>>, "abs does not work on vectors");
static_assert(has_abs<std::valarray<double>>, "abs should work on valarray");

static_assert(has_maxCoeff<Eigen::Array<double, 1, 1>>,
              "eigen should have maxCoeff");
static_assert(has_sum<Eigen::Array<double, 1, 1>>, "eigen should have sum");

/*
 * using SFINAE to generate the correct vector functions.
 */

/*template <typename V, class Enable = void>
constexpr inline auto vsum(const V& v, index_type start, index_type size) {
  double sum = 0.0;
  const double* p = &v[start + size];
  while (p-- != &v[start]) {
    sum += *p;
  }

  return sum;
}
*/
template <typename V>
typename std::enable_if_t<!has_sum<V>, double> constexpr inline vsum(
    const V& v,
    index_type start,
    index_type size) {
  double sum = 0.0;
  const double* p = &v[start + size];
  while (p-- != &v[start]) {
    sum += *p;
  }
  return sum;
}

template <typename V>
typename std::enable_if_t<has_sum<V>, double> constexpr inline vsum(
    const V& v,
    index_type,
    index_type) {
  return v.sum();
}

/*
template <typename V, class Enable = void>
constexpr inline auto vmax(const V& v, index_type start, index_type size) {
  double maxv = std::numeric_limits<double>::lowest();
  const double* p = &v[start + size];

  while (p-- != &v[start]) {
    maxv = std::max(maxv, *p);
  }
  return maxv;
}
*/
template <typename V>
typename std::enable_if_t<!has_maxCoeff<V> && !has_max<V>,
                          double> constexpr inline vmax(const V& v,
                                                        index_type start,
                                                        index_type size) {
  double maxv = std::numeric_limits<double>::lowest();
  const double* p = &v[start + size];

  while (p-- != &v[start]) {
    maxv = std::max(maxv, *p);
  }
  return maxv;
}
/*template <typename V, typename std::enable_if_t<has_maxCoeff<V>>>
constexpr inline auto vmax(const V& v, index_type = -1, index_type = -1) {
  return v.maxCoeff();
}*/
template <typename V>
typename std::enable_if_t<has_maxCoeff<V>, double> constexpr inline vmax(
    const V& v,
    index_type = -1,
    index_type = -1) {
  return v.maxCoeff();
}

template <typename V>
typename std::enable_if_t<has_max<V>, double> constexpr inline vmax(
    const V& v,
    index_type = -1,
    index_type = -1) {
  return v.max();
}

template <typename V>
typename std::enable_if_t<!has_abs<V>, V> constexpr inline vabs(
    const V& v,
    index_type start,
    index_type size) {
  V rv(v);
  double* p = &rv[start + size];

  while (p-- != &v[start]) {
    *p = std::abs(*p);
  }
  return std::move(rv);
}

template <typename V>
typename std::enable_if_t<has_abs<V>, V> constexpr inline vabs(
    const V& v,
    index_type = -1,
    index_type = -1) {
  return abs(v);
}

/*
template <typename V, typename std::enable_if_t<has_max<V>>>
constexpr inline auto vmax(const V& v, index_type = -1, index_type = -1) {
  return v.max();
}
*/
/*
template <typename V, class Enable = void>
constexpr inline auto vbas(const V& v, index_type start, index_type size) {
  V rv(v);
  double* p = &rv[start + size];

  while (p-- != &v[start]) {
    *p = std::abs(*p);
  }
  return rv;
}

template <typename V, typename std::enable_if_t<has_abs<V>>>
constexpr inline auto vabs(const V& v, index_type = -1, index_type = -1) {
  return abs(v);
}
*/
}  // namespace vec_adt

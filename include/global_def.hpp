#pragma once

#include <array>
#include <experimental/type_traits>
#include <utility>
#include <valarray>
#include <vector>

// namespace detail {
template <typename V, typename... Vs>
using resize_t = decltype(std::declval<V>().resize(std::declval<Vs>()...));
template <typename V>
inline constexpr bool has_resize =
    std::experimental::is_detected_v<resize_t, V, size_t>;

template <typename V>
using sum_t = decltype(std::declval<V&>().sum());
template <typename V>
inline constexpr bool has_sum = std::experimental::is_detected_v<sum_t, V>;

template <typename V>
using max_t = decltype(std::declval<V&>().max());
template <typename V>
inline constexpr bool has_max = std::experimental::is_detected_v<max_t, V>;
//}  // namespace detail

static_assert(!has_resize<std::array<double, 1>>, "its broke");
static_assert(has_resize<std::vector<double>>, "its broke");
static_assert(has_resize<std::valarray<double>>, "its broke");

static_assert(!has_sum<std::array<double, 1>>, "its broke");
static_assert(!has_sum<std::vector<double>>, "its broke");
static_assert(!has_max<std::vector<double>>, "its broke");
static_assert(has_max<std::valarray<double>>, "its broke");
static_assert(has_sum<std::valarray<double>>, "its broke");

/*#define METHOD_CHK(f)                            \
  template <class V>                             \
  using f##_t = decltype(std::declval<V>().f()); \
  template <typename V>                          \
  inline constexpr bool has_##f = std::experimental::is_detected_v<f##_t, V>;
*/
// METHOD_CHK(resize)
// METHOD_CHK(sum)
// METHOD_CHK(max)
// index type
using index_type = size_t;

// vector storage type
// template <index_type N>
// using fixed_storage_type = std::array<double, N>;

// fluid variable IDs
enum FLUID_VI { RHO, PRS, XVL, YVL, ZVL, ENT, EIN, FLT };

// IDs of the primative variables (technically prs is not primative but easier
// to use. these define the states of the riemann problem during ppm
constexpr std::array<FLUID_VI, 3> PRIMATIVE_VARS{RHO, PRS, XVL};

// conserved quantities. note these are 'specific', i.e. per-unit-mass.
constexpr std::array<FLUID_VI, 6> CONSERVATIVE_VARS{RHO, XVL, YVL,
                                                    ZVL, EIN, ENT};

// total number of advected quantities
constexpr index_type N_FLUIDV = 8;

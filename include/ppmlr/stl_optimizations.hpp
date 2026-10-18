#pragma once

#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>

// STL algorithm optimizations for PPM computations
// Provides optimized implementations using parallel execution policies where beneficial

namespace ppm {
namespace stl_opt {

// Execution policy selection based on problem size
// Small problems: sequential (avoid threading overhead)
// Large problems: parallel (leverage multiple cores)
constexpr size_t PARALLEL_THRESHOLD = 1000;

template<typename T>
inline auto get_execution_policy(size_t size) {
  if constexpr (std::is_same_v<T, std::execution::parallel_policy>) {
    return (size > PARALLEL_THRESHOLD) ? std::execution::par : std::execution::seq;
  }
  return std::execution::seq;
}

// Optimized vector operations using STL algorithms

// Transform with binary operation (optimized for SIMD)
template<typename InputIt1, typename InputIt2, typename OutputIt, typename BinaryOp>
inline void transform_binary(InputIt1 first1, InputIt1 last1,
                             InputIt2 first2, OutputIt d_first,
                             BinaryOp op, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    std::transform(std::execution::par_unseq, first1, last1, first2, d_first, op);
  } else {
    std::transform(first1, last1, first2, d_first, op);
  }
}

// Transform with unary operation
template<typename InputIt, typename OutputIt, typename UnaryOp>
inline void transform_unary(InputIt first, InputIt last, OutputIt d_first,
                            UnaryOp op, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    std::transform(std::execution::par_unseq, first, last, d_first, op);
  } else {
    std::transform(first, last, d_first, op);
  }
}

// Reduction operation (sum, max, min)
template<typename InputIt, typename T, typename BinaryOp>
inline T reduce(InputIt first, InputIt last, T init, BinaryOp op, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    return std::reduce(std::execution::par, first, last, init, op);
  } else {
    return std::reduce(first, last, init, op);
  }
}

// Transform-reduce (fused operation for better performance)
template<typename InputIt1, typename InputIt2, typename T>
inline T inner_product(InputIt1 first1, InputIt1 last1,
                       InputIt2 first2, T init, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    return std::transform_reduce(std::execution::par, first1, last1, first2, init);
  } else {
    return std::inner_product(first1, last1, first2, init);
  }
}

// For_each with execution policy
template<typename InputIt, typename UnaryFunction>
inline void for_each(InputIt first, InputIt last, UnaryFunction f, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    std::for_each(std::execution::par_unseq, first, last, f);
  } else {
    std::for_each(first, last, f);
  }
}

// Optimized min/max element finding
template<typename ForwardIt>
inline ForwardIt max_element(ForwardIt first, ForwardIt last, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    return std::max_element(std::execution::par_unseq, first, last);
  } else {
    return std::max_element(first, last);
  }
}

template<typename ForwardIt>
inline ForwardIt min_element(ForwardIt first, ForwardIt last, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    return std::min_element(std::execution::par_unseq, first, last);
  } else {
    return std::min_element(first, last);
  }
}

// Optimized fill operations
template<typename ForwardIt, typename T>
inline void fill(ForwardIt first, ForwardIt last, const T& value, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    std::fill(std::execution::par_unseq, first, last, value);
  } else {
    std::fill(first, last, value);
  }
}

// Optimized copy operations
template<typename InputIt, typename OutputIt>
inline OutputIt copy(InputIt first, InputIt last, OutputIt d_first, size_t hint_size = 0) {
  if (hint_size > PARALLEL_THRESHOLD) {
    return std::copy(std::execution::par_unseq, first, last, d_first);
  } else {
    return std::copy(first, last, d_first);
  }
}

// Fused multiply-add operation (optimized for FMA instructions)
inline void fused_multiply_add(const double* __restrict__ a,
                               const double* __restrict__ b,
                               const double* __restrict__ c,
                               double* __restrict__ result,
                               size_t n) {
  // a * b + c, optimized for FMA
  for (size_t i = 0; i < n; ++i) {
    result[i] = std::fma(a[i], b[i], c[i]);
  }
}

// Vector magnitude squared (avoid sqrt when possible)
inline double magnitude_squared(const double* __restrict__ x,
                               const double* __restrict__ y,
                               const double* __restrict__ z,
                               size_t i) {
  return x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
}

// Clamp to range [min, max]
template<typename T>
inline constexpr T clamp(const T& value, const T& min_val, const T& max_val) {
  return std::clamp(value, min_val, max_val);
}

// Linear interpolation (lerp) - C++20
template<typename T>
inline constexpr T lerp(const T& a, const T& b, const T& t) {
  return std::lerp(a, b, t);
}

// Vectorized operations using STL algorithms
namespace vectorized {

// Compute: result = a * b (element-wise multiply)
inline void multiply(const double* __restrict__ a,
                    const double* __restrict__ b,
                    double* __restrict__ result,
                    size_t n) {
  std::transform(a, a + n, b, result, std::multiplies<>());
}

// Compute: result = a + b (element-wise add)
inline void add(const double* __restrict__ a,
               const double* __restrict__ b,
               double* __restrict__ result,
               size_t n) {
  std::transform(a, a + n, b, result, std::plus<>());
}

// Compute: result = a - b (element-wise subtract)
inline void subtract(const double* __restrict__ a,
                    const double* __restrict__ b,
                    double* __restrict__ result,
                    size_t n) {
  std::transform(a, a + n, b, result, std::minus<>());
}

// Compute: result = a / b (element-wise divide)
inline void divide(const double* __restrict__ a,
                  const double* __restrict__ b,
                  double* __restrict__ result,
                  size_t n) {
  std::transform(a, a + n, b, result, std::divides<>());
}

// Compute: result = scalar * a (scalar multiply)
inline void scale(const double* __restrict__ a,
                 double scalar,
                 double* __restrict__ result,
                 size_t n) {
  std::transform(a, a + n, result, [scalar](double x) { return scalar * x; });
}

// Compute: result = sqrt(a) (element-wise sqrt)
inline void sqrt_vec(const double* __restrict__ a,
                    double* __restrict__ result,
                    size_t n) {
  std::transform(a, a + n, result, [](double x) { return std::sqrt(x); });
}

// Compute: result = abs(a) (element-wise absolute value)
inline void abs_vec(const double* __restrict__ a,
                   double* __restrict__ result,
                   size_t n) {
  std::transform(a, a + n, result, [](double x) { return std::abs(x); });
}

// Compute: result = max(a, b) (element-wise maximum)
inline void max_vec(const double* __restrict__ a,
                   const double* __restrict__ b,
                   double* __restrict__ result,
                   size_t n) {
  std::transform(a, a + n, b, result,
                [](double x, double y) { return std::max(x, y); });
}

// Compute: result = min(a, b) (element-wise minimum)
inline void min_vec(const double* __restrict__ a,
                   const double* __restrict__ b,
                   double* __restrict__ result,
                   size_t n) {
  std::transform(a, a + n, b, result,
                [](double x, double y) { return std::min(x, y); });
}

} // namespace vectorized

// Cache-friendly operations
namespace cache_opt {

// Block size for cache-friendly operations (typically L1 cache line)
constexpr size_t BLOCK_SIZE = 64 / sizeof(double); // 8 doubles per cache line

// Prefetch hint for next iteration
template<typename T>
inline void prefetch(const T* ptr) {
  __builtin_prefetch(ptr, 0, 3); // Read, high temporal locality
}

// Blocked loop for better cache utilization
template<typename Func>
inline void blocked_loop(size_t start, size_t end, Func&& func) {
  for (size_t block = start; block < end; block += BLOCK_SIZE) {
    size_t block_end = std::min(block + BLOCK_SIZE, end);
    for (size_t i = block; i < block_end; ++i) {
      func(i);
    }
  }
}

} // namespace cache_opt

} // namespace stl_opt
} // namespace ppm

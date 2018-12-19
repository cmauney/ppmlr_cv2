#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

namespace ppm {

template <typename V>
inline void print_and_deal(const size_t lo,
                           const size_t hi,
                           const size_t i,
                           const V& v,
                           std::string var,
                           std::string msg,
                           bool deal = false) {
  std::cout << "i = " << i << "\n";
  std::cout << var << " [ " << i << " ] = " << v[i] << "\n";
  if (i > lo)
    std::cout << var << " [ " << i - 1 << " ] = " << v[i - 1] << std::endl;
  if (i < hi - 1)
    std::cout << var << " [ " << i + 1 << " ] = " << v[i + 1] << std::endl;

  std::cout << msg << std::endl;
  if (deal)
    std::exit(1);
}

template <typename V>
inline void chknan(const size_t lo,
                   const size_t hi,
                   const V& v,
                   std::string var,
                   std::string msg,
                   bool deal = false) {
  for (size_t i = lo; i < hi; ++i) {
    if (std::isnan(v[i])) {
      std::cout << "NAN " << var << std::endl;
      print_and_deal(lo, hi, i, v, var, msg, deal);
    }
  }
}

template <typename V>
inline void chkzero(const size_t lo,
                    const size_t hi,
                    const V& v,
                    std::string var,
                    std::string msg,
                    bool deal = false) {
  for (size_t i = lo; i < hi; ++i) {
    if (v[i] == 0.0) {
      std::cout << "BAD0 " << var << std::endl;
      print_and_deal(lo, hi, i, v, var, msg, deal);
    }
  }
}

template <typename V>
inline void chkneg(const size_t lo,
                   const size_t hi,
                   const V& v,
                   std::string var,
                   std::string msg,
                   bool deal = false) {
  for (size_t i = lo; i < hi; ++i) {
    if (v[i] < 0.0) {
      std::cout << " NEG " << var << std::endl;
      print_and_deal(lo, hi, i, v, var, msg, deal);
    }
  }
}

template <typename I, typename... Vs>
inline void dump_vec_and_keep_calm(const I lo, const I hi, Vs&&... vs) {
  for (I i = lo; i < hi; ++i) {
    auto s = std::to_string(i) + " " + std::to_string(i + 1);
    s += (" " + ... + (" " + std::to_string(vs[i])));
    std::cout << s << std::endl;
  }
  std::exit(1);
}

}  // namespace ppm


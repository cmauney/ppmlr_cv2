#pragma once

#include <algorithm>
#include <execution>
#include <cmath>
#include "ppm_constant.hpp"
#include "ppm_def.hpp"

// Computational kernels prepared for device offloading
// These functions are structured to be easily parallelized on GPU
// with OpenMP target directives, CUDA, or other offloading technologies
// Optimized with STL algorithms and parallel execution policies where beneficial

namespace ppm {
namespace kernels {

// Use parallel execution for large arrays
constexpr size_t PARALLEL_THRESHOLD = 500;

// Kinetic energy computation kernel
// OFFLOAD TARGET: This loop is embarrassingly parallel
// OPTIMIZED: Manual loop with restrict for best vectorization
inline void compute_kinetic_energy(
    const double* __restrict__ xvl,
    const double* __restrict__ yvl,
    const double* __restrict__ zvl,
    double* __restrict__ ekin,
    const index_type i0,
    const index_type iN) {

  // Manual loop with restrict pointers allows compiler to:
  // - Auto-vectorize with SIMD (AVX, NEON, etc.)
  // - Apply FMA (fused multiply-add) optimizations
  // - Eliminate redundant loads/stores
  // Perfect for GPU offload: independent iterations, regular access pattern
  // Future: #pragma omp target teams distribute parallel for simd
  #pragma GCC ivdep  // Tell compiler loop iterations are independent
  for (index_type i = i0; i < iN; ++i) {
    const double vx = xvl[i];
    const double vy = yvl[i];
    const double vz = zvl[i];
    ekin[i] = 0.5 * (vx * vx + vy * vy + vz * vz);
  }
}

// Sound speed and CFL computation kernel
// OFFLOAD TARGET: Parallel reduction for max_cdx
inline double compute_sound_speed(
    const double* __restrict__ prs,
    const double* __restrict__ rho,
    const double* __restrict__ dx,
    double* __restrict__ cdx,
    const double gamma,
    const double radius,
    const index_type i0,
    const index_type iN) {

  double max_cdx = 0.0;

  // Future: #pragma omp target teams distribute parallel for simd reduction(max:max_cdx)
  for (index_type i = i0; i < iN; ++i) {
    cdx[i] = std::sqrt(gamma * prs[i] / rho[i]) / (dx[i] * radius);
    max_cdx = std::max(max_cdx, cdx[i]);
  }

  return max_cdx;
}

// Entropy computation kernel
// OFFLOAD TARGET: Embarrassingly parallel
inline void compute_entropy(
    const double* __restrict__ prs,
    const double* __restrict__ rho,
    const double* __restrict__ ekin,
    double* __restrict__ ent,
    const double gammam1,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    ent[i] = prs[i] / (rho[i] * gammam1) + ekin[i];
  }
}

// Density remap kernel
// OFFLOAD TARGET: Independent density updates
inline void remap_density(
    double* __restrict__ rho,
    const double* __restrict__ flux_rho,
    const double* __restrict__ dvol0,
    const double* __restrict__ dvol1,
    double* __restrict__ dm,
    double* __restrict__ dm0,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    dm[i] = rho[i] * dvol1[i];
    dm0[i] = dm[i] + flux_rho[i] - flux_rho[i + 1];
    rho[i] = dm0[i] / dvol0[i];
    rho[i] = std::max(smallr, rho[i]);
    dm0[i] = 1.0 / (rho[i] * dvol0[i]);
  }
}

// Velocity/scalar remap kernel
// OFFLOAD TARGET: Independent scalar updates
inline void remap_scalar(
    double* __restrict__ var,
    const double* __restrict__ flux_var,
    const double* __restrict__ dm,
    const double* __restrict__ dm0,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    var[i] = (var[i] * dm[i] + flux_var[i] - flux_var[i + 1]) * dm0[i];
  }
}

// Evolve density kernel (Lagrangian step)
// OFFLOAD TARGET: Independent density evolution
inline void evolve_density(
    double* __restrict__ rho,
    const double* __restrict__ dvol0,
    const double* __restrict__ dvol1,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    rho[i] *= (dvol0[i] / dvol1[i]);
  }
}

// Evolve velocity kernel (momentum update)
// OFFLOAD TARGET: Independent velocity updates
inline void evolve_velocity(
    double* __restrict__ xvl,
    double* __restrict__ uold,
    const double* __restrict__ pmid,
    const double* __restrict__ amid,
    const double* __restrict__ dtdm,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    uold[i] = xvl[i];
    xvl[i] -= dtdm[i] * (pmid[i + 1] - pmid[i]) * 0.5 * (amid[i + 1] + amid[i]);
  }
}

// Evolve energy kernel
// OFFLOAD TARGET: Independent energy updates
inline void evolve_energy(
    double* __restrict__ ent,
    const double* __restrict__ amid,
    const double* __restrict__ upmid,
    const double* __restrict__ dtdm,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    ent[i] -= dtdm[i] * (amid[i + 1] * upmid[i + 1] - amid[i] * upmid[i]);
  }
}

// Internal energy from total energy kernel
// OFFLOAD TARGET: Independent energy decomposition
inline void compute_internal_energy(
    double* __restrict__ ein,
    const double* __restrict__ ent,
    const double* __restrict__ ekin,
    const double* __restrict__ rho,
    const double gammam1,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    ein[i] = ent[i] - ekin[i];
    ein[i] = std::max(ein[i], smallp / (gammam1 * rho[i]));
  }
}

// Pressure from internal energy kernel
// OFFLOAD TARGET: Equation of state evaluation
inline void compute_pressure(
    double* __restrict__ prs,
    const double* __restrict__ rho,
    double* __restrict__ ein,
    const double* __restrict__ ent,
    const double* __restrict__ ekin,
    const double gammam1,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    // Supersonic flow check
    if (ekin[i] / ein[i] < 100.0) {
      ein[i] = ent[i] - ekin[i];
    }
    prs[i] = gammam1 * rho[i] * ein[i];
    prs[i] = std::max(smallp, prs[i]);
  }
}

// Grid volume and mass kernel
// OFFLOAD TARGET: Independent grid calculations
inline void compute_grid_mass(
    double* __restrict__ dm,
    double* __restrict__ dtdm,
    const double* __restrict__ rho,
    const double* __restrict__ dvol,
    const double dt,
    const index_type i0,
    const index_type iN) {

  // Future: #pragma omp target teams distribute parallel for simd
  for (index_type i = i0; i < iN; ++i) {
    dm[i] = rho[i] * dvol[i];
    dtdm[i] = dt / dm[i];
  }
}

} // namespace kernels
} // namespace ppm

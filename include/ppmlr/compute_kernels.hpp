#pragma once

#include <Kokkos_Core.hpp>
#include "ppm_constant.hpp"
#include "ppm_def.hpp"

// Computational kernels for GPU offloading with Kokkos
// All kernels use Kokkos::parallel_for or Kokkos::parallel_reduce
// for portable performance across CPU (Serial/OpenMP) and GPU (CUDA/HIP/SYCL)

namespace ppm {
namespace kernels {

// Kinetic energy computation kernel
// GPU-accelerated: Embarrassingly parallel, coalesced memory access
KOKKOS_INLINE_FUNCTION
void compute_kinetic_energy(
    const double* __restrict__ xvl,
    const double* __restrict__ yvl,
    const double* __restrict__ zvl,
    double* __restrict__ ekin,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("kinetic_energy",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        const double vx = xvl[i];
        const double vy = yvl[i];
        const double vz = zvl[i];
        ekin[i] = 0.5 * (vx * vx + vy * vy + vz * vz);
      });
}

// Sound speed and CFL computation kernel with parallel reduction
// GPU-accelerated: Parallel reduction for maximum value
KOKKOS_INLINE_FUNCTION
double compute_sound_speed(
    const double* __restrict__ prs,
    const double* __restrict__ rho,
    const double* __restrict__ dx,
    double* __restrict__ cdx,
    const double gamma,
    const double radius,
    const index_type i0,
    const index_type iN) {

  double max_cdx = 0.0;

  Kokkos::parallel_reduce("sound_speed",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i, double& local_max) {
        const double c = Kokkos::sqrt(gamma * prs[i] / rho[i]) / (dx[i] * radius);
        cdx[i] = c;
        local_max = Kokkos::fmax(local_max, c);
      }, Kokkos::Max<double>(max_cdx));

  return max_cdx;
}

// Entropy computation kernel
// GPU-accelerated: Embarrassingly parallel
KOKKOS_INLINE_FUNCTION
void compute_entropy(
    const double* __restrict__ prs,
    const double* __restrict__ rho,
    const double* __restrict__ ekin,
    double* __restrict__ ent,
    const double gammam1,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("compute_entropy",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        ent[i] = prs[i] / (rho[i] * gammam1) + ekin[i];
      });
}

// Density remap kernel
// GPU-accelerated: Conservative remap with independent updates
KOKKOS_INLINE_FUNCTION
void remap_density(
    double* __restrict__ rho,
    const double* __restrict__ flux_rho,
    const double* __restrict__ dvol0,
    const double* __restrict__ dvol1,
    double* __restrict__ dm,
    double* __restrict__ dm0,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("remap_density",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        dm[i] = rho[i] * dvol1[i];
        dm0[i] = dm[i] + flux_rho[i] - flux_rho[i + 1];
        rho[i] = dm0[i] / dvol0[i];
        rho[i] = Kokkos::fmax(smallr, rho[i]);
        dm0[i] = 1.0 / (rho[i] * dvol0[i]);
      });
}

// Velocity/scalar remap kernel
// GPU-accelerated: Independent scalar advection
KOKKOS_INLINE_FUNCTION
void remap_scalar(
    double* __restrict__ var,
    const double* __restrict__ flux_var,
    const double* __restrict__ dm,
    const double* __restrict__ dm0,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("remap_scalar",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        var[i] = (var[i] * dm[i] + flux_var[i] - flux_var[i + 1]) * dm0[i];
      });
}

// Evolve density kernel (Lagrangian step)
// GPU-accelerated: Volume ratio scaling
KOKKOS_INLINE_FUNCTION
void evolve_density(
    double* __restrict__ rho,
    const double* __restrict__ dvol0,
    const double* __restrict__ dvol1,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("evolve_density",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        rho[i] *= (dvol0[i] / dvol1[i]);
      });
}

// Evolve velocity kernel (momentum update)
// GPU-accelerated: Pressure gradient force
KOKKOS_INLINE_FUNCTION
void evolve_velocity(
    double* __restrict__ xvl,
    double* __restrict__ uold,
    const double* __restrict__ pmid,
    const double* __restrict__ amid,
    const double* __restrict__ dtdm,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("evolve_velocity",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        uold[i] = xvl[i];
        xvl[i] -= dtdm[i] * (pmid[i + 1] - pmid[i]) * 0.5 * (amid[i + 1] + amid[i]);
      });
}

// Evolve energy kernel
// GPU-accelerated: PdV work term
KOKKOS_INLINE_FUNCTION
void evolve_energy(
    double* __restrict__ ent,
    const double* __restrict__ amid,
    const double* __restrict__ upmid,
    const double* __restrict__ dtdm,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("evolve_energy",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        ent[i] -= dtdm[i] * (amid[i + 1] * upmid[i + 1] - amid[i] * upmid[i]);
      });
}

// Internal energy from total energy kernel
// GPU-accelerated: Energy decomposition with floor
KOKKOS_INLINE_FUNCTION
void compute_internal_energy(
    double* __restrict__ ein,
    const double* __restrict__ ent,
    const double* __restrict__ ekin,
    const double* __restrict__ rho,
    const double gammam1,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("compute_internal_energy",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        ein[i] = ent[i] - ekin[i];
        ein[i] = Kokkos::fmax(ein[i], smallp / (gammam1 * rho[i]));
      });
}

// Pressure from internal energy kernel
// GPU-accelerated: Equation of state with supersonic check
KOKKOS_INLINE_FUNCTION
void compute_pressure(
    double* __restrict__ prs,
    const double* __restrict__ rho,
    double* __restrict__ ein,
    const double* __restrict__ ent,
    const double* __restrict__ ekin,
    const double gammam1,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("compute_pressure",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        // Supersonic flow check - branch acceptable on modern GPUs
        if (ekin[i] / ein[i] < 100.0) {
          ein[i] = ent[i] - ekin[i];
        }
        prs[i] = gammam1 * rho[i] * ein[i];
        prs[i] = Kokkos::fmax(smallp, prs[i]);
      });
}

// Grid volume and mass kernel
// GPU-accelerated: Mass and timestep factor computation
KOKKOS_INLINE_FUNCTION
void compute_grid_mass(
    double* __restrict__ dm,
    double* __restrict__ dtdm,
    const double* __restrict__ rho,
    const double* __restrict__ dvol,
    const double dt,
    const index_type i0,
    const index_type iN) {

  Kokkos::parallel_for("compute_grid_mass",
      Kokkos::RangePolicy<>(i0, iN),
      KOKKOS_LAMBDA(const index_type i) {
        dm[i] = rho[i] * dvol[i];
        dtdm[i] = dt / dm[i];
      });
}

} // namespace kernels
} // namespace ppm

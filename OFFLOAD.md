# Device Offloading Guide

This codebase has been structured to facilitate GPU/accelerator offloading. The computational kernels have been extracted into separate functions with clean data dependencies, making them ready for parallelization.

## Architecture

### 1. Compute Kernels (`include/ppmlr/compute_kernels.hpp`)

All computationally intensive loops have been extracted into standalone kernels in the `ppm::kernels` namespace. These kernels:

- **Use raw pointers** (`double* __restrict__`) for optimal device performance
- **Have clear loop bounds** (i0, iN) passed as parameters
- **Are embarrassingly parallel** - no cross-iteration dependencies
- **Use restrict keyword** to enable compiler vectorization/optimization
- **Are marked with comments** indicating they're offload targets

### 2. Fluid Data Structure (`include/ppmlr/fluid_data.hpp`)

The `fluid_data_t` template provides:

- **Contiguous memory layout** for all fluid variables
- **Raw pointer access** via `.data()` method for device transfers
- **Enum-based indexing** to maintain type safety
- **Ready for managed memory** or explicit device allocation

### 3. Modified Main Solver (`include/ppmlr/ppmlr.hpp`)

The main PPM solver now:

- **Calls kernel functions** instead of inline loops
- **Passes raw pointers** to kernels for device compatibility
- **Groups operations** into offloadable units
- **Separates host-only logic** (grid updates, boundary conditions) from device logic

## Offloading Strategies

### Option 1: OpenMP Target Directives (Recommended)

The easiest path to GPU offloading. Add target directives to kernel loops:

```cpp
// In compute_kernels.hpp, modify compute_kinetic_energy:
inline void compute_kinetic_energy(...) {
  #pragma omp target teams distribute parallel for simd \
      map(to: xvl[i0:iN-i0], yvl[i0:iN-i0], zvl[i0:iN-i0]) \
      map(from: ekin[i0:iN-i0])
  for (index_type i = i0; i < iN; ++i) {
    ekin[i] = 0.5 * (xvl[i] * xvl[i] + yvl[i] * yvl[i] + zvl[i] * zvl[i]);
  }
}
```

**Compile flags:**
```bash
# NVIDIA GPUs
cmake -DCMAKE_CXX_FLAGS="-fopenmp -fopenmp-targets=nvptx64" ..

# AMD GPUs  
cmake -DCMAKE_CXX_FLAGS="-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa" ..

# Intel GPUs
cmake -DCMAKE_CXX_FLAGS="-fopenmp -fopenmp-targets=spir64" ..
```

### Option 2: Kokkos Abstraction Layer

For portable performance across multiple architectures:

1. Wrap kernels in Kokkos parallel_for
2. Use Kokkos::View for data management
3. Build with appropriate backend (CUDA, HIP, SYCL, OpenMP)

```cpp
Kokkos::parallel_for("kinetic_energy", 
    Kokkos::RangePolicy<>(i0, iN),
    KOKKOS_LAMBDA(const int i) {
      ekin[i] = 0.5 * (xvl[i]*xvl[i] + yvl[i]*yvl[i] + zvl[i]*zvl[i]);
    });
```

### Option 3: Direct CUDA/HIP

For maximum control:

1. Convert kernels to CUDA `__global__` functions
2. Replace function calls with kernel launches
3. Manage device memory explicitly with cudaMalloc/cudaMemcpy

```cpp
__global__ void compute_kinetic_energy_kernel(
    const double* __restrict__ xvl,
    const double* __restrict__ yvl, 
    const double* __restrict__ zvl,
    double* __restrict__ ekin,
    const int N) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    ekin[i] = 0.5 * (xvl[i]*xvl[i] + yvl[i]*yvl[i] + zvl[i]*zvl[i]);
  }
}

// Launch: compute_kinetic_energy_kernel<<<blocks, threads>>>(...)
```

## Performance Optimization Tips

### 1. Data Locality

**Current structure supports:**
- All fluid variables are contiguous in memory
- Use `.data()` to get raw pointers for GPU transfer
- Consider pinned/managed memory for automatic transfers

**Optimization:**
```cpp
// Allocate device memory once at initialization
cudaMalloc(&d_rho, N * sizeof(double));
cudaMalloc(&d_prs, N * sizeof(double));
// ... etc

// Transfer at beginning of timestep
cudaMemcpy(d_rho, rho.data(), N * sizeof(double), cudaMemcpyHostToDevice);

// Compute on device
// ...

// Transfer back at end of timestep  
cudaMemcpy(rho.data(), d_rho, N * sizeof(double), cudaMemcpyDeviceToHost);
```

### 2. Kernel Fusion

The current structure allows easy fusion of related kernels to reduce transfers:

```cpp
// Instead of:
compute_kinetic_energy(...);
cudaDeviceSynchronize();
compute_entropy(...);

// Fuse into single kernel:
__global__ void compute_kinetic_and_entropy(...) {
  // Compute ekin and ent in one pass
}
```

### 3. Asynchronous Execution

Use streams to overlap computation and communication:

```cpp
cudaStream_t stream1, stream2;
cudaStreamCreate(&stream1);
cudaStreamCreate(&stream2);

// Overlap H2D transfer of next timestep with computation
cudaMemcpyAsync(d_data_next, h_data_next, size, 
                cudaMemcpyHostToDevice, stream1);
compute_kernel<<<grid, block, 0, stream2>>>(...);
```

### 4. Loop Bounds

All kernels use explicit loop bounds (i0, iN) rather than entire array:
- Reduces transfer sizes (only active zones)
- Avoids ghost cell computations on device
- Maps naturally to GPU thread blocks

### 5. Reduction Operations

For operations with reductions (e.g., `compute_sound_speed` with max):

**OpenMP:**
```cpp
#pragma omp target teams distribute parallel for simd reduction(max:max_cdx)
for (index_type i = i0; i < iN; ++i) {
  cdx[i] = ...;
  max_cdx = std::max(max_cdx, cdx[i]);
}
```

**CUDA:**
Use CUB or Thrust for efficient reductions:
```cpp
#include <thrust/reduce.h>
max_cdx = thrust::reduce(thrust::device, d_cdx, d_cdx + N, 
                         0.0, thrust::maximum<double>());
```

## Kernels Ready for Offload

| Kernel | Location | Parallelism | Data Dependencies |
|--------|----------|-------------|-------------------|
| `compute_kinetic_energy` | compute_kernels.hpp:20 | Embarrassingly parallel | None |
| `compute_sound_speed` | compute_kernels.hpp:33 | Parallel + reduction | None |
| `compute_entropy` | compute_kernels.hpp:52 | Embarrassingly parallel | None |
| `remap_density` | compute_kernels.hpp:67 | Embarrassingly parallel | Flux stencil (i±1) |
| `remap_scalar` | compute_kernels.hpp:85 | Embarrassingly parallel | Flux stencil (i±1) |
| `evolve_density` | compute_kernels.hpp:97 | Embarrassingly parallel | None |
| `evolve_velocity` | compute_kernels.hpp:109 | Embarrassingly parallel | Pmid stencil (i±1) |
| `evolve_energy` | compute_kernels.hpp:123 | Embarrassingly parallel | Amid stencil (i±1) |
| `compute_internal_energy` | compute_kernels.hpp:136 | Embarrassingly parallel | None |
| `compute_pressure` | compute_kernels.hpp:151 | Embarrassingly parallel | None |
| `compute_grid_mass` | compute_kernels.hpp:169 | Embarrassingly parallel | None |

## Current Limitations

1. **Parabola construction** - Complex with internal state, stays on host
2. **Riemann solver** - Conditional logic, may not benefit from offload
3. **Boundary conditions** - Small, host-only
4. **Grid operations** - Sequential dependencies, host-only

## Next Steps to Enable Offloading

### Phase 1: Validation (Current State)
✅ Kernels extracted and tested on CPU
✅ Data access patterns cleaned up
✅ Loop bounds parameterized

### Phase 2: Add OpenMP Target
1. Add `#pragma omp target` directives to kernels
2. Update CMakeLists.txt with offload flags
3. Test on target GPU architecture
4. Profile and optimize data transfers

### Phase 3: Optimize
1. Measure kernel timings
2. Fuse kernels where beneficial
3. Use asynchronous transfers
4. Consider managed memory

### Phase 4: Multi-GPU (Optional)
1. Domain decomposition
2. MPI + GPU for each rank
3. Overlap communication with computation

## Example: Enabling OpenMP Offload

**1. Update CMakeLists.txt:**
```cmake
option(ENABLE_GPU_OFFLOAD "Enable GPU offloading via OpenMP" OFF)

if(ENABLE_GPU_OFFLOAD)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        # NVIDIA
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=nvptx64")
        # AMD
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=amdgcn-amd-amdhsa")
    endif()
endif()
```

**2. Update compute_kernels.hpp:**
```cpp
#ifdef _OPENMP
  #pragma omp target teams distribute parallel for simd \
      map(to: xvl[i0:iN-i0], yvl[i0:iN-i0], zvl[i0:iN-i0]) \
      map(from: ekin[i0:iN-i0])
#endif
  for (index_type i = i0; i < iN; ++i) {
    ekin[i] = 0.5 * (xvl[i]*xvl[i] + yvl[i]*yvl[i] + zvl[i]*zvl[i]);
  }
```

**3. Build and run:**
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_GPU_OFFLOAD=ON ..
cmake --build .
./bin/ppm1d
```

## Performance Expectations

For a typical shock tube problem (N=250):
- **CPU only**: ~6-8 µs per timestep
- **GPU offload (naive)**: May be slower due to transfer overhead
- **GPU offload (optimized)**: 2-3x speedup for N>1000
- **GPU multi-physics (3D, N>10000)**: 10-100x speedup

**Note**: Small 1D problems don't benefit from GPU. Real gains come with:
- 2D/3D simulations (N > 10^6 zones)
- Multiple species/physics
- Long time integration

## References

- OpenMP 5.0 Target Directives: https://www.openmp.org/spec-html/5.0/openmpsu60.html
- Kokkos Programming Guide: https://kokkos.github.io/kokkos-core-wiki/
- CUDA C++ Programming Guide: https://docs.nvidia.com/cuda/cuda-c-programming-guide/
- HIP Programming Guide: https://rocmdocs.amd.com/en/latest/Programming_Guides/HIP-GUIDE.html

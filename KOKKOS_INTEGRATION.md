# Kokkos Integration Guide

## Overview
Kokkos has been integrated as a **required** dependency (git submodule) to provide performance portability across different hardware architectures (CPUs, GPUs, etc.). The runtime is automatically initialized and finalized in the main function.

## Submodule Location
- **Path**: `external/kokkos`
- **Repository**: https://github.com/kokkos/kokkos.git
- **Version**: 5.2.99 (development)

## CMake Configuration

### Basic Usage
Kokkos is **required** and enabled by default with the Serial backend:
```bash
cmake -B build
cmake --build build
```

**Note:** Kokkos cannot be disabled - it is a required dependency of this project.

### Backend Selection
Enable different backends by setting CMake options:

#### Serial Backend (Default)
```bash
cmake -B build -DKokkos_ENABLE_SERIAL=ON
```

#### OpenMP Backend
```bash
cmake -B build -DKokkos_ENABLE_OPENMP=ON
```

#### CUDA Backend (NVIDIA GPUs)
```bash
cmake -B build -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON
```
Replace `AMPERE80` with your GPU architecture:
- `VOLTA70` for V100
- `AMPERE80` for A100
- `ADA89` for RTX 4090
- `HOPPER90` for H100

#### HIP Backend (AMD GPUs)
```bash
cmake -B build -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_VEGA90A=ON
```
Replace with your AMD architecture as needed.

#### SYCL Backend (Intel GPUs)
```bash
cmake -B build -DKokkos_ENABLE_SYCL=ON
```

### Architecture Optimization
Enable native CPU architecture optimizations:
```bash
cmake -B build -DKokkos_ARCH_NATIVE=ON
```

## Environment Notes

### Cray Systems
On Cray systems, you may see a warning about `CRAYPE_LINK_TYPE`. If linking fails, set:
```bash
export CRAYPE_LINK_TYPE=dynamic
```

## Cloning the Repository with Submodules

### New Clone
When cloning this repository, initialize submodules:
```bash
git clone <repository-url>
cd ppmlr_cv2
git submodule update --init --recursive
```

### Updating Kokkos
To update Kokkos to the latest version:
```bash
cd external/kokkos
git pull origin develop  # or 'master' for stable releases
cd ../..
git add external/kokkos
git commit -m "Update Kokkos submodule"
```

## Build Configuration Summary

Current configuration shows:
- **Kokkos Version**: 5.2.99
- **Default Backend**: SERIAL
- **C++ Standard**: C++20 (inherited from project)
- **Compiler**: GCC 12.3.0

## Integration in Code

Kokkos has been integrated into `src/main.cpp`:
- `Kokkos::initialize(argc, argv)` is called at the start of main
- All application code runs within a Kokkos scope
- `Kokkos::finalize()` is called before program exit

The integration allows access to command-line options:
```bash
./ppm1d --kokkos-help              # Show Kokkos options
./ppm1d --kokkos-print-configuration  # Show runtime configuration
./ppm1d --kokkos-num-threads=8     # Set thread count
```

To add Kokkos parallel constructs in your code:
```cpp
#include <Kokkos_Core.hpp>

// In your function
Kokkos::parallel_for("MyKernel", N, KOKKOS_LAMBDA(int i) {
    // kernel body
});
```

## Next Steps

1. **Choose appropriate backend**: Decide which execution space fits your hardware
2. **Set architecture flags**: Optimize for your specific CPU/GPU
3. **Port code**: Begin migrating performance-critical loops to Kokkos parallel constructs
4. **Benchmark**: Compare performance across different backends

## References
- [Kokkos Documentation](https://kokkos.org/kokkos-core-wiki/)
- [Kokkos Tutorials](https://github.com/kokkos/kokkos-tutorials)
- [Programming Guide](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/ProgrammingModel.html)

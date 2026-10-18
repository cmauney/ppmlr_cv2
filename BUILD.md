# Build Instructions

This project uses CMake as the build system.

## Quick Start

### Debug Build (default)
```bash
mkdir -p build && cd build
cmake ..
cmake --build .
```

### Optimized Build (equivalent to old 'opt' target)
```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
cmake --build .
```

### Fast/Release Build (equivalent to old 'fast' target)
```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

## Build Targets

- **Debug**: Full debug symbols, sanitizers enabled, no optimization
- **RelWithDebInfo**: `-O2` optimization with debug info (good for profiling)
- **Release**: `-O3` with aggressive optimizations and loop transformations

## Advanced Options

### Enable Advanced Loop Optimizations
If your system has ISL library support, you can enable advanced loop optimizations:
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_LOOP_OPTS=ON ..
```

### Parallel Build
```bash
cmake --build . -j$(nproc)
```

### Clean Build
```bash
rm -rf build
mkdir build && cd build
cmake ..
cmake --build .
```

### Install
```bash
cmake --build . --target install
```

### Change Compiler
```bash
cmake -DCMAKE_CXX_COMPILER=clang++ ..
```

## Output

The executable `ppm1d` will be located in `build/bin/ppm1d`.

## Migration from Makefile

| Old Command | New Command |
|------------|-------------|
| `make` or `make debug` | `cmake -DCMAKE_BUILD_TYPE=Debug ..` |
| `make opt` | `cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..` |
| `make fast` | `cmake -DCMAKE_BUILD_TYPE=Release ..` |
| `make clean` | `rm -rf build` |

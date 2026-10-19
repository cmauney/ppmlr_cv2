cmake_minimum_required(VERSION 3.15)
project(ppm1d VERSION 1.0 LANGUAGES CXX)

# Set C++ standard (GCC 12.3 has partial C++23, using C++20 for full support)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Kokkos configuration (required)
# Kokkos backend selection
option(Kokkos_ENABLE_SERIAL "Enable Kokkos Serial backend" ON)
option(Kokkos_ENABLE_OPENMP "Enable Kokkos OpenMP backend" OFF)
option(Kokkos_ENABLE_CUDA "Enable Kokkos CUDA backend" OFF)
option(Kokkos_ENABLE_HIP "Enable Kokkos HIP backend for AMD GPUs" OFF)
option(Kokkos_ENABLE_SYCL "Enable Kokkos SYCL backend for Intel GPUs" OFF)

# Kokkos architecture flags (uncomment as needed)
# set(Kokkos_ARCH_NATIVE ON CACHE BOOL "Enable native CPU architecture")
# For NVIDIA: set(Kokkos_ARCH_AMPERE80 ON) or Kokkos_ARCH_VOLTA70, etc.
# For AMD: set(Kokkos_ARCH_VEGA90A ON) or appropriate architecture

# Add Kokkos subdirectory (required)
add_subdirectory(external/kokkos)
message(STATUS "Kokkos enabled (required)")

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Build type configuration
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Common compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color")

# LTO flags (optional, commented out by default as in Makefile)
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# OpenMP support (optional)
option(ENABLE_OPENMP "Enable OpenMP threading" OFF)
option(ENABLE_GPU_OFFLOAD "Enable GPU offloading via OpenMP target" OFF)

if(ENABLE_OPENMP OR ENABLE_GPU_OFFLOAD)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        message(STATUS "OpenMP enabled")

        if(ENABLE_GPU_OFFLOAD)
            # GPU offload target - uncomment appropriate line for your hardware
            # NVIDIA GPUs:
            # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=nvptx64")
            # AMD GPUs:
            # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=amdgcn-amd-amdhsa")
            # Intel GPUs:
            # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=spir64")
            message(STATUS "GPU offloading flags NOT set - edit CMakeLists.txt to uncomment target architecture")
        endif()
    endif()
endif()

# Advanced loop optimization flags (requires ISL library)
# These flags require Graphite support (-fgraphite, -floop-*), which needs ISL
# Disable them if ISL is not available on your system
option(ENABLE_LOOP_OPTS "Enable advanced loop optimizations (requires ISL)" OFF)
if(ENABLE_LOOP_OPTS)
    set(ILS_FLAGS "-ftree-vectorize -floop-interchange -floop-strip-mine -floop-block -ftree-loop-distribution -ftree-loop-distribute-patterns")
else()
    set(ILS_FLAGS "-ftree-vectorize")
endif()

# Debug build flags
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -Og -ggdb -fmax-errors=5 -fno-omit-frame-pointer -fno-inline")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize=alignment -fstack-check -fsanitize=undefined -fsanitize=float-divide-by-zero")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address -fsanitize=alignment -fstack-check -fsanitize=undefined -fsanitize=float-divide-by-zero")

# Release with debug info (equivalent to 'opt' target)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -march=native -Winline -g")

# Release build flags (equivalent to 'fast' target)
# Added: -ffast-math for aggressive FP optimizations (safe for this application)
#        -funroll-loops for loop unrolling
#        -fno-math-errno for faster math functions
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native ${ILS_FLAGS} -Winline -ffast-math -funroll-loops -fno-math-errno")

# Collect source files (excluding main.cpp for library)
file(GLOB LIBRARY_SOURCES
    ${CMAKE_SOURCE_DIR}/src/print.cpp
)

# Create library target for ppmlr
add_library(ppmlr ${LIBRARY_SOURCES})

# Add alias for modern CMake usage
add_library(ppmlr::ppmlr ALIAS ppmlr)

# Target-specific include directories
target_include_directories(ppmlr
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Target-specific link libraries
target_link_libraries(ppmlr
    PUBLIC
        Kokkos::kokkos
    PRIVATE
        m
)

# Create executable
add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)

# Link executable against ppmlr library
target_link_libraries(${PROJECT_NAME} PRIVATE ppmlr::ppmlr)

message(STATUS "Created ppmlr::ppmlr library target")
message(STATUS "Linking ${PROJECT_NAME} with ppmlr::ppmlr")

# Installation rules (optional)
install(TARGETS ${PROJECT_NAME} ppmlr
    EXPORT ppmlrTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

# Install export targets
install(EXPORT ppmlrTargets
    FILE ppmlrTargets.cmake
    NAMESPACE ppmlr::
    DESTINATION lib/cmake/ppmlr
)

# Print build configuration
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug Flags: ${CMAKE_CXX_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Release Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(STATUS "RelWithDebInfo Flags: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

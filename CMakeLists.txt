cmake_minimum_required(VERSION 3.15)
project(ppm1d VERSION 1.0 LANGUAGES CXX)

# Set C++ standard (GCC 12.3 has partial C++23, using C++20 for full support)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# OpenMP (optional, commented out by default as in Makefile)
# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()

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
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native ${ILS_FLAGS} -Winline")

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/gsl
    $ENV{HOME}/.local/include
    /usr/include
    /usr/local/include
)

# Collect all source files
file(GLOB_RECURSE SOURCES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE m)

# Installation rules (optional)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

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

project(ppmlr_cxx)
cmake_minimum_required(VERSION 3.30)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_language(CXX)
include(CMakeDetermineCXXCompiler)
include(CheckCXXCompilerFlag)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_subdirectory(external/fmt)
add_subdirectory(external/taskflow)
add_subdirectory(external/mdspan)

add_subdirectory(ppmlr)
add_subdirectory(app)

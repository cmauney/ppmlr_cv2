
add_library(ppmlr INTERFACE)
add_library(ppmlr::ppmlr ALIAS ppmlr)

target_compile_features(ppmlr INTERFACE cxx_std_23)

target_include_directories(ppmlr
INTERFACE
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ppmlr INTERFACE
  mdspan::mdspan
  Taskflow::Taskflow
  fmt::fmt
)

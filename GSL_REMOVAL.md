# GSL Library Removal

## Summary
The Microsoft GSL (Guidelines Support Library) has been completely removed from this project.

## Rationale
1. **Not Used**: GSL was not actively used in the codebase (only a commented-out include existed)
2. **Conflict with Kokkos**: GSL's `span` implementation conflicted with Kokkos's mdspan library
3. **C++20 Standard**: C++20 provides `std::span` natively, making GSL's span unnecessary
4. **Maintenance**: Removes an unused dependency, simplifying the build

## Changes Made

### Files Removed
- `include/gsl/gsl`
- `include/gsl/gsl_algorithm`
- `include/gsl/gsl_assert`
- `include/gsl/gsl_byte`
- `include/gsl/gsl_util`
- `include/gsl/multi_span`
- `include/gsl/pointers`
- `include/gsl/span`
- `include/gsl/string_span`

### Files Modified
1. **CMakeLists.txt**
   - Removed `${CMAKE_SOURCE_DIR}/include/gsl` from include_directories
   - Removed related comments about GSL

2. **include/ppmlr/boundary_conditions.hpp**
   - Removed commented-out line: `//#include "gsl/gsl"`

3. **src/main.cpp**
   - No changes needed (GSL was never included here)

## Verification
- ✓ Clean Debug build successful
- ✓ Clean Release build successful
- ✓ No compilation errors
- ✓ Runtime execution verified
- ✓ No GSL references remain in source code

## Migration Notes
If you need span functionality:
- Use `std::span` from `<span>` (C++20 standard)
- Kokkos provides its own view types for multi-dimensional arrays
- Modern C++20 features are available as replacements for most GSL utilities

## Date
October 19, 2026

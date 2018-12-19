#include <fstream>
#include <iostream>
#include <netcdf>
#include <string>
#include <valarray>

#include "ppm_constant.hpp"
#include "ppm_def.hpp"

template <typename V, typename... Outs>
void hdf5_out(const V& xc, const V& yc, Outs&&... os) {
  constexpr size_t N_VAR = sizeof...(os);
}


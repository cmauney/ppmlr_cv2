#pragma once

#include <vector>
#include "ppm_def.hpp"

namespace ppm {

// Device-friendly fluid data structure
// Groups all fluid variables with contiguous memory layout
// for efficient host-device transfers
template <index_type N>
struct fluid_data_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  // Primary fluid variables - all contiguous in memory
  vector_type rho;   // density
  vector_type prs;   // pressure
  vector_type xvl;   // x-velocity
  vector_type yvl;   // y-velocity
  vector_type zvl;   // z-velocity
  vector_type ent;   // total energy (entropy)
  vector_type ein;   // internal energy
  vector_type flt;   // flattening parameter

  // Constructor - initialize all vectors
  fluid_data_t()
      : rho(base::make_vector()), prs(base::make_vector()),
        xvl(base::make_vector()), yvl(base::make_vector()),
        zvl(base::make_vector()), ent(base::make_vector()),
        ein(base::make_vector()), flt(base::make_vector()) {}

  // Get raw pointer to specific variable for device offload
  inline double* data(FLUID_VI var) {
    switch(var) {
      case RHO: return rho.data();
      case PRS: return prs.data();
      case XVL: return xvl.data();
      case YVL: return yvl.data();
      case ZVL: return zvl.data();
      case ENT: return ent.data();
      case EIN: return ein.data();
      case FLT: return flt.data();
      default: return nullptr;
    }
  }

  // Get const raw pointer
  inline const double* data(FLUID_VI var) const {
    switch(var) {
      case RHO: return rho.data();
      case PRS: return prs.data();
      case XVL: return xvl.data();
      case YVL: return yvl.data();
      case ZVL: return zvl.data();
      case ENT: return ent.data();
      case EIN: return ein.data();
      case FLT: return flt.data();
      default: return nullptr;
    }
  }

  // Size of each vector
  inline size_t size() const { return rho.size(); }

  // Get reference to vector by enum
  inline vector_type& operator[](FLUID_VI var) {
    switch(var) {
      case RHO: return rho;
      case PRS: return prs;
      case XVL: return xvl;
      case YVL: return yvl;
      case ZVL: return zvl;
      case ENT: return ent;
      case EIN: return ein;
      case FLT: return flt;
    }
    return rho; // fallback
  }

  inline const vector_type& operator[](FLUID_VI var) const {
    switch(var) {
      case RHO: return rho;
      case PRS: return prs;
      case XVL: return xvl;
      case YVL: return yvl;
      case ZVL: return zvl;
      case ENT: return ent;
      case EIN: return ein;
      case FLT: return flt;
    }
    return rho; // fallback
  }
};

} // namespace ppm

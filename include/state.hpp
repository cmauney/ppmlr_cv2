#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <valarray>

#include "hydro_physics.hpp"
#include "parabola.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"
#include "ppm_util.hpp"

namespace ppm
{
template <typename PPMVector>
struct state_domain_t {
  using vector_type = PPMVector;
  using index_type = typename PPMVector::index_type;

  vector_type Cdtdx, fCdtdx;

  state_domain_t() = default;

  inline void build(const hydro_physics_t& hydro, const vector_type& cdx) {
    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      Cdtdx[i] = cdx[i] * hydro.hdt;
      fCdtdx[i] = 1.0 - fourthird * Cdtdx[i];
    }
  }
};

template <typename ReconstructType, typename DomainType>
struct state_t {
  using domain_type = DomainType;
  using reconstruct_type = ReconstructType;

  using vector_type = typename domain_type::vector_type;
  using index_type = typename vector_type::index_type;

  vector_type l, r;
  state_t() = default;

  inline void build_state(const reconstruct_type& reco,
                          const domain_type& domain) {
    for (index_type i = vector_type::i0 + 2; i < vector_type::iN - 2; ++i) {
      l[i + 1] = reco.rc_left(i, domain.Cdtdx[i], domain.fCdtdx[i]);
      r[i] = reco.rc_right(i, domain.Cdtdx[i], domain.fCdtdx[i]);
    }
  }
};

template <typename StateType, size_t NSt>
struct state_set_t {
  using state_type = StateType;

  using vector_type = typename state_type::vector_type;
  using index_type = typename vector_type::index_type;
  using reconstruct_type = typename state_type::reconstruct_type;
  using domain_type = typename state_type::domain_type;

  domain_type domain;
  std::array<state_type, NSt> states;
  std::array<FLUID_VI, NSt> fluid_ids;

  state_set_t() = default;

  state_set_t(std::array<FLUID_VI, NSt> rcis) : fluid_ids(std::move(rcis)) {}

  inline void prepare(const hydro_physics_t& hydro, const vector_type& cdx) {
    domain.build(hydro, cdx);
  }

  template <size_t Sp>
  inline void build(const std::array<reconstruct_type, Sp>& recos) {
    for (index_type i = 0; i < NSt; ++i) {
      states[i].build_state(recos[fluid_ids[i]], domain);
    }
  }

  inline const state_type& operator[](const FLUID_VI fi) const {
    for (index_type i = 0; i < NSt; ++i) {
      if (fi == fluid_ids[i])
        return states[i];
    }
  }
};

  }  // namespace ppm

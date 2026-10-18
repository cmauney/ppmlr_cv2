#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

#include "hydro_physics.hpp"
#include "parabola.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"

namespace ppm
{
template <index_type N>
struct state_domain_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  vector_type Cdtdx, fCdtdx;

  state_domain_t() : Cdtdx(base::make_vector()), fCdtdx(base::make_vector()) {}

  constexpr inline void build(const hydro::hydro_physics_t& hydro,
                              const vector_type& cdx) {
    for (index_type i = bp::i0 + 2; i < bp::iN - 2; ++i) {
      Cdtdx[i] = cdx[i] * hydro.hdt;
      fCdtdx[i] = 1.0 - fourthird * Cdtdx[i];
    }
  }
};

template <typename ReconstructType, index_type N>
struct state_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using bp = typename base::bp;
  using vector_type = typename base::vector_type;

  using domain_type = state_domain_t<N>;
  using reconstruct_type = ReconstructType;

  vector_type l, r;
  state_t() : l(base::make_vector()), r(base::make_vector()) {}

  constexpr inline void build_state(const reconstruct_type& reco,
                                    const domain_type& domain) {
    for (index_type i = bp::i0 + 2; i < bp::iN - 2; ++i) {
      l[i + 1] = reco.rc_left(i, domain.Cdtdx[i], domain.fCdtdx[i]);
      r[i] = reco.rc_right(i, domain.Cdtdx[i], domain.fCdtdx[i]);
    }
  }
};

template <typename StateType, index_type N, index_type NSt>
struct state_set_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using vector_type = typename base::vector_type;

  using state_type = StateType;
  using domain_type = typename state_type::domain_type;
  using reconstruct_type = typename state_type::reconstruct_type;

  domain_type domain;
  std::array<state_type, NSt> states;
  std::array<FLUID_VI, NSt> fluid_ids;

  state_set_t() = default;

  state_set_t(std::array<FLUID_VI, NSt> rcis) : fluid_ids(std::move(rcis)) {}

  constexpr inline void prepare(const hydro::hydro_physics_t& hydro,
                                const vector_type& cdx) {
    domain.build(hydro, cdx);
  }

  template <size_t Sp>
  constexpr inline void build(const std::array<reconstruct_type, Sp>& recos) {
    for (index_type i = 0; i < NSt; ++i) {
      states[i].build_state(recos[fluid_ids[i]], domain);
    }
  }

  // Build from individual parabola instances (no array indexing)
  template <typename ParabolaSet>
  constexpr inline void build(const ParabolaSet& parabolas) {
    for (index_type i = 0; i < NSt; ++i) {
      const auto& reco = parabolas.get(fluid_ids[i]);
      states[i].build_state(reco, domain);
    }
  }

  constexpr inline const state_type& operator[](const FLUID_VI fi) const {
    for (index_type i = 0; i < NSt; ++i) {
      if (fi == fluid_ids[i])
        return states[i];
    }
  }
};

  }  // namespace ppm

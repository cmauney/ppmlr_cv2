#pragma once

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

// state_set_t with individual state members (no FLUID_VI enum, no std::array)
template <typename StateType, index_type N>
struct state_set_t : ppmlr_base_t<N> {
  using base = ppmlr_base_t<N>;
  using vector_type = typename base::vector_type;

  using state_type = StateType;
  using domain_type = typename state_type::domain_type;
  using reconstruct_type = typename state_type::reconstruct_type;

  domain_type domain;

  // Individual state members for each variable (RHO, PRS, XVL)
  state_type rho_state;
  state_type prs_state;
  state_type xvl_state;

  state_set_t() = default;

  constexpr inline void prepare(const hydro::hydro_physics_t& hydro,
                                const vector_type& cdx) {
    domain.build(hydro, cdx);
  }

  // Build from individual parabola instances (direct access, no enum)
  template <typename ParabolaSet>
  constexpr inline void build(const ParabolaSet& parabolas) {
    rho_state.build_state(parabolas.para_rho, domain);
    prs_state.build_state(parabolas.para_prs, domain);
    xvl_state.build_state(parabolas.para_xvl, domain);
  }

  // Direct accessors (replaces operator[] with FLUID_VI)
  constexpr inline const state_type& rho() const { return rho_state; }
  constexpr inline const state_type& prs() const { return prs_state; }
  constexpr inline const state_type& xvl() const { return xvl_state; }

  constexpr inline state_type& rho() { return rho_state; }
  constexpr inline state_type& prs() { return prs_state; }
  constexpr inline state_type& xvl() { return xvl_state; }
};

  }  // namespace ppm

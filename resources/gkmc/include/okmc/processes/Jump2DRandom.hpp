// Jump2DRandom.hpp

#ifndef GKMC_Jump2DRandom_H
#define GKMC_Jump2DRandom_H

#include <helper/invars.hpp>

namespace gkmc {

template <class Sys, class T> class Jump2DRandom {
public:
  using System = Sys;
  using Species = typename Sys::Species;
  using ParticleType = T;
  using rate_type = typename System::value_type;

private:
  rate_type _rate;
  rate_type _dist;
  // RandReal<> _rnd{-invars::pi, invars::pi};
  RandReal<> _rnd{0.0, 1.0};
  Species _species;

public:
  Jump2DRandom(System &sys, const Species &sp, const rate_type &Em,
               const rate_type &v, const rate_type &dist)
      : _dist{dist}, _species{sp} {
    _rate = v * exp(-Em / (sys.props().temperature() * invars::kB));
  }

  const auto &species() const { return _species; }

  auto executeSelected(const size_t &pIndex, System &sys,
                       const rate_type subDice) {
    auto pid = sys.cookPid(_species, pIndex);
    const auto &particle = sys.template particle<T>(pid);

    auto th = 2 * invars::pi * subDice;
    typename System::Coords coords{
        {_dist * std::cos(th), _dist * std::sin(th)}};

    std::transform(std::begin(coords), std::end(coords),
                   std::begin(particle.coords()), std::begin(coords),
                   std::plus<rate_type>());
    sys.template changeCoords<T>(pid, coords);
  }

  auto executeSelected(const size_t &pIndex, System &sys) {
    executeSelected(pIndex, sys, _rnd());
  }

  const rate_type &ratePerParticle(const System &s) const { return _rate; }
};
}

#endif // !GKMC_Jump2DRandomPROCESS_H

/*
    auto th = _rnd();
    auto z = _dist * subDice;
    auto temp = std::sqrt(_dist*_dist - z*z);
    auto x = temp * std::cos(th);
    auto y = temp * std::sin(th);
    typename System::Coords coords{{x, y, z}};
*/

/*
                                    */
// typename System::Coords coords{{1.0, 1.0, 1.0}};

// Jump3DRandom.hpp

#ifndef GKMC_JUMP3DRandom_H
#define GKMC_JUMP3DRandom_H

#include <helper/RandReal.hpp>
#include <helper/invars.hpp>

namespace gkmc {

template <class Sys, class T> class Jump3DRandom {
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
  Jump3DRandom(const Species &sp, const rate_type &Em, const rate_type &v,
               const rate_type &dist, rate_type temperature)
      : _dist{dist}, _species{sp} {
    _rate = v * exp(-Em / (temperature * invars::kB));
  }

  const auto &species() const { return _species; }

  auto executeSelected(const size_t &pIndex, System &sys,
                       const rate_type subDice) {
    auto pid = sys.cookPid(_species, pIndex);
    const auto& particle = sys.template particle<T>(pid);
    auto coords = rndOnSphere(subDice, _rnd(), _dist);
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

#endif // !GKMC_JUMP3DRandomPROCESS_H

// Remove.hpp

#ifndef GKMC_Remove_H
#define GKMC_Remove_H

#include <helper/invars.hpp>

namespace gkmc {

template <class Sys, class T> class Remove {
public:
  using System = Sys;
  using Species = typename Sys::Species;
  using ParticleType = T;
  using rate_type = typename System::value_type;

private:
  rate_type _rate;
  Species _species;

public:
  Remove(rate_type temperature, const Species &sp, const rate_type &Em,
               const rate_type &v)
      : _species{sp} {
    _rate = v * exp(-Em / (temperature * invars::kB));
  }

  const auto &species() const { return _species; }

  auto executeSelected(const size_t &pIndex, System &sys,
                       const rate_type subDice = 0.0) {
    auto pid = sys.cookPid(_species, pIndex);
    sys.template remove<T>(pid);
  }

  const rate_type &ratePerParticle(const System &s) const { return _rate; }
};
}

#endif // !GKMC_Remove_H

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

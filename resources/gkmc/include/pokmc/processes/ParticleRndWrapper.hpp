// ParticleRndWrapper.hpp

#ifndef GKMC_ParticleRndWrapper_H
#define GKMC_ParticleRndWrapper_H

#include <helper/RandReal.hpp>

namespace gkmc {

template <class Proc, bool isRndRequired> class ParticleRndWrapper {
public:
  using System = typename Proc::System;
  using rate_type = typename System::value_type;
  ParticleRndWrapper(Proc process) : _process{process} {}
  void executeSelected(size_t pIndex, System &sys) {
    _process.executeSelected(pIndex, sys, 0.0);
  }
  decltype(auto) ratePerParticle(System &s) {
    return _process.ratePerParticle(s);
  }
  const auto &proc() const { return _process; }
private:
  Proc _process;
};

template <class Proc> class ParticleRndWrapper<Proc, true> {
public:
  using System = typename Proc::System;
  using rate_type = typename System::value_type;
  ParticleRndWrapper(Proc process) : _process{process} {}
  void executeSelected(size_t pIndex, System &sys) {
    _process.executeSelected(pIndex, sys, _rnd());
  }
  decltype(auto) ratePerParticle(System &s) {
    return _process.ratePerParticle(s);
  }
  const auto &proc() const { return _process; }
private:
  Proc _process;
  RandReal<> _rnd{0.0, 1.0};
};

}

#endif

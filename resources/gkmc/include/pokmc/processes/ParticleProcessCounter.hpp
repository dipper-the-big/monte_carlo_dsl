// ParticleProcessCounter.hpp

#ifndef GKMC_ParticleProcessCounter_H
#define GKMC_ParticleProcessCounter_H

/**
 * Process counter wrapper to know how many times a process
 * has been executed. Can be used for debugging and logging.
 * */
namespace gkmc {
template <class Proc> class ParticleProcessCounter {
public:
  using System = typename Proc::System;
  using rate_type = typename System::value_type;
  using ParticleType = typename Proc::ParticleType;

  ParticleProcessCounter(Proc process) : _process{process} {}
  void executeSelected(size_t pIndex, System &sys, rate_type dice) {
    ++_count;
    _process.executeSelected(pIndex, sys, dice);
  }
  void executeSelected(size_t pIndex, System &sys) {
    ++_count;
    _process.executeSelected(pIndex, sys);
  }
  decltype(auto) ratePerParticle(System &s) {
    return _process.ratePerParticle(s);
  }
  const auto &species() const { return _process.species(); }
  const auto &count() const { return _count; }
  const auto &proc() const { return _process; }
private:
  Proc _process;
  size_t _count = 0;
};
}

#endif

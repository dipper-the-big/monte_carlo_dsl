// VirtualParticleProcess.hpp

#ifndef GKMC_VirtualParticleProcess_H
#define GKMC_VirtualParticleProcess_H

#include <algorithm>
#include <array>
#include <assert.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include <helper/meta.hpp>

namespace gkmc {

template <class Sys, class Particle>
class IParticleProcess {
public:
  using System = Sys;
  using Species = typename Sys::Species;
  using ParticleType = Particle;
  using rate_type = typename System::value_type;

  virtual void executeSelected(size_t i, System &sys, rate_type dice) = 0;
  virtual void executeSelected(size_t i, System &sys) = 0;
  virtual const rate_type ratePerParticle(System &sys) = 0;
  virtual const Species &species(System &sys) = 0;
};

template <class P>
class VirtualizeParticleProc : public IParticleProcess<typename P::System, typename P::ParticleType> {
public:
  using System = typename P::System;
  using Species = typename System::Species;
  using ParticleType = typename P::ParticleType;
  using rate_type = typename System::value_type;
  VirtualizeParticleProc(P proc) : _proc{proc} {}
  void executeSelected(size_t i, System &sys, rate_type dice) override {
    _proc.executeSelected(i, sys, dice);  
  }
  void executeSelected(size_t i, System &sys) override {
    _proc.executeSelected(i, sys);  
  }
  const rate_type ratePerParticle(System &sys) override { return _proc.ratePerParticle(sys); }
  const Species &species(System &sys) override { return _proc.species(); }
private:
  P _proc;
};

} // namespace gkmc

#endif
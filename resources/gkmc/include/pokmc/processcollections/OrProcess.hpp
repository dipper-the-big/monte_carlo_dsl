// OrProcess.hpp

#ifndef GKMC_OrProcess_H
#define GKMC_OrProcess_H

#include <algorithm>
#include <memory>
#include <vector>
#include <random>

namespace gkmc {

template <class Proc>
class OrProcess {
public:
  using System = typename Proc::System;
  using ParticleType = typename Proc::ParticleType;

  void add(Proc proc) { _procs.emplace_back(std::forward<Proc>(proc)); }

  size_t size() const { return _procs.size(); }

  const auto& procs() const { return _procs; }

  auto isDecide(bool f) { _isDecide = f; }
  auto isDecide() const { return _isDecide; }

  auto orRate(System& sys) {
    auto res = 0.0;
    for (auto &it : this->_procs) {
      res += it.ratePerParticle(sys);
    }
    return res;
  }

  auto executeFor(System& sys, double deltaT) {
    if (this->_procs.empty()) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    auto nSystemTuple = totalSystemTuple(sys);
    for (auto &proc : this->_procs) {
      auto d = std::poisson_distribution<> {deltaT * proc.ratePerParticle(sys)};
      for (size_t i = 0; i < nSystemTuple; ++i) {
        auto times = d(gen);
        while (times--) proc.executeSelected(i, sys);
      }
    }
  }
  
  size_t totalSystemTuple(System &sys) {
    if (this->_procs.empty()) return 0;
    return sys.template particles<ParticleType>(this->_procs[0].species()).size();
  }

protected:
  bool _isDecide = true;
  std::vector<Proc> _procs;
  //std::vector<Species> _species;
};

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
  virtual const Species &species() = 0;
};

template <class P>
class VirtualizeParticleProc : public IParticleProcess<typename P::System, typename P::ParticleType> {
public:
  using System = typename P::System;
  using Species = typename System::Species;
  using ParticleType = typename P::ParticleType;
  using rate_type = typename System::value_type;
  VirtualizeParticleProc(P proc) : _proc{std::forward<P>(proc)} {}
  void executeSelected(size_t i, System &sys, rate_type dice) override {
    _proc.executeSelected(i, sys, dice);  
  }
  void executeSelected(size_t i, System &sys) override {
    _proc.executeSelected(i, sys);  
  }
  const rate_type ratePerParticle(System &sys) override { return _proc.ratePerParticle(sys); }
  const Species &species() override { return _proc.species(); }
private:
  P _proc;
};



template <class Sys, class Particle>
class OrProcessVirtual {
public:
  using System = Sys;
  using ParticleType = Particle;

  template <class Proc>
  void add(Proc p) { 
    _procs.emplace_back(std::make_unique<VirtualizeParticleProc<Proc>>(std::forward<Proc>(p)));
  }

  size_t size() const { return _procs.size(); }

  const auto& procs() const { return _procs; }

  auto isDecide(bool f) { _isDecide = f; }
  auto isDecide() const { return _isDecide; }

  auto orRate(System& sys) {
    auto res = 0.0;
    for (auto &it : this->_procs) {
      res += it->ratePerParticle(sys);
    }
    return res;
  }

  auto executeFor(System& sys, double deltaT) {
    if (this->_procs.empty()) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    auto nSystemTuple = totalSystemTuple(sys);
    for (auto &proc : this->_procs) {
      auto d = std::poisson_distribution<> {deltaT * proc->ratePerParticle(sys)};
      for (size_t i = 0; i < nSystemTuple; ++i) {
        auto times = d(gen);
        while (times--) proc->executeSelected(i, sys);
      }
    }
  }
  
  size_t totalSystemTuple(System &sys) {
    if (this->_procs.empty()) return 0;
    return sys.template particles<ParticleType>(this->_procs[0]->species()).size();
  }

protected:
  bool _isDecide = true;
  std::vector<std::unique_ptr<IParticleProcess<Sys, Particle>>> _procs;
};


} // !namespace gkmc

#endif //

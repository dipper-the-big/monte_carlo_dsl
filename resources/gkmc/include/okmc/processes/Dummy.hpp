// DummyProcess.hpp

#ifndef GKMC_DummyProcess_H
#define GKMC_DummyProcess_H

#include <map>

namespace gkmc {

template <class Sys, class T> class DummyParticleProcess {
  using Species = typename Sys::Species;
  double _rate;
  Species _sp;

public:
  using System = Sys;
  using ParticleType = T;
  DummyParticleProcess(double rate, const Species &sp) : _rate{rate}, _sp{sp} {}

  auto executeSelected(size_t particle, Sys &sys, double subDice) {
    return std::make_pair(particle, subDice);
  }

  auto execute(Sys &sys, double dice) { return rate(sys); }

  double rate(Sys &sys) {
    return _rate * sys.template particles<T>(_sp).size();
  }

  double ratePerParticle(Sys &sys) { return _rate; }

  const auto &species() const { return _sp; }
};

class DummySysProcess {
  double _rate;

public:
  DummySysProcess(double rate) : _rate{rate} {}

  template <class Sys> auto execute(Sys &sys, double dice) { return rate(sys); }

  template <class Sys> double rate(const Sys &s) const { return _rate; }
};

template <class Sys> class DummyBaseProcess {
public:
  virtual double execute(Sys &sys, double dice) = 0;
  virtual double rate(Sys &sys) = 0;
};

template <class Sys>
class DummyConcreteSysProcess : public DummyBaseProcess<Sys> {
public:
  DummyConcreteSysProcess(double rate) : _rate{rate} {}
  virtual double execute(Sys &sys, double dice) override { return _rate; };
  virtual double rate(Sys &sys) override { return _rate; };

private:
  double _rate;
};

template <class Sys, class T>
class DummyConcreteParticleProcess1 : public DummyBaseProcess<Sys> {
using Species = typename Sys::Species;
public:
  DummyConcreteParticleProcess1(double rate, Species &sp)
      : _rate{rate}, _sp{sp} {}

  virtual double execute(Sys &sys, double dice) override { return rate(sys); }

  virtual double rate(Sys &sys) override {
    return _rate * sys.template particles<T>(_sp).size();
  }

private:
  double _rate;
  Species _sp;
};

template <class Sys, class T>
class DummyConcreteParticleProcess2 : public DummyBaseProcess<Sys> {
using Species = typename Sys::Species;
public:
  DummyConcreteParticleProcess2(double rate, Species &sp)
      : _rate{rate}, _sp{sp} {}

  virtual double execute(Sys &sys, double dice) override { return rate(sys); }

  virtual double rate(Sys &sys) override {
    return _rate * sys.template particles<T>(_sp).size();
  }

private:
  double _rate;
  Species _sp;
};
}

#endif // !GKMC_JUMP3DRandomPROCESS_H

// SpawnSingle.hpp

#ifndef GKMC_SpawnSingle_H
#define GKMC_SpawnSingle_H

#include <helper/RandReal.hpp>

namespace gkmc {

/**
 * Process to add one particle at a time with a given rate
 * */
template <class T, class Species>
class SpawnSingle {
private:
  double _rate;
  Species _sp;
  RandReal<> _rnd {0., 1.};

public:
  SpawnSingle(double rate, Species sp) : _rate {rate}, _sp{sp} {}

  template <class System>
  void execute(System &sys, double dice = 0.0) {
    // TODO: use dice for rnd
    typename System::Coords newCoords;
    for (size_t i = 0; i < newCoords.size(); ++i) {
      newCoords[i] = _rnd() * sys.props().boxDimensions()[i];
    }
    sys.template add<T>(_sp, T{newCoords});
  }

  template <class System>
  const double& rate(const System &s) const { return _rate; }
};

template <class Sys, class T>
class SpawnSingleOfEach {
public:
  using System = Sys;
  using Species = typename Sys::Species;
private:
  double _rate;
  std::vector<Species> _sp;
  RandReal<> _rnd {0., 1.};

public:
  SpawnSingleOfEach(double rate, std::vector<Species> sp) : _rate {rate}, _sp{sp} {}

  void execute(System &sys, double dice = 0.0) {
    // TODO: use dice for rnd
    typename System::Coords newCoords;
    for (auto& it : _sp) {
      for (size_t i = 0; i < newCoords.size(); ++i) {
        newCoords[i] = _rnd() * sys.props().boxDimensions()[i];
      }
      sys.template add<T>(it, T{newCoords});
    }
  }

  const double& rate(const System &s) const { return _rate; }
};

}

#endif // !GKMC_SpawnSingleProcess_H

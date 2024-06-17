// Emission.hpp

#ifndef GKMC_Emission_H
#define GKMC_Emission_H

#include <iostream>

#include <helper/RandReal.hpp>

namespace gkmc {

template <class Sys, class P1, class P2, class SizeCalc> class Emission {
public:
  using System = Sys;
  using Species = typename Sys::Species;
  using ParticleType = P1;
  using rate_type = typename System::value_type;

private:
  SizeCalc _sizeCalc;
  rate_type _dist;
  rate_type _rate;
  RandReal<> _rnd{0.0, 1.0};
  Species _spOn;
  Species _spEmit;
  Species _spShift;
  int _sizeShift;

public:
  Emission(SizeCalc sizeCalc, rate_type dist, rate_type rate, Species spOn, Species spEmit,  Species spShift, int sizeShift = 1)
      : _sizeCalc{sizeCalc}, _dist{dist}, _rate{rate}, _spOn{spOn}, _spEmit{spEmit}, _spShift{spShift}, _sizeShift{sizeShift}  {}

  auto executeSelected(const size_t &pIndex, System &sys,
                       const rate_type subDice) {
    auto pid = sys.cookPid(_spOn, pIndex);
    auto& particle = sys.template particle<P1>(pid);
    /*
    std::cout << "\n ================== emitting ================== \n";
    std::cout << particle.size() << " size\n";
    */
    auto away = _sizeCalc(sys, particle, pid) + _dist;
    //  http://corysimon.github.io/articles/uniformdistn-on-sphere/
    auto th = 2 * invars::pi * _rnd(); // subDice; TODO
    auto phi = acos(1. - 2. * _rnd());
    typename System::Coords coords{{away * std::sin(phi) * std::cos(th),
                                    away * std::sin(phi) * std::sin(th),
                                    away * std::cos(phi)}};

    std::transform(std::begin(coords), std::end(coords),
                   std::begin(particle.coords()), std::begin(coords),
                   std::plus<rate_type>());
    if (particle.size() - _sizeShift <= _sizeShift) {
      sys.template remove<P1>(pid);
      sys.template add<P2>(_spShift, P2{particle.coords(), particle.size() - _sizeShift});
    } else {
      sys.template changeSize <P1>(pid, particle.size() - _sizeShift);
    }
    sys.template add<P2>(_spEmit, P2{coords, _sizeShift});
  }

  const rate_type &ratePerParticle(const System &s) const { return _rate; }

  const Species &species() const { return _spOn; }
};

}

#endif // !GKMC_Emission_H

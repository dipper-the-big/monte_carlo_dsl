// Jump3DBcc.hpp

#ifndef GKMC_JUMP3DBCC_H
#define GKMC_JUMP3DBCC_H

#include <assert.h>
#include <helper/RandReal.hpp>
#include <helper/invars.hpp>

namespace gkmc {

template <class Sys, class T>
class Jump3DBcc {
public:
  using System = Sys;
  using Species = typename System::Species;
  using ParticleType = T;
  using Coords = typename System::Coords;
  using rate_type = typename System::value_type;

private:
  rate_type _rate;
  std::array<rate_type, 8> _prob;
  Coords _jump;
  std::array<Coords, 8> _signs;
  Species _species;
  //RandReal<> _rnd{0., 1.};

public:
  Jump3DBcc(System& sys, const Species &sp, const rate_type &Em, const rate_type &v) : _species{sp} {
    auto jumpFrac = _prob[0] = 1. / 8;
    assert(sys.props().latticeConstant() > 0.00006);
    _jump[0] = _jump[1] = _jump[2] = sys.props().latticeConstant() * 0.5;
    _signs[0] = Coords{{-1, -1, -1}};
    for (auto i = 1; i < 8; ++i) {
      _prob[i] = _prob[i - 1] + jumpFrac;
      _signs[i] = Coords{{ (2 * (i / 4)) - 1.,
                                   (2 * ((i / 2) % 2)) - 1.,
                                   (2 * (i % 2)) - 1.}};
    }
    _rate = v * exp(-Em / (sys.props().temperature() * invars::kB));
  }

  const auto &species() const { return _species; }

  auto executeSelected(size_t pIndex, System &sys, rate_type subDice) {
    auto urn = subDice;//_rnd();
    auto pid = sys.cookPid(_species, pIndex);
    auto particle = sys.template particle<T>(pid);
    auto it = std::lower_bound(std::begin(_prob), std::end(_prob), urn);
    auto index = std::distance(std::begin(_prob), it);
    const auto &sign = _signs[index];
    //std::cerr<<index<<std::endl;
    auto coords = _jump;
    std::transform(std::begin(coords), std::end(coords), std::begin(sign),
                   std::begin(coords), std::multiplies<rate_type>());
    //std::cout << coords[0] << '\t' << coords[1] << '\t' << coords[2] << '\n';
    std::transform(std::begin(coords), std::end(coords),
                   std::begin(particle.coords()), std::begin(coords),
                   std::plus<rate_type>());
    //std::cout << coords[0] << '\t' << coords[1] << '\t' << coords[2] << '\n';
    sys.template changeCoords<T>(pid, coords);
  }

  const rate_type& ratePerParticle(const System &s) const { return _rate; }
};

template <class Sys, class T>
class Jump1DBcc {
public:
  using System = Sys;
  using Species = typename System::Species;
  using ParticleType = T;
  using Coords = typename System::Coords;
  using rate_type = typename System::value_type;


private:
  static constexpr double kB = 8.6173324e-5;
  double _rate;
  std::array<double, 8> _prob;
  Coords _jump;
  Species _species;
  RandReal<> _rnd{0., 1.};

auto _areSame(const std::array<double, 3>& prev, const std::array<double, 3>& next) {
  for (auto i : {0, 1, 2}) {
    if (fabs(next[i] - prev[i]) > 0.00001) return false;
  }
  return true;
}

void _jumpRnd(System &sys, T& p, typename System::Pid &pid, double urn) {
   auto sign = Coords{};
   for (auto i : {0, 1, 2}) {
     auto urn = _rnd();
     if (urn < 0.5) sign[i] = -1;
     else sign[i] = 1;
   }
   auto coords = _jump;
   std::transform(std::begin(coords), std::end(coords), std::begin(sign),
                  std::begin(coords), std::multiplies<double>());
   std::transform(std::begin(coords), std::end(coords),
                  std::begin(p.coords()), std::begin(coords),
                  std::plus<double>());
   sys.template changeCoords<T>(pid, coords);
}
public:
  Jump1DBcc(System &sys, const Species &sp, double Em, double v)
      : _species{sp} {
    _jump[0] = _jump[1] = _jump[2] = sys.props().latticeConstant() * 0.5;
    _rate = v * exp(-Em / (sys.props().temperature() * kB));
  }

  const auto &species() const { return _species; }

  auto executeSelected(size_t pIndex, System &sys, double subDice) {
    auto pid = sys.cookPid(_species, pIndex);
    auto particle = sys.template particle<T>(pid);
    if (_areSame(particle.coords(), particle.prevCoords())) {
      _jumpRnd(sys, particle, pid, subDice);
      return;
    }
    auto urn = _rnd();//subDice;
    if (urn < 0.5) {
      sys.template changeCoords<T>(pid, particle.prevCoords());
    } else {
      Coords newCoords{};
      auto temp = particle.coords();
      for(auto i : {0, 1, 2}) {
        newCoords[i] = temp[i] + (temp[i] - particle.prevCoords()[i]);
      }
      sys.template changeCoords<T>(pid, newCoords);
    }
  }

  const rate_type& ratePerParticle(const System &s) const { return _rate; }
};

}

#endif // !GKMC_JUMP3DBCCPROCESS_H

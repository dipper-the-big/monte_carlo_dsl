// diffuse.hpp

#ifndef GKMC_Diffuse_H
#define GKMC_Diffuse_H

#include <helper/RandReal.hpp>
#include <helper/invars.hpp>

namespace gkmc {

/**
 * Simple 3D random diffusion for MCD algo.
 */
template <class Particle, class Species>
class Diffuse3D {
  Species _sp;
  double _dist;
  RandReal<> _rnd{0.0, 1.0};
public:
  Diffuse3D(Species sp, double dist) : _sp{sp}, _dist{dist} {}
  template <class System>
  auto operator () (System &sys, double delT) {
    auto& particles = sys.template particles<Particle>(_sp);
    size_t i = 0;
    for (auto& particle : particles) {
      auto coords = rndOnSphere(_rnd(), _rnd(), _dist);
      std::transform(std::begin(coords), std::end(coords),
                     std::begin(particle.coords()), std::begin(coords),
                     std::plus<double>());
      sys.template changeCoords<Particle>(sys.cookPid(_sp, i++), coords);
    }
  }
};

}

#endif // !GKMC_Diffuse_H

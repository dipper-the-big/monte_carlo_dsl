// Absorb.hpp

#ifndef GKMC_Absorb_H
#define GKMC_Absorb_H

#include <cmath>
#include <helper/dist.hpp>

namespace gkmc {

/**
 * Binary reaction that creates a particle with the given species
 * and size as the addition of the two. The original particles are
 * merged and hence removed.
 * */
template <class Species>
struct AbsorbEq {
  AbsorbEq(Species sp) : _sp{sp} {}
  template <class Sys, class Particle>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, Particle &p, Particle &n) {
    sys.template remove<Particle>(pid);
    sys.template remove<Particle>(nid);
    sys.template add<Particle>(_sp, Particle{p.coords(), p.size() + n.size()});
    return true;
  }
private:
  Species _sp;
};

/**
 * Binary reaction that creates a particle with the species that
 * has bigger size. The size of the resulting particle is same 
 * as the addition of the two. 
 * */
class Absorb {
public:
  template <class Sys, class P, class N>
  auto operator()(Sys &ps, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, N &n) {
    if (p.size() > n.size()) {
      ps.template changeSize<P>(pid, p.size() + n.size());
      ps.template remove<N>(nid);
    } else {
      ps.template changeSize<N>(nid, n.size() + p.size());
      ps.template remove<P>(pid);
    }
    return true;
  }
};

/**
 * Binary reaction that creates a particle with the species that
 * has strictly bigger size. If the sizes are equal then the 
 * new particle has species that is specified in ctor.
 * */
template <class P2, class Sp>
class AbsorbNew {
  Sp _sp;
  const int _shiftSize;
public:
  AbsorbNew(const Sp& sp, int shiftSize = 2) : _sp{sp}, _shiftSize{shiftSize} {}
  template <class Sys, class P, class N>
  auto operator()(Sys &ps, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, N &n) {
    if (p.size() + n.size() == _shiftSize) {
      ps.template remove<P>(pid);
      ps.template remove<N>(nid);
      ps.template add<P2>(_sp, P2{p.coords(), _shiftSize});
    } else if (p.size() > n.size()) {
      ps.template remove<N>(nid);
      ps.template changeSize<P>(pid, p.size() + n.size());
    } else {
      ps.template remove<P>(pid);
      ps.template changeSize<N>(nid, n.size() + p.size());
    }
    return true;
  }
};

}
#endif

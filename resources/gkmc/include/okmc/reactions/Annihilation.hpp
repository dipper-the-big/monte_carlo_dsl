// Annihilation.hpp

#ifndef GKMC_Annihilation_H
#define GKMC_Annihilation_H

#include <cmath>
#include <helper/dist.hpp>

namespace gkmc {

/**
 * Binary reaction that removes both the particles
 * */
struct AnnihilationEq {
  template <class Sys, class Particle>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, Particle &p, Particle &n) {
    sys.template remove<Particle>(pid);
    sys.template remove<Particle>(nid);
    return true;
  }
};

/**
 * Binary reaction that removes the particle with smaller size
 * and reduces the size of the bigger one.
 * */
class Annihilation {
public:
  template <class Sys, class P, class N>
  auto operator()(Sys &SystemTuple, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, N &n) {
    if (p.size() > n.size()) {
      p.size(p.size() - n.size());
      SystemTuple.template edited<P>(pid);
      SystemTuple.template remove<N>(nid);
    } else if (p.size() < n.size()) {
      n.size(n.size() - p.size());
      SystemTuple.template edited<N>(nid);
      SystemTuple.template remove<P>(pid);
    } else {
      SystemTuple.template remove<P>(pid);
      SystemTuple.template remove<N>(nid);
    }
    return true;
  }
};

/**
 * Binary reaction that calls a function which might add particle
 * or edit system and return true in which case the reaction removes
 * both the particles. If the function returns false then it acts like
 * the above Annihilation reaction.
 * */
template <class Fn>
class AnnihilationFn {
  Fn _fn;
public:
  AnnihilationFn(Fn fn) : _fn{fn} {}
  template <class Sys, class P, class N>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, N &n) {
    auto hasAdded = _fn(sys, pid, nid, p, n);
    if (hasAdded || p.size() == n.size()) {
      sys.template remove<P>(pid);
      sys.template remove<N>(nid);
    } else if (p.size() > n.size()) {
      sys.template changeSize<P>(pid, p.size() - n.size());
      sys.template remove<N>(nid);
    } else { // if (p.size() < n.size())
      sys.template changeSize<N>(nid, n.size() - p.size());
      sys.template remove<P>(pid);
    }
    return true;
  }
};

/**
 * Binary reaction that removes the particle with smaller size
 * and introduces a particle with the given species if the reduced
 * size of the particle is same as the given size in ctor. If both the sizes
 * are equal then it removes both the particles (given that the size given
 * in ctor is not zero.)
 * */
template <class P2, class Sp>
class AnnihilationNew {
  Sp _sp;
  const int _shiftSize = 1;
public:
  AnnihilationNew(Sp sp, int shiftSize = 1) : _sp{sp}, _shiftSize{1} {}
  template <class Sys, class P, class N>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, N &n) {
    auto newSize = abs(p.size() - n.size());
    if (newSize == _shiftSize) {
      auto& coords = (p.size() > n.size()) ? p.coords() : n.coords();
      sys.template add<P2>(_sp, P2{coords, newSize});
      sys.template remove<N>(nid);
      sys.template remove<P>(pid);
    } else if (p.size() > n.size()) {
      p.size(p.size() - n.size());
      sys.template edited<P>(pid);
      sys.template remove<N>(nid);
    } else if (p.size() < n.size()) {
      n.size(n.size() - p.size());
      sys.template edited<N>(nid);
      sys.template remove<P>(pid);
    } else {
      sys.template remove<P>(pid);
      sys.template remove<N>(nid);
    }
    return true;
  }
};

}
#endif

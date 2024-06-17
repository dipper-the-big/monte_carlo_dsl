// Binary.hpp

#ifndef GKMC_Binary_H
#define GKMC_Binary_H

namespace gkmc {

/**
 * The input binary reaction only initiates if the type of 
 * particles are same and equal to the input particle type.
 * */
template <class Reaction, class Neigh, class Tol, class P>
class ReflexiveType {
public:
  ReflexiveType(Reaction r, Neigh neigh, Tol tol)
    : _reaction{r}, _neigh{neigh}, _tol{tol} {}
  template <class SystemTuple>
  auto operator()(SystemTuple &ps, const typename SystemTuple::Pid &pid, P &p) {
    auto neigh = _neigh.template neighbour<P, P>(ps, pid, _tol);
    if (neigh.first == true) {
      return _reaction(ps, pid, neigh.second, p,
                       ps.template particle<P>(neigh.second));
    }
    return false;
  }

  auto &reaction() { return _reaction; }
  auto &neigh() { return _neigh; }
  auto &distFn() { return _tol; }

private:
  Reaction _reaction;
  Neigh _neigh;
  Tol _tol;
};

/**
 * The input binary reaction only initiates between the input
 * pair of types and species given.
 * The specialized version for the same particle type follows.
 * */
template <class Reaction, class Neigh, class Tol, class Sp, class P1, class P2>
class Symmetric {
public:
  using Species = Sp;
  template <class SpRange>
  Symmetric(Reaction r, SpRange &&sp1, SpRange &&sp2, Neigh neigh, Tol tol)
    : _reaction{r}, _sp1{begin(sp1), end(sp1)}, _sp2{begin(sp2), end(sp2)},
      _neigh{neigh}, _tol{tol} {}
  template <class SystemTuple>
  auto operator()(SystemTuple &ps, const typename SystemTuple::Pid &pid,
                  P1 &p) {
    if (_sp1.find(ps.template getSpecies<P1>(pid)) == std::end(_sp1))
      return false;
    for (auto &s : _sp2) {
      auto neigh = _neigh.template neighbour<P1, P2>(ps, pid, _tol, s);
      if (neigh.first == true) {
        return _reaction(ps, pid, neigh.second, p,
                         ps.template particle<P2>(neigh.second));
      }
    }
    return false;
  }

  template <class SystemTuple>
  auto operator()(SystemTuple &ps, const typename SystemTuple::Pid &pid,
                  P2 &p) {
    if (_sp2.find(ps.template getSpecies<P2>(pid)) == std::end(_sp2))
      return false;
    for (auto &s : _sp1) {
      auto neigh = _neigh.template neighbour<P2, P1>(ps, pid, _tol, s);
      if (neigh.first == true) {
        return _reaction(ps, pid, neigh.second, p,
                         ps.template particle<P1>(neigh.second));
      }
    }
    return false;
  }

  template <class SystemTuple, class T>
  auto operator()(SystemTuple &ps, const typename SystemTuple::Pid &pid, T &p) {
    // TODO: Warning on unintended calls
    return false;
  }

  auto &reaction() { return _reaction; }
  auto &neigh() { return _neigh; }
  auto &distFn() { return _tol; }

private:
  Reaction _reaction;
  std::set<Species> _sp1;
  std::set<Species> _sp2;
  Neigh _neigh;
  Tol _tol;
};

/**
 * The specialized template used when the particle type is same.
 * */
template <class Reaction, class Neigh, class Tol, class Sp, class P>
class Symmetric<Reaction, Neigh, Tol, Sp, P, P> {
public:
  using Species = Sp;
  Symmetric(Reaction r, const Sp &sp, Neigh neigh, Tol tol)
      : _reaction{r}, _sp1{sp}, _sp2{sp}, _neigh{neigh}, _tol{tol} {}
  Symmetric(Reaction r, const Sp &sp1, const Sp &sp2, Neigh neigh, Tol tol)
      : _reaction{r}, _sp1{sp1}, _sp2{sp2}, _neigh{neigh}, _tol{tol} {}
  template <class SpRange>
  Symmetric(Reaction r, SpRange sp1, SpRange sp2, Neigh neigh, Tol tol)
      : _reaction{r}, _sp1{std::begin(sp1), std::end(sp1)}, _sp2{std::begin(sp2), std::end(sp2)},
        _neigh{neigh}, _tol{tol} {}
  template <class Particles>
  auto operator()(Particles &ps, const typename Particles::Pid &pid, P &p) {
    std::set<Species> sp;
    if (_sp1.find(ps.template getSpecies<P>(pid)) != std::end(_sp1)) sp = _sp2;
    else if (_sp2.find(ps.template getSpecies<P>(pid)) != std::end(_sp2)) sp = _sp1;
    else return false;
    for (auto& s : sp) {
      auto neigh = _neigh.template neighbour<P, P>(ps, pid, _tol, s);
      if (neigh.first == true) {
        return _reaction(ps, pid, neigh.second, p,
                         ps.template particle<P>(neigh.second));
      }
    }
    return false;
  }

  template <class Particles, class T>
  auto operator()(Particles &ps, const typename Particles::Pid &pid,
                  T &p) {
    // TODO: Warning on unintended calls
    return false;
  }

  auto& reaction() { return _reaction; }
  auto& neigh() { return _neigh; }
  auto& distFn() { return _tol; }
private:
  Reaction _reaction;
  std::set<Species> _sp1;
  std::set<Species> _sp2;
  Neigh _neigh;
  Tol _tol;
};

// alias
// TODO: make it not with all but only if species are equal
template <class Reaction, class Neigh, class Tol, class Sp, class P>
using Reflexive = Symmetric<Reaction, Neigh, Tol, Sp, P, P>;

/**
 * The reaction initiates between any of the input particle types and
 * species.
 * */
template <class Reaction, class Neigh, class Tol, class Sp, class... Ps>
class Any {
public:
  template <class SpRange>
  Any(Reaction r, SpRange &&sp, Neigh neigh, Tol tol)
    : _reaction{r}, _sp{begin(sp), end(sp)}, _neigh{neigh}, _tol{tol} {}

  template <class Sys, class P, class T, class... Ts>
  auto _helper(Sys &ps, const typename Sys::Pid &pid, P &p) {
    for (auto &s : _sp) {
      if (!ps.template isSpecies<T>(s)) continue;
      auto neigh = _neigh.template neighbour<P, T>(ps, pid, _tol, s);
      if (std::get<0>(neigh) == true) {
        return _reaction(ps, pid, std::get<1>(neigh), p,
                         ps.template particle<T>(neigh.second));
      }
    }
    if (sizeof...(Ts) == 0)
      return false;
    else
      return _helper<Sys, P, Ts...>(ps, pid, p);
  }

  template <class Sys, class P>
  auto _helper(Sys &ps, const typename Sys::Pid &pid, P &p) {
    return false;
  }

  template <class Sys, class P>
  auto operator()(Sys &ps, const typename Sys::Pid &pid, P &p) {
    if (_sp.find(ps.template getSpecies<P>(pid)) == std::end(_sp)) return false;
    return _helper<Sys, P, Ps...>(ps, pid, p);
  }

  auto &reaction() const { return _reaction; }
  auto &neigh() const { return _neigh; }
  auto &distFn() const { return _tol; }

private:
  Reaction _reaction;
  std::set<Sp> _sp;
  Neigh _neigh;
  Tol _tol;
};

/**
 * The reaction initiates for all the particles of the given types.
 * */
template <class Reaction, class Neigh, class Tol, class... Ps> class All {
public:
  All(Reaction r, Neigh neigh, Tol tol)
    : _reaction{r}, _neigh{neigh}, _tol{tol} {}

  template <class Sys, class P, class T, class... Ts>
  auto _helper(Sys &ps, const typename Sys::Pid &pid, P &p) {
    auto neigh = _neigh.template neighbour<P, T>(ps, pid, _tol);
    if (std::get<0>(neigh) == false) {
      if (sizeof...(Ts) == 0)
        return false;
      else
        return _helper<Sys, P, Ts...>(ps, pid, p);
    }
    return _reaction(ps, pid, std::get<1>(neigh), p,
                     ps.template particle<T>(neigh.second));
  }

  template <class Sys, class P>
  auto _helper(Sys &ps, const typename Sys::Pid &pid, P &p) {
    return false;
  }

  template <class Sys, class P>
  auto operator()(Sys &ps, const typename Sys::Pid &pid, P &p) {
    return _helper<Sys, P, Ps...>(ps, pid, p);
  }

  auto &reaction() const { return _reaction; }
  auto &neigh() const { return _neigh; }
  auto &distFn() const { return _tol; }

private:
  Reaction _reaction;
  Neigh _neigh;
  Tol _tol;
};
}
#endif

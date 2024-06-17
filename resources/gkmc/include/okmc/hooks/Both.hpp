// Both.hpp

#ifndef GKMC_Both_H
#define GKMC_Both_H

namespace gkmc {

/**
 * Combiner hook for two hooks, multiple hooks can also be combined
 **/
template <class A, class B> class Both {
public:
  Both(A a, B b) : _a{a}, _b{b} {}

  template <class Sys, class Particle>
  void added(Sys &sys, Particle &p, const typename Sys::Pid &pid) {
    _a.added(sys, p, pid);
    _b.added(sys, p, pid);
  }
  template <class Sys, class Particle>
  void edited(Sys &sys, Particle &p, const typename Sys::Pid &pid) {
    _a.edited(sys, p, pid);
    _b.edited(sys, p, pid);
  }
  template <class Sys, class Particle>
  void editedCoords(Sys &sys, Particle &p, const typename Sys::Pid &pid,
                    const typename Sys::Coords &c) {
    _a.editedCoords(sys, p, pid, c);
    _b.editedCoords(sys, p, pid, c);
  }
  template <class Sys, class Particle>
  void editedSize(Sys &sys, Particle &p, const typename Sys::Pid &pid,
                    const int &c) {
    _a.editedSize(sys, p, pid, c);
    _b.editedSize(sys, p, pid, c);
  }

  template <class Sys, class Particle>
  void removed(Sys &sys, Particle &p, const typename Sys::Pid &pid) {
    _a.removed(sys, p, pid);
    _b.removed(sys, p, pid);
  }
  template <class Sys, class Particle>
  void erased(Sys &sys, Particle &p, const typename Sys::Pid &pid) {
    _a.erased(sys, p, pid);
    _b.erased(sys, p, pid);
  }
  template <class Sys> void update(Sys &s, double t, double dt) {
    _a.update(s, t, dt);
    _b.update(s, t, dt);
  }

  auto& first() { return _a; }
  auto& second() { return _b; }
private:
  A _a;
  B _b;
};

} // namespace !gkmc

#endif // !GKMC_Both_H

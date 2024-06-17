// Eager.hpp

#ifndef GKMC_Eager_H
#define GKMC_Eager_H

#include <array>
#include <assert.h>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <helper/meta.hpp>

namespace gkmc {

/**
 * A hook that calls reactions as soon as a particle is edited
 * */
template <class SysProps, class Reaction> struct Eager {
  static constexpr int dims = SysProps::dims;
  using Coords = typename SysProps::Coords;
  using Pid = typename SysProps::Pid;
  using Species = typename SysProps::Species;

  Eager(Reaction r) : _reaction{r} {}

  template <class Sys, class T> void added(Sys &sys, T &p, const Pid &pid) {
    _execute(sys, p, pid);
  }

  template <class Sys, class T>
  void editedCoords(Sys &sys, T &p, const Pid &pid, const Coords &pPrior) {
    _execute(sys, p, pid);
  }

  template <class Sys, class T>
  void editedSize(Sys &sys, T &p, const Pid &pid, const int &pPrior) {
    _execute(sys, p, pid);
  }

  template <class Sys, class T>
  void edited(Sys &sys, T &p, const Pid &pid, T &pPrior) {
    _execute(sys, p, pid);
  }

  template <class Sys, class T>
  void edited(Sys &sys, T &p, const Pid &pid) {
    _execute(sys, p, pid);
  }

  template <class Sys, class T> void removed(Sys &, T &p, const Pid &pid) {}

  template <class Sys, class T> void erased(Sys &,T &p, const Pid &pid) {}

  template <class Sys> void update(Sys &sys, double, double) {
    sys.commit();
  }
private:
  template <class Sys, class T> auto _execute(Sys &sys, T&p, const Pid &pid) {
    auto res = _reaction(sys, pid, p);
    sys.commit();
    return res;
  }

private:
  Reaction _reaction;
};
} // namespace !gkmc

#endif // !GKMC_Eager_H
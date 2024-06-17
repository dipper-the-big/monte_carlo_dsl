#include <cmath>

#ifndef GKMC_PBC_H
#define GKMC_PBC_H

namespace gkmc {

/**
 * Periodic boundary conditions
 * */
struct PBC {
template <class Sys, class P>
bool operator()(Sys& sys, const typename Sys::Pid& pid, const P &p) {
  auto& coords = p.coords();
  auto& box = sys.props().boxDimensions();
  auto ret = false;
  for (size_t i = 0; i < coords.size(); ++i) {
    if (coords[i] > box[i] || coords[i] < 0.0) {
      auto rem = std::remainder(coords[i], box[i]);
      if (rem < 0.0)  rem += box[i];
      auto quo = int(coords[i] / box[i]);
      if (quo <= 0) --quo;
      sys.template changeCoordsAdd<P>(pid, rem, quo, i);
      ret = true;
    }
  }
  return ret;
}
};

/**
 * Absorbing boundary conditions
 * */
struct ABC {
template <class Sys, class P>
bool operator()(Sys& sys, const typename Sys::Pid& pid, const P &p) {
  auto& coords = p.coords();
  auto& box = sys.props().boxDimensions();
  for (size_t i = 0; i < coords.size(); ++i) {
    if (coords[i] > box[i] || coords[i] < 0.0) {
      sys.template remove<P>(pid);
      return true;
    }
  }
  return false;
}
};

}

#endif

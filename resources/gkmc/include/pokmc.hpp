// pokmc.hpp

#ifndef GKMC_POKMC_H
#define GKMC_POKMC_H

#include <algorithm>
#include <cmath>
#include <map>

namespace gkmc {
/*!
 * Based on the paper - A GPU-based parallel Object kinetic Monte Carlo
 * algorithm for the evolution of defects in irradiated materials
 * F. Jim ́enez, C.J. Ortiz
 * */
template <class ProcessCollection, class Predicate, class System>
// requires Predicate to return a boolean checkable when passed a sys.
// requires System to have update method defined.
auto pokmc(ProcessCollection &&procs,
         Predicate &&predicate = Predicate(), System &&sys = System(), double w = 1.0) {
  double rMax = procs.maxAndRate(sys);
  double delT = w / rMax;
  auto systemTime = 0.0;
  while(predicate(sys, systemTime)) {
    procs.executeFor(sys, delT);
    sys.update(systemTime, delT);
    systemTime += delT;
  }
  return systemTime;
}

}  // namespace gkmc

#endif // !GKMC_POKMC_H

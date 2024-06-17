// mcd.hpp

#ifndef GKMC_MCD_H
#define GKMC_MCD_H

#include <algorithm>
#include <cmath>
#include <map>

namespace gkmc {

/* 
* Monte Carlo Diffusive Reactive basic code with fixed delT
*/
template <class Diffusion, class Predicate, class System>
// requires callable diffusion Process
// requires Predicate that return a boolean checkable
// requires System to have update method defined.
auto mcd(Diffusion diffuse, Predicate &&predicate = Predicate(),
           System &&sys = System(), double delT = 1.0) {
  auto systemTime = 0.0;
  while (predicate(sys, systemTime)) {
    diffuse(sys, delT);
    systemTime += delT;
    sys.update(delT, systemTime);
  }
  return systemTime;
}

} // namespace gkmc

#endif // !GKMC_MCD_H

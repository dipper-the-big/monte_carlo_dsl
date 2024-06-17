// bkl.hpp

#ifndef GKMC_BKL_H
#define GKMC_BKL_H

#include <algorithm>
#include <cmath>
#include <map>

namespace gkmc {
enum class bklStatus : bool { NoProcesses, Finished };

/*!
 * Based on wikipedia algo for OKMC
 * link- https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC
 * */
template <class ProcessCollection, class Random, class Predicate, class System>
// requires ProcessCollection to have cumulativeRates and execute.
// requires Random object to return a double between 0, 1
// requires Predicate to return a boolean checkable when passed a sys.
// requires System to have update method defined.
auto bkl(ProcessCollection &&procs, Random &&rnd = Random(),
         Predicate &&predicate = Predicate(), System &&sys = System()) {
  using std::begin; using std::end; using std::make_pair;
  using std::distance; using std::lower_bound;
  constexpr double epsilon = 1e-6;
  // Step 1
  auto systemTime = 0.0;
  // Step 2, also with system ctor call
  while (predicate(sys, systemTime)) {
    // Step 3, 4. cumulative rates are a fn. of system and procs. If there are
    // two processes and rate of first process is 10, and that of 2nd is 5;
    // rates will be [10, 15]
    auto rates = procs.cumulativeRates(sys);
    if (rates.empty() || rates[rates.size() - 1] < epsilon) {
      return make_pair(systemTime, bklStatus::NoProcesses);
    }
    // Step 5, 6.
    auto dice = rnd() * rates[rates.size() - 1];
    auto picked = lower_bound(begin(rates), end(rates), dice);
    auto processIndex = distance(begin(rates), picked);
    // step 7
    procs.executeSelected(processIndex, sys, dice);
    double delT = std::log(1.0 / rnd()) / (rates[rates.size() - 1]);
    sys.update(systemTime, delT);
    // step 8, 9
    systemTime += delT;
  }
  return make_pair(systemTime, bklStatus::Finished);
}

}  // namespace gkmc

#endif // !GKMC_BKL_H

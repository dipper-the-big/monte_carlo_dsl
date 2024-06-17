// PrototypeProcessCollection.hpp

#ifndef GKMC_PrototypeProcessCollection_H
#define GKMC_PrototypeProcessCollection_H

#include <algorithm>
#include <vector>
#include <memory>
#include <okmc/processes/VirtualProcess.hpp>

namespace gkmc {

template <class System>
class PrototypeProcessCollection {
  using Species = typename System::Species;
  using rate_type = typename System::value_type;
public:
  auto executeSelected(size_t i, System &sys, rate_type dice) {
    auto prevCumulativeRate = 0.0;
    if (i != 0) prevCumulativeRate = _rateSum[i - 1];
    auto subDice = dice - prevCumulativeRate;
    return _procs[i]->execute(sys, subDice);
  }

  auto execute(System &sys, rate_type dice) {
    auto pickedCumulative =
        std::lower_bound(begin(_rateSum), end(_rateSum), dice);
    auto indexProcess = std::distance(begin(_rateSum), pickedCumulative);
    return executeSelected(indexProcess, sys, dice);
  }

  template <class Proc>
  void add(Proc p) {
    _procs.emplace_back(std::make_unique<Virtualize<Proc>>(std::forward<Proc>(p)));
  }

  size_t size() const { return _procs.size(); }

  rate_type rate(const System &sys) {
    cumulativeRate(sys);
    return _rateSum[size() - 1];
  }

  const std::vector<rate_type> &cumulativeRates(System &sys) {
    rate_type cumRate = 0.0;
    _rateSum.resize(_procs.size());
    for (size_t i = 0; i < _procs.size(); ++i) {
      cumRate += _procs[i]->rate(sys);
      _rateSum[i] = cumRate;
    }
    return _rateSum;
  }

  const auto& procs() const { return _procs; }

protected:
  std::vector<std::unique_ptr<IProcess<Sys>>> _procs;
  std::vector<rate_type> _rateSum;
};

} // !namespace gkmc

#endif //

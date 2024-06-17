// TupleProcessCollection.hpp

#ifndef GKMC_TupleProcessCollection_H
#define GKMC_TupleProcessCollection_H

#include <array>
#include <tuple>

namespace gkmc {

// TODO: add compile check for if all procs are different
template <class System, class... Procs> class TupleProcessCollection {
public:
  using rate_type = typename System::value_type;
  constexpr static std::size_t size = sizeof...(Procs);
  using seq = std::make_index_sequence<size>;
  TupleProcessCollection(Procs... procs) : _procs{procs...} {}
  auto executeSelected(size_t i, System &sys, rate_type dice) {
    auto prevCumulativeRate = 0.0;
    if (i != 0) prevCumulativeRate = _rateSum[i - 1];
    auto subDice = dice - prevCumulativeRate;
    return _executeHelper<size - 1>(i, sys, subDice);
  }
  auto execute(System &sys, rate_type dice) {
    auto pickedCumulative =
        std::lower_bound(begin(_rateSum), end(_rateSum), dice);
    auto indexProcess = std::distance(begin(_rateSum), pickedCumulative);
    return executeSelected(indexProcess, sys, dice);
  }

  const std::array<rate_type, size> &cumulativeRates(System &sys) {
    _cumulativeHelper(sys, seq{});
    return _rateSum;
  }
  rate_type rate(System &sys) {
    cumulativeRates(sys);
    if (size == 0) return 0.0;
    return _rateSum[size - 1];
  }
  const auto &procs() const { return _procs; }

protected:
  template <size_t index>
  auto _executeHelper(size_t cur, System &sys, rate_type dice,
                      typename std::enable_if<index != 0>::type *dummy = 0) {
    if (cur == index) return std::get<index>(_procs).execute(sys, dice);
    return _executeHelper<index - 1>(cur, sys, dice);
  }
  template <size_t index>
  auto _executeHelper(size_t cur, System &sys, rate_type dice,
                      typename std::enable_if<index == 0>::type *dummy = 0) {
    return std::get<0>(_procs).execute(sys, dice);
  }

  template <size_t... index>
  auto _cumulativeHelper(System &sys, std::index_sequence<index...>) {
    _rateSum = {{std::get<index>(_procs).rate(sys)...}};
    for (size_t i = 1; i < size; ++i) { // calculating cumulative hence i = 1
      _rateSum[i] += _rateSum[i - 1];
    }
  }
  std::tuple<Procs...> _procs;
  std::array<rate_type, size> _rateSum;
};

} // !namespace gkmc

#endif //

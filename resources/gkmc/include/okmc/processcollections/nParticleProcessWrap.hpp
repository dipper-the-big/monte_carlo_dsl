// nParticleProcessWrap.hpp

#ifndef GKMC_nParticleProcessWrap_H
#define GKMC_nParticleProcessWrap_H

namespace gkmc {

template <class Proc> class nParticleProcessWrap {
public:
  using System = typename Proc::System;
  using ParticleType = typename Proc::ParticleType;
  using rate_type = typename Proc::System::value_type;

  nParticleProcessWrap(Proc proc) : _proc{proc} {}

  auto execute(System &sys, rate_type dice) {
    auto temp = (dice / _proc.ratePerParticle(sys));
    auto i = (size_t)(temp);
    auto subDice = temp - i;
    return _proc.executeSelected(i, sys, subDice);
  }

  auto rate(System &sys) {
    return _proc.ratePerParticle(sys) * totalSystemTuple(sys);
  }
  
  auto totalSystemTuple(System &sys) {
    return sys.template particles<ParticleType>(_proc.species()).size();
  }
  const auto &proc() const { return _proc; }

private:
  Proc _proc;
};

template <class Proc, class Recalc> class ParticleRerateWrapper {
public:
  using System = typename Proc::System;
  using ParticleType = typename Proc::ParticleType;

  ParticleRerateWrapper(Proc proc, Recalc recalc)
      : _proc{proc}, _recalc{recalc} {}

  double rate(System &sys) {
    cumulativeRates(sys);
    if (_rateSum.size() == 0) return 0.0;
    return _rateSum[_rateSum.size() - 1];
  }

  const std::vector<double> &cumulativeRates(System &sys) {
    auto it = sys.template begin<ParticleType>(_proc.species());
    auto last = sys.template end<ParticleType>(_proc.species());
    auto &ps = sys.template particles<ParticleType>(_proc.species()); // TODO: clean
    double cumRate = 0.0;
    _rateSum.resize(ps.size());
    size_t i = 0;
    for (; it != last; ++it) {
      cumRate += _proc.ratePerParticle(sys) * _recalc(*it);
      _rateSum[i++] = cumRate;
    }
    return _rateSum;
  }

  auto executeSelected(size_t i, System &sys, double dice) {
    auto prevCumulativeRate = 0.0;
    if (i != 0)
      prevCumulativeRate = _rateSum[i - 1];
    auto subDice = dice - prevCumulativeRate;
    return _proc.executeSelected(i, sys, subDice);
  }

  auto execute(System &sys, const double &dice) {
    auto pickedCumulative =
        std::lower_bound(begin(_rateSum), end(_rateSum), dice);
    auto indexParticle = std::distance(begin(_rateSum), pickedCumulative);
    return executeSelected(indexParticle, sys, dice);
  }

  const auto &proc() const { return _proc; }

private:
  std::vector<double> _rateSum;
  Proc _proc;
  Recalc _recalc;
};
} // !namespace gkmc

#endif //

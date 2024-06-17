// SystemProcessWrap.hpp

#ifndef GKMC_SystemProcessWrap_H
#define GKMC_SystemProcessWrap_H

#include <random>

namespace gkmc {

template <class Proc>
class SystemProcessWrap {
public:
  using System = typename Proc::System;

  SystemProcessWrap(Proc proc) : _proc{std::forward<Proc>(proc)} { }

  auto isDecide(bool f) { _isDecide = f; }
  auto isDecide() const { return _isDecide; }

  auto executeFor(System& sys, double deltaT) {
    std::random_device rd;
    std::mt19937 gen(rd());
    auto d = std::poisson_distribution<> {deltaT * _proc.rate(sys)};
    auto times = d(gen);
    while (times--) {
      _proc.execute(sys);
    }
  }
 
  auto orRate(System& sys) { return _proc.rate(sys); }

  const auto &proc() const { return _proc; }
private:
  Proc _proc;
  bool _isDecide = true;
};

} // !namespace gkmc

#endif //

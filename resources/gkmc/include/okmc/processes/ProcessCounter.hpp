// ProcessCounter.hpp

#ifndef GKMC_ProcessCounter_H
#define GKMC_ProcessCounter_H

namespace gkmc {
template <class Proc> class ProcessCounter {
public:
  using System = typename Proc::System;
  using rate_type = typename System::value_type;
  ProcessCounter(Proc process) : _process{process} {}
  void execute(System &sys, rate_type dice) {
    ++_count;
    _process.execute(sys, dice);
  }
  decltype(auto) rate(System &s) {
    return _process.rate(s);
  }
  const auto &count() const { return _count; }
  const auto &proc() const { return _process; }
private:
  Proc _process;
  size_t _count = 0;
};
}

#endif
// AndProcess.hpp

#ifndef GKMC_AndProcess_H
#define GKMC_AndProcess_H

#include <algorithm>
#include <vector>
#include <memory>

namespace gkmc {

/**
 * Collection for the processes that execute independently in a single time-step
 * */
template <class Proc>
class AndProcess {
public:
  using System = typename Proc::System;

  void add(Proc proc) { _procs.emplace_back(std::forward<Proc>(proc)); }

  size_t size() const { return _procs.size(); }

  const auto& procs() const { return _procs; }

  auto maxAndRate(System& sys) {
    auto max = 0.0;
    for (auto &it : _procs) {
      auto temp = it.orRate(sys);
      if (temp > max) max = temp;
    }
    return max;
  }

  auto executeFor(System& sys, double deltaT) {
    for (auto &proc : _procs) {
      proc.executeFor(sys, deltaT); 
    }
  }
protected:
  std::vector<Proc> _procs;
};

template <class Sys>
class ISysProcess {
public:
  using System = Sys;
  using rate_type = typename System::value_type;

  virtual void executeFor(System &sys, rate_type delT) = 0;
  virtual const rate_type orRate(System &sys) = 0;
  virtual bool isDecide() const = 0;
};

template <class P>
class VirtualizeSysProc : public ISysProcess<typename P::System> {
public:
  using System = typename P::System;
  using rate_type = typename System::value_type;
  VirtualizeSysProc(P proc) : _proc{std::forward<P>(proc)} {}
  void executeFor(System &sys, rate_type dice) override {
    _proc.executeFor(sys, dice);  
  }
  const rate_type orRate(System &sys) override { return _proc.orRate(sys); }
  bool isDecide() const override { return _proc.isDecide(); }
private:
  P _proc;
};

template <class Sys>
class AndProcessVirtual {
public:
  using System = Sys;

  template <class Proc>
  void add(Proc p) { 
    _procs.emplace_back(std::make_unique<VirtualizeSysProc<Proc>>(std::forward<Proc>(p))); 
  }

  size_t size() const { return _procs.size(); }

  const auto& procs() const { return _procs; }

  auto maxAndRate(System& sys) {
    auto max = 0.0;
    for (auto &proc : _procs) {
      if (!proc->isDecide()) continue;
      auto temp = proc->orRate(sys);
      if (temp > max) max = temp;
    }
    return max;
  }

  auto executeFor(System& sys, double deltaT) {
    for (auto &proc : _procs) {
      proc->executeFor(sys, deltaT); 
    }
  }
protected:
  std::vector<std::unique_ptr<ISysProcess<Sys>>> _procs;
};

} // !namespace gkmc

#endif //

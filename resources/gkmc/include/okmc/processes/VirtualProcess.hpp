// VirtualProcess.hpp

#ifndef GKMC_VirtualProcess_H
#define GKMC_VirtualProcess_H

#include <algorithm>
#include <array>
#include <assert.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include <helper/meta.hpp>

namespace gkmc {

/**
 * interface for runtime polymorphism in Process
 * */
template <class Sys>
class IProcess {
public:
  using System = Sys;
  using rate_type = typename System::value_type;
  virtual void execute(System &sys, rate_type dice) = 0;
  virtual const rate_type rate(System &sys) = 0;
};

/**
 * interface for prototype based runtime polymorphism as
 * described by Sean. The Process class does not need to
 * inherit from an interface but using this can be added
 * in a single container with other processes.
 * */
template <class P, class System = typename P::System>
class Virtualize : public IProcess<System> {
public:
  using rate_type = typename System::value_type;
  Virtualize(P proc) : _proc{proc} {}
  void execute(System &sys, rate_type dice) override {
    _proc.execute(sys, dice);  
  }
  const rate_type rate(System &sys) override { return _proc.rate(sys); }
private:
  P _proc;
};

} // namespace gkmc

#endif
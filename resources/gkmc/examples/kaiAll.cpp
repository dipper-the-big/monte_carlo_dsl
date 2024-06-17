#include <array>
#include <chrono>
#include <iostream>
#include <random>

#include <bkl.hpp>
#include <helper/RandReal.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/ProcessCounter.hpp>
#include <okmc/reactions/Annihilation.hpp>
#include <okmc/reactions/ReactionInactiveCounter.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/particles/SystemTuple.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/SystemProps.hpp>
#include <okmc/hooks/Eager.hpp>
#include <okmc/hooks/Both.hpp>
#include <okmc/processes/VirtualProcess.hpp>
#include <okmc/reactions/ReactionCollection.hpp>
#include <okmc/particles/SystemSingle.hpp>
#include <okmc/particles/SystemPtr.hpp>
#include <okmc/predicates/benchPreds.hpp>

template <class Gen> auto givePoint(Gen &gen) {
  using Coords = std::array<double, 3>;
  Coords c;
  std::generate(begin(c), end(c), gen);
  return c;
}

enum class Species : int { I, V };

namespace probInvars {
  constexpr double ws[] = {1.717e15, 0.001282e15};
  constexpr double ems[] = {1.37, 0.1};
  constexpr int counts[] = {150, 150};
  constexpr double initRadius[] = {60.0, 20.0};
  constexpr double jump = 2.35;
}

struct CapRadius {
  template <class Sys, class Particle>
  double operator()(const Sys &, const Particle &, const Particle &,
                       const typename Sys::Pid &, const typename Sys::Pid &) {
  return 4.0; }
};

using namespace gkmc;

using Props = SystemProps<3, Species>;

template <class Props, class Particle>
struct SysHook {
  template <class Hook>
  using Sys = SystemTuple<Props, Hook, Particle>;
};

template <class Props, class Particle>
struct SysHookVirtual {
  template <class Hook>
  using Sys = SystemPtr<Props, Hook, Particle>;
};

auto sysProps(double temperature) {
  auto props = Props{};
  props.temperature(temperature);
  return props;
}

auto neighBasic() { return NeighBasic{true}; }

template <class Particle, class Sys> auto processes(Sys &sys) {
  using Proc = Jump3DRandom<Sys, Particle>;
  using ProcWrap = nParticleProcessWrap<Proc>;
  using ProcCounter = ProcessCounter<ProcWrap>;
  using Procs = VectorProcessCollection<Sys, ProcCounter>;
  using namespace probInvars;
  Procs procs;
  procs.add(ProcCounter{ProcWrap{Proc{Species::I, ems[0], ws[0], jump, sys.props().temperature()}}});
  procs.add(ProcCounter{ProcWrap{Proc{Species::V, ems[1], ws[1], jump, sys.props().temperature()}}});
  return procs;
} 

template <class Particle, class Sys> auto virtualProcesses(Sys &sys) {
  using Proc = Jump3DRandom<Sys, Particle>;
  using ProcWrap = nParticleProcessWrap<Proc>;
  using ProcCounter = ProcessCounter<ProcWrap>;
  using Procs = VectorProcessCollection<Sys, std::unique_ptr<IProcess<Sys>>>;
  using namespace probInvars;
  Procs procs;
  procs.add(std::make_unique<Virtualize<ProcCounter>>(ProcCounter{ProcWrap{Proc{Species::I, ems[0], ws[0], jump, sys.props().temperature()}}}));
  procs.add(std::make_unique<Virtualize<ProcCounter>>(ProcCounter{ProcWrap{Proc{Species::V, ems[1], ws[1], jump, sys.props().temperature()}}}));
  return procs;
}

template <class Particle, class Neigh> auto reactions(Neigh &neigh) {
  using Reaction = Reflexive<Annihilation, Neigh, CapRadius, Species, Particle>;
  return Reaction{Annihilation{}, Species::I, Species::V, neigh, CapRadius{}};
}

template <class Props, class Particle, class Reaction>
auto dirty(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, Particle>;
  auto dirt = Dirt{reaction};
  return dirt;
}

template <class Props, class Reaction> 
auto eager(Reaction &reaction) {
  using Hook = Eager<Props, Reaction &>;
  return Hook{reaction};
}

template <class Uno, class Dos> auto both(Uno &uno, Dos &dos) {
  return Both<Uno &, Dos &>{uno, dos};
}

template <class Particle, class Props, class Hook>
auto system(Props &props, Hook &hook) {
  using Sys = SystemTuple<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Species::V, "V");
  sys.template addSpecies<Particle>(Species::I, "I");
  return sys;
}

template <class Particle, class Props, class Hook>
auto systemUniquePtr(Props &props, Hook &hook) {
  using Sys = SystemPtr<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Species::V, "V");
  sys.template addSpecies<Particle>(Species::I, "I");
  return sys;
}

template <class Particle, class Props, class Hook>
auto systemSingleType(Props &props, Hook &hook) {
  using Sys = SystemSingle<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Species::V, "V");
  sys.template addSpecies<Particle>(Species::I, "I");
  return sys;
}

template <class Particle, class Sys>
auto fillSystemTuple(Sys &sys) {
  using std::normal_distribution;
  using namespace probInvars;
  for (auto _ = 0; _ <= counts[0]; ++_) {
    RandReal<normal_distribution<>> rnd{-initRadius[0],
                                        initRadius[0]}; // TODO: outside
    auto particle = Particle{givePoint(rnd)};
    sys.template add<Particle>(Species::I, particle);
  }
  for (auto _ = 0; _ <= counts[1]; ++_) {
    RandReal<normal_distribution<>> rnd{-initRadius[1],
                                        initRadius[1]}; // TODO: outside
    auto particle = Particle{givePoint(rnd)};
    sys.template add<Particle>(Species::V, particle);
  }
  sys.update(0.0,0.0);
}

template <class Sys, class Particle> struct Logger {
  Logger(size_t dumpInterval) : _dumpInterval{dumpInterval} {}
  auto operator()(Sys &p, double time, size_t steps) {
    using std::chrono::system_clock;
    if (!steps || steps % _dumpInterval) return;
    std::chrono::duration<double> elapsed = system_clock::now() - start;
    auto isize = p.template particles<Particle>(Species::I).size();
    auto vsize = p.template particles<Particle>(Species::V).size();
    std::cout << elapsed.count() << '\t' << time << '\t' << steps << '\t'
              << isize << '\t' << vsize << '\n';
    start = system_clock::now();
  }
private:
  std::chrono::system_clock::time_point start =
      std::chrono::system_clock::now();
  size_t _dumpInterval;
};

template <class Sys, class Particle>
auto predicate(size_t maxSteps, size_t dumpInterval) {
  using Log = Logger<Sys, Particle>;
  return StepsPredicate<Sys, Log>{maxSteps, Log{dumpInterval}};
}

auto runBasic(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = system<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count(); 
}

auto runEager(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = eager<Props>(reacts);
  auto sys = system<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count(); 
}

auto runVirtualProcess(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = system<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  virtualProcesses<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
    // printSystemTuple(sys);
  }
  return elapsed.count();
}

auto runVirtualReaction(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle>;
  auto hook = Hook{};
  using Reaction = Reflexive<Annihilation, NeighBasic, CapRadius, Species, Particle>;
  auto reaction = Reaction{Annihilation{}, Species::I, Species::V, neigh, CapRadius{}};
  hook.addReaction(reaction, Species::I);
  hook.addReaction(reaction, Species::V);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
    // printSystemTuple(sys);
  }
  return elapsed.count();
}

auto runVirtualReactionProcess(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle>;
  auto hook = Hook{};
  using Reaction = Reflexive<Annihilation, NeighBasic, CapRadius, Species, Particle>;
  auto reaction = Reaction{Annihilation{}, Species::I, Species::V, neigh, CapRadius{}};
  hook.addReaction(reaction, Species::I);
  hook.addReaction(reaction, Species::V);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
  fillSystemTuple<Particle>(sys);
  auto procs =  virtualProcesses<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count();
}

auto runSingleType(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = systemSingleType<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count(); 
}

auto runUniquePtr(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = systemUniquePtr<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count(); 
}

auto runUniquePtrEager(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticle<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = eager<Props>(reacts);
  auto sys = systemUniquePtr<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count(); 
}

auto runUniquePtrVirtual(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticleVirtual<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = systemUniquePtr<Particle>(props, hook);
  fillSystemTuple<Particle>(sys);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count(); 
}

auto runAllVirtual(double temperature, size_t maxSteps, size_t dumpInterval) {
  using Particle = SimpleParticleVirtual<typename Props::Coords>;
  auto props = sysProps(temperature);
  NeighBasic neigh = neighBasic();
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle>;
  auto hook = Hook{};
  using Reaction = Reflexive<Annihilation, NeighBasic, CapRadius, Species, Particle>;
  auto reaction = Reaction{Annihilation{}, Species::I, Species::V, neigh, CapRadius{}};
  hook.addReaction(reaction, Species::I);
  hook.addReaction(reaction, Species::V);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
  fillSystemTuple<Particle>(sys);
  auto procs =  virtualProcesses<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval);
  auto rnd = RandReal<>{0.0, 1.0};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  }
  return elapsed.count();
}

auto timeAll(double temperature, size_t maxSteps, size_t dumpInterval) {
  std::cout << "\n==============\n";
  auto t = runBasic(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in Basic: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runEager(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in Eager: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualProcess(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in VirtualProcess: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualReaction(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in VirtualReaction: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualReactionProcess(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in VirtualReactionProcess: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runSingleType(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in SingleType: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtr(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in UniquePtr: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtrEager(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in UniquePtrEager: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtrVirtual(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in UniquePtrVirtual: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runAllVirtual(temperature, maxSteps, dumpInterval);
  std::cout << "total time takes in AllVirtual: " << t;
  std::cout << "\n==============\n";
}

int main(int argc, char *argv[]) {
  const auto temperature = std::array<double, 3>{{500, 1500, 2500}};
  const auto maxSteps = std::array<size_t, 3>{{8000000, 31000000, 170000000}};
  for (auto index = 0; index < 3; ++index) {
    for (auto i = 0; i < 5; ++i) {
      timeAll(temperature[index], maxSteps[index] + 1, maxSteps[index] / 10);
    }
  }
  return 0;
}

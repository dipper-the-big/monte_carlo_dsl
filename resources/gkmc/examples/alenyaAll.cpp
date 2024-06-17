#include <chrono>
#include <iostream>
#include <random>
#include <string>

#include <bkl.hpp>
#include <helper/RandReal.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/Eager.hpp>
#include <okmc/hooks/NeighCell2D.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/hooks/Both.hpp>
#include <okmc/processcollections/TupleProcessCollection.hpp>
#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <okmc/processes/Jump2DRandom.hpp>
#include <okmc/processes/SpawnSingle.hpp>
#include <okmc/processes/Remove.hpp>
#include <okmc/processes/VirtualProcess.hpp>
#include <okmc/reactions/Fallback.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/reactions/pbc.hpp>
#include <okmc/reactions/ReactionCollection.hpp>
#include <okmc/particles/SystemTuple.hpp>
#include <okmc/particles/SystemSingle.hpp>
#include <okmc/particles/SystemPtr.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/predicates/benchPreds.hpp>
#include <okmc/SystemProps.hpp>

enum class Species : int { H, H2 };

using real_type = double;

using namespace gkmc;

using Props = SystemProps<2, Species, real_type>;

namespace probInvars {
constexpr real_type temperature = 500;
constexpr real_type boxDim = 1000.0;
constexpr real_type wForAll = 1e13;
constexpr real_type spawnRateH = 1e24 * 1e-20 * boxDim * boxDim;
constexpr real_type emJumpH = 0.9;
constexpr real_type lenJumpH = 34.4;
constexpr real_type emDestroyH = 1.9;
constexpr real_type emDestroyH2 = 0.06;
constexpr size_t maxSteps = 1e6;
}

auto sysProps() {
  auto p = Props{};
  p.temperature(probInvars::temperature);
  p.boxDimensions(
  Props::Coords{{probInvars::boxDim, probInvars::boxDim}});
  return p;
}

template <class Particle> struct Recombine {
  template <class Sys>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, Particle &p, Particle &n) {
    sys.template remove<Particle>(pid);
    sys.template remove<Particle>(nid);
    sys.template add<Particle>(Species::H2, p);
    return true;
  }
};

template <class Particle> struct RecombineVirtual {
  template <class Sys>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, SimpleParticleVirtual<typename Props::Coords, 1> &p, 
                  SimpleParticleVirtual<typename Props::Coords, 1> &n) {
    sys.template remove<Particle>(pid);
    sys.template remove<Particle>(nid);
    sys.template add<Particle>(Species::H2, (Particle&)p);
    return true;
  }
};

struct Dist {
  template <class Sys, class Particle>
  real_type operator()(const Sys &, const Particle &, const Particle &,
                       const typename Sys::Pid &, const typename Sys::Pid &) {
    return 2.0;
  }
};

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

auto neighBasic() { return NeighBasic{true}; }

template <class Props> auto neighCell(const Props &p, double cellLen) {
  return NeighCell2D<Props>{p.boxDimensions(), cellLen, true};
}

template <class Particle, class Sys> auto processes(Sys &sys) {
  using Spawn = SpawnSingle<Particle, Species>;
  using Jump = Jump2DRandom<Sys, Particle>;
  using Destroy = Remove<Sys, Particle>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using DestroyWrap = nParticleProcessWrap<Destroy>;
  using Destroys = VectorProcessCollection<Sys, DestroyWrap>;
  using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap, Destroys>;
  using namespace probInvars;
  auto spawnH = Spawn{spawnRateH, Species::H};
  auto jumpH = JumpWrap{Jump{sys, Species::H, emJumpH, wForAll, lenJumpH}};
  auto destroyH = DestroyWrap{
      Destroy{sys.props().temperature(), Species::H, emDestroyH, wForAll}};
  auto destroyH2 = DestroyWrap{
      Destroy{sys.props().temperature(), Species::H2, emDestroyH2, wForAll}};
  Destroys ds{};
  ds.add(destroyH2);
  ds.add(destroyH);
  return Procs{spawnH, jumpH, ds};
}

template <class Particle, class Sys> auto virtualProcesses(Sys &sys) {
  using Spawn = SpawnSingle<Particle, Species>;
  using Jump = Jump2DRandom<Sys, Particle>;
  using Destroy = Remove<Sys, Particle>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using DestroyWrap = nParticleProcessWrap<Destroy>;
  using Procs = VectorProcessCollection<Sys, std::unique_ptr<IProcess<Sys>>>;
  using namespace probInvars;
  auto spawnH = std::make_unique<Virtualize<Spawn>>(Spawn{spawnRateH, Species::H});
  auto jumpH = std::make_unique<Virtualize<JumpWrap>>(JumpWrap{Jump{sys, Species::H, emJumpH, wForAll, lenJumpH}});
  auto destroyH = std::make_unique<Virtualize<DestroyWrap>>(DestroyWrap{
      Destroy{sys.props().temperature(), Species::H, emDestroyH, wForAll}});
  auto destroyH2 = std::make_unique<Virtualize<DestroyWrap>>(DestroyWrap{
      Destroy{sys.props().temperature(), Species::H2, emDestroyH2, wForAll}});
  Procs procs{};
  procs.add(std::move(spawnH));
  procs.add(std::move(jumpH));
  procs.add(std::move(destroyH2));
  procs.add(std::move(destroyH));
  return procs;
}

template <class Particle, class Neigh> auto reactions(Neigh &neigh) {
  using Add = Reflexive<Recombine<Particle>, Neigh &, Dist, Species, Particle>;
  using Reaction = Fallback<PBC, Add>;
  using ReactionSp = SpeciesSelect<Reaction, NoReact, Species>;
  auto add = Add{Recombine<Particle>{}, Species::H, neigh, Dist{}};
  return ReactionSp{Reaction{PBC{}, add}, NoReact{}, Species::H};
}

template <class Particle, class Neigh> auto virtualTypeReactions(Neigh &neigh) {
  using Add = Reflexive<RecombineVirtual<Particle>, Neigh &, Dist, Species, SimpleParticleVirtual<typename Props::Coords, 1>>;
  using Reaction = Fallback<PBC, Add>;
  using ReactionSp = SpeciesSelect<Reaction, NoReact, Species>;
  auto add = Add{RecombineVirtual<Particle>{}, Species::H, neigh, Dist{}};
  return ReactionSp{Reaction{PBC{}, add}, NoReact{}, Species::H};
}


template <class Props, class Particle, class Reaction>
auto dirty(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, Particle>;
  auto dirt = Dirt{reaction};
  dirt.blackList(Species::H2);
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
  sys.template addSpecies<Particle>(Species::H, "H");
  sys.template addSpecies<Particle>(Species::H2, "H2");
  return sys;
}

template <class Particle, class Props, class Hook>
auto systemUniquePtr(Props &props, Hook &hook) {
  using Sys = SystemPtr<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Species::H, "H");
  sys.template addSpecies<Particle>(Species::H2, "H2");
  return sys;
}

template <class Particle, class Props, class Hook>
auto systemSingleType(Props &props, Hook &hook) {
  using Sys = SystemSingle<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Species::H, "H");
  sys.template addSpecies<Particle>(Species::H2, "H2");
  return sys;
}

template <class Sys, class Particle> struct Logger {
  Logger(size_t dumpInterval) : _dumpInterval{dumpInterval} {}
  auto operator()(Sys &p, double time, size_t steps) {
    using std::chrono::system_clock;
    if (!steps || steps % _dumpInterval) return;
    std::chrono::duration<double> elapsed = system_clock::now() - start;
    auto& h = p.template particles<Particle>(Species::H);
    auto hsize = h.size();
    auto h2size = p.template particles<Particle>(Species::H2).size();
    std::cout << elapsed.count() << '\t' << time << '\t' << steps << '\t'
              << hsize << '\t' << h2size << '\n';
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

auto runBasic(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = system<Particle>(props, hook);
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

auto runCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions<Particle>(neigh);
  auto dirt = dirty<Props, Particle>(reacts);
  auto hook = both(dirt, neigh);
  auto sys = system<Particle>(props, hook);
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

auto runEager(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = eager<Props>(reacts);
  auto sys = system<Particle>(props, hook);
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

auto runEagerCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions<Particle>(neigh);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);  // TODO: why the opposite is taking too long
  auto sys = system<Particle>(props, hook);
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

auto runVirtualProcess(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = system<Particle>(props, hook);
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

auto runVirtualProcessEagerCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions<Particle>(neigh);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = system<Particle>(props, hook);
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

auto runVirtualReaction(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle>;
  auto hook = Hook{};
  using Add = Reflexive<Recombine<Particle>, NeighBasic &, Dist, Species, Particle>;
  auto add = Add{Recombine<Particle>{}, Species::H, neigh, Dist{}};
  hook.addReaction(PBC{}, Species::H);
  hook.addReaction(add, Species::H);
  hook.first().blackList(Species::H2);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
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

auto runVirtualReactionCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle, NeighCell2D<Props>&>;
  auto hook = Hook{neigh};
  using Add = Reflexive<Recombine<Particle>, NeighCell2D<Props>&, Dist, Species, Particle>;
  auto add = Add{Recombine<Particle>{}, Species::H, neigh, Dist{}};
  hook.addReaction(PBC{}, Species::H);
  hook.addReaction(add, Species::H);
  hook.first().blackList(Species::H2);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
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

auto runVirtualReactionProcess(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle>;
  auto hook = Hook{};
  using Add = Reflexive<Recombine<Particle>, NeighBasic &, Dist, Species, Particle>;
  auto add = Add{Recombine<Particle>{}, Species::H, neigh, Dist{}};
  hook.addReaction(PBC{}, Species::H);
  hook.addReaction(add, Species::H);
  hook.first().blackList(Species::H2);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
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

auto runVirtualReactionProcessCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  using CurSysHook = SysHook<Props, Particle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, Particle, NeighCell2D<Props>&>;
  auto hook = Hook{neigh};
  using Add = Reflexive<Recombine<Particle>, NeighCell2D<Props>&, Dist, Species, Particle>;
  auto add = Add{Recombine<Particle>{}, Species::H, neigh, Dist{}};
  hook.addReaction(PBC{}, Species::H);
  hook.addReaction(add, Species::H);
  hook.first().blackList(Species::H2);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
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

auto runSingleType(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = systemSingleType<Particle>(props, hook);
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

auto runSingleTypeEagerCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions<Particle>(neigh);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = systemSingleType<Particle>(props, hook);
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

auto runUniquePtr(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = dirty<Props, Particle>(reacts);
  auto sys = systemUniquePtr<Particle>(props, hook);
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

auto runUniquePtrCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions<Particle>(neigh);
  auto dirt = dirty<Props, Particle>(reacts);
  auto hook = both(dirt, neigh);
  auto sys = systemUniquePtr<Particle>(props, hook);
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

auto runUniquePtrEager(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  auto hook = eager<Props>(reacts);
  auto sys = systemUniquePtr<Particle>(props, hook);
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

auto runUniquePtrEagerCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions<Particle>(neigh);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag); // TODO: why runtime err in other way
  auto sys = systemUniquePtr<Particle>(props, hook);
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

auto runUniquePtrVirtual(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticleVirtual<typename Props::Coords, 1>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords, 1>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto reacts = virtualTypeReactions<Particle>(neigh);
  auto hook = dirty<Props, VirtualParticle>(reacts);
  auto sys = systemUniquePtr<VirtualParticle>(props, hook);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), VirtualParticle>(maxSteps, dumpInterval);
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

auto runUniquePtrVirtualEagerCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticleVirtual<typename Props::Coords, 1>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords, 1>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto reacts = virtualTypeReactions<Particle>(neigh);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = systemUniquePtr<VirtualParticle>(props, hook);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), VirtualParticle>(maxSteps, dumpInterval);
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

auto runAllVirtual(size_t maxSteps, size_t dumpInterval) {
  using Particle = PbcParticle<SimpleParticleVirtual<typename Props::Coords, 1>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords, 1>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  using CurSysHook = SysHookVirtual<Props, VirtualParticle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, VirtualParticle>;
  auto hook = Hook{};
  using Add = Reflexive<RecombineVirtual<Particle>, NeighBasic &, Dist, Species, VirtualParticle>;
  auto add = Add{RecombineVirtual<Particle>{}, Species::H, neigh, Dist{}};
  hook.addReaction(PBC{}, Species::H);
  hook.addReaction(add, Species::H);
  hook.first().blackList(Species::H2);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
  auto procs =  virtualProcesses<Particle>(sys);
  auto pred = predicate<decltype(sys), VirtualParticle>(maxSteps, dumpInterval);
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

auto runAllVirtualCell(size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticleVirtual<typename Props::Coords, 1>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords, 1>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  using CurSysHook = SysHookVirtual<Props, VirtualParticle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, VirtualParticle, NeighCell2D<Props>&>;
  auto hook = Hook{neigh};
  using Add = Reflexive<RecombineVirtual<Particle>, NeighCell2D<Props> &, Dist, Species, VirtualParticle>;
  auto add = Add{RecombineVirtual<Particle>{}, Species::H, neigh, Dist{}};
  hook.addReaction(PBC{}, Species::H);
  hook.addReaction(add, Species::H);
  hook.first().blackList(Species::H2);
  using Sys = Hook::Sys;
  auto sys = Sys{props, hook};
  auto procs =  virtualProcesses<Particle>(sys);
  auto pred = predicate<decltype(sys), VirtualParticle>(maxSteps, dumpInterval);
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

auto timeAll(double cellLen) {
  constexpr auto maxSteps = probInvars::maxSteps;
  constexpr auto maxStepsBig = probInvars::maxSteps * 4;
  constexpr auto dumpSteps = 50000;
  constexpr auto dumpStepsBig = dumpSteps * 4;

  std::cout << "\n==============\n";
  auto t = runBasic(maxSteps, dumpSteps);
  std::cout << "total time takes in Basic: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in Cell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runEager(maxSteps, dumpSteps);
  std::cout << "total time takes in Eager: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runEagerCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in EagerCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runVirtualProcess(maxSteps, dumpSteps);
  std::cout << "total time takes in VirtualProcess: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runVirtualProcessEagerCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in VirtualProcessEagerCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runVirtualReaction(maxSteps, dumpSteps);
  std::cout << "total time takes in VirtualReaction: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runVirtualReactionCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in VirtualReactionCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runVirtualReactionProcess(maxSteps, dumpSteps);
  std::cout << "total time takes in VirtualReactionProcess: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runVirtualReactionProcessCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in VirtualReactionProcessCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runSingleType(maxSteps, dumpSteps);
  std::cout << "total time takes in SingleType: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runSingleTypeEagerCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in SingleTypeEagerCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runUniquePtr(maxSteps, dumpSteps);
  std::cout << "total time takes in UniquePtr: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runUniquePtrCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in UniquePtrCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runUniquePtrEager(maxSteps, dumpSteps);
  std::cout << "total time takes in UniquePtrEager: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runUniquePtrEagerCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in UniquePtrEagerCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runUniquePtrVirtual(maxSteps, dumpSteps);
  std::cout << "total time takes in UniquePtrVirtual: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runUniquePtrVirtualEagerCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in UniquePtrVirtualEagerCell: " << t;
  std::cout << "\n==============\n";


  std::cout << "\n==============\n";
  t = runAllVirtual(maxSteps, dumpSteps);
  std::cout << "total time takes in AllVirtual: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runAllVirtualCell(maxStepsBig, dumpStepsBig, cellLen);
  std::cout << "total time takes in AllVirtualCell: " << t;
  std::cout << "\n==============\n";
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Please provide cell length as command line argument.\n";
    return 1;
  }
  auto cellLen = 50.0;
  try {
    cellLen = std::stod(argv[1]);
  } catch (...) {
    std::cerr << "couldn't convert first argument to double.\n";
    return 1;
  }
  timeAll(cellLen);
  return 0;
}

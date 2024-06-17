// HeAll.cpp

#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <string>

#include <bkl.hpp>
#include <helper/RandReal.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/Eager.hpp>
#include <okmc/hooks/NeighCell.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/hooks/Both.hpp>
#include <okmc/processcollections/TupleProcessCollection.hpp>
#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/SpawnSingle.hpp>
#include <okmc/processes/Remove.hpp>
#include <okmc/processes/Emission.hpp>
#include <okmc/processes/VirtualProcess.hpp>
#include <okmc/reactions/Fallback.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/reactions/pbc.hpp>
#include <okmc/reactions/ReactionCollection.hpp>
#include <okmc/reactions/Absorb.hpp>
#include <okmc/particles/SystemTuple.hpp>
#include <okmc/particles/SystemSingle.hpp>
#include <okmc/particles/SystemPtr.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/predicates/benchPreds.hpp>
#include <okmc/SystemProps.hpp>

enum class Species : int { He, HeC };

using namespace gkmc;

using Props = SystemProps<3, Species>;
using P1 = PbcParticle<SimpleParticle<Props::Coords>>;
using P2 = SizeParticle<PbcParticle<SimpleParticle<Props::Coords>>>;

namespace probInvars {
constexpr double spawnRate = 22;
constexpr double jumpEm = 1.257;
constexpr double jumpW = 2.0e12;
constexpr double jumpLen = 9.0;
constexpr double emitDist = 10.0;
constexpr double minDist = 5.0;
constexpr double boxDim = 1e4;
constexpr double temperature = 800;
}

struct Dist {
  struct ClusterRadiusHe {
    double operator()(std::size_t s) const {
      return pow(8.37 * pow(s, 1.02), 0.333333);
    }
    template<class Sys, class P>
    auto operator()(const Sys &, const P &p, const typename Sys::Pid &) {
      return this->operator()(p.size());
    }
  };
  Dist() {}
  using RadiusCalc = SizeBuffer<ClusterRadiusHe>;
  template<class Sys, class P, class N>
  double operator()(const Sys&, const P &p, const N &n, typename Sys::Pid, typename Sys::Pid) {
    return radiusCalc()(p.size()) + radiusCalc()(n.size()) + probInvars::minDist;
  };
  auto& radiusCalc() { return _sizeBuffer; }
private:
  RadiusCalc _sizeBuffer = RadiusCalc{ClusterRadiusHe{}};
};

struct DistNoBuffer {
  struct ClusterRadiusHe {
    double operator()(int s) const {
      return pow(8.37 * pow(s, 1.02), 0.333333);
    }
    template<class Sys, class P>
    auto operator()(const Sys &, const P &p, const typename Sys::Pid &) {
      return this->operator()(p.size());
    }

  };
  using RadiusCalc = ClusterRadiusHe;
  DistNoBuffer() {}
  template<class Sys, class P, class N>
  double operator()(const Sys&, const P &p, const N &n, typename Sys::Pid, typename Sys::Pid) {
    return radiusCalc()(p.size()) + radiusCalc()(n.size()) + probInvars::minDist;
  };
  auto& radiusCalc() { return _sizeCalc; }
private:
  ClusterRadiusHe _sizeCalc = ClusterRadiusHe{};
};

auto calcHeEmissionRate(size_t size, double temperature) {
  auto nrg = [&](int i) -> double { return 1.863 * std::pow(i, 0.8998425); };
  auto tm = [&](int i) -> double {
    return 5e10 * std::exp(-(nrg(i) - nrg(i - 1)) / (invars::kB * temperature));
  };
  return 1.0 / tm(size);
}

class HeEmissionRate {
public:
  HeEmissionRate(double temperature) : _temperature{temperature} {}
  auto operator()(size_t sz) const { return calcHeEmissionRate(sz, _temperature); }

private:
  double _temperature;
};

template <class Proc> class HeEmissionWrapper {
public:
  using System = typename Proc::System;
  using ParticleType = typename Proc::ParticleType;

  HeEmissionWrapper(Proc proc, double temperature)
      : _proc{proc}, _calcRate{HeEmissionRate{temperature}} {}

  double rate(System &sys) {
    cumulativeRates(sys);
    if (_rateSum.size() == 0)
      return 0.0;
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
      cumRate += _calcRate(sys, *it, sys.cookPid(_proc.species(), i));
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
  SizeBuffer<HeEmissionRate> _calcRate;
};

auto sysProps() {
  auto p = Props{};
  p.temperature(probInvars::temperature);
  p.boxDimensions(Props::Coords{
      {probInvars::boxDim, probInvars::boxDim, probInvars::boxDim}});
  return p;
}

template <class Neigh, class Cap> auto reactions(Neigh &neigh, Cap &cap) {
  using AbsorbT = All<AbsorbNew<P2, Species>, Neigh &, Cap &, P1, P2>;
  using TypePBC = TypeSelect<PBC, NoReact, P1>;
  using Reaction = Fallback<TypePBC, AbsorbT>;
  return Reaction{
      TypePBC{PBC{}, NoReact{}},
      AbsorbT{AbsorbNew<P2, Species>{Species::HeC, 2}, neigh, cap}};
}

template <class Neigh, class Cap> auto reactionsSingleType(Neigh &neigh, Cap &cap) {
  using AbsorbT = All<AbsorbNew<P2, Species>, Neigh &, Cap &, P2>;
  using SpeciesPBC = SpeciesSelect<PBC, NoReact, Species>;
  using Reaction = Fallback<SpeciesPBC, AbsorbT>;
  return Reaction{
      SpeciesPBC{PBC{}, NoReact{}, Species::He},
      AbsorbT{AbsorbNew<P2, Species>{Species::HeC, 2}, neigh, cap}};
}

template <class Neigh, class Cap> auto reactionsVirtual(Neigh &neigh, Cap &cap) {
  using Particle = SizeParticle<PbcParticle<SimpleParticleVirtual<typename Props::Coords>>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords>;
  using AbsorbT = All<AbsorbNew<Particle, Species>, Neigh &, Cap &, VirtualParticle>;
  using SpeciesPBC = SpeciesSelect<PBC, NoReact, Species>;
  using Reaction = Fallback<SpeciesPBC, AbsorbT>;
  return Reaction{
      SpeciesPBC{PBC{}, NoReact{}, Species::He},
      AbsorbT{AbsorbNew<Particle, Species>{Species::HeC, 2}, neigh, cap}};
}

template <class Sys, class Cap> auto processesVirtualType(Sys& sys, Cap &cap) {
  using T1 = PbcParticle<SimpleParticleVirtual<Props::Coords>>;
  using T2 = SizeParticle<PbcParticle<SimpleParticleVirtual<Props::Coords>>>;
  using Spawn = SpawnSingle<T1, Species>;
  using Jump = Jump3DRandom<Sys, T1>;
  using Emission = Emission<Sys, T2, T1, typename Cap::RadiusCalc &>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using EmissionWrap = HeEmissionWrapper<Emission>;
  using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap, EmissionWrap>;
  using namespace probInvars;
  auto spawn = Spawn{spawnRate, Species::He};
  auto jumpWrap =
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}};
  auto emission = EmissionWrap{
      Emission{cap.radiusCalc(), emitDist, 1.0, Species::HeC, Species::He, Species::He},
      sys.props().temperature()};
  return Procs{spawn, jumpWrap, emission};
}

template <class Sys, class Cap> auto processes(Sys& sys, Cap &cap) {
  using Spawn = SpawnSingle<P1, Species>;
  using Jump = Jump3DRandom<Sys, P1>;
  using Emission = Emission<Sys, P2, P1, typename Cap::RadiusCalc &>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using EmissionWrap = HeEmissionWrapper<Emission>;
  using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap, EmissionWrap>;
  using namespace probInvars;
  auto spawn = Spawn{spawnRate, Species::He};
  auto jumpWrap =
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}};
  auto emission = EmissionWrap{
      Emission{cap.radiusCalc(), emitDist, 1.0, Species::HeC, Species::He, Species::He},
      sys.props().temperature()};
  return Procs{spawn, jumpWrap, emission};
}

template <class Sys, class Cap> auto processesSingleType(Sys& sys, Cap &cap) {
  using Spawn = SpawnSingle<P2, Species>;
  using Jump = Jump3DRandom<Sys, P2>;
  using Emission = Emission<Sys, P2, P2, typename Cap::RadiusCalc &>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using EmissionWrap = HeEmissionWrapper<Emission>;
  using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap, EmissionWrap>;
  using namespace probInvars;
  auto spawn = Spawn{spawnRate, Species::He};
  auto jumpWrap =
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}};
  auto emission = EmissionWrap{
      Emission{cap.radiusCalc(), emitDist, 1.0, Species::HeC, Species::He, Species::He},
      sys.props().temperature()};
  return Procs{spawn, jumpWrap, emission};
}

template <class Sys, class Cap> auto virtualProcesses(Sys& sys, Cap &cap) {
  using Spawn = SpawnSingle<P1, Species>;
  using Jump = Jump3DRandom<Sys, P1>;
  using Emission = Emission<Sys, P2, P1, typename Cap::RadiusCalc &>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using EmissionWrap = HeEmissionWrapper<Emission>;
  using Procs = VectorProcessCollection<Sys, std::unique_ptr<IProcess<Sys>>>;
  using namespace probInvars;
  auto spawn = std::make_unique<Virtualize<Spawn, Sys>>(Spawn{spawnRate, Species::He});
  auto jumpWrap = std::make_unique<Virtualize<JumpWrap>>(
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}});
  auto emission = std::make_unique<Virtualize<EmissionWrap>>(EmissionWrap{
      Emission{cap.radiusCalc(), emitDist, 1.0, Species::HeC, Species::He, Species::He},
      sys.props().temperature()});
  Procs procs{};
  procs.add(std::move(spawn));
  procs.add(std::move(jumpWrap));
  procs.add(std::move(emission));
  return procs;
}

template <class Sys, class Cap> auto processesSingleTypeVirtual(Sys& sys, Cap &cap) {
  using Spawn = SpawnSingle<P2, Species>;
  using Jump = Jump3DRandom<Sys, P2>;
  using Emission = Emission<Sys, P2, P2, typename Cap::RadiusCalc &>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using EmissionWrap = HeEmissionWrapper<Emission>;
  using Procs = VectorProcessCollection<Sys, std::unique_ptr<IProcess<Sys>>>;
  using namespace probInvars;
  auto spawn = std::make_unique<Virtualize<Spawn, Sys>>(Spawn{spawnRate, Species::He});
  auto jumpWrap = std::make_unique<Virtualize<JumpWrap>>(
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}});
  auto emission = std::make_unique<Virtualize<EmissionWrap>>(EmissionWrap{
      Emission{cap.radiusCalc(), emitDist, 1.0, Species::HeC, Species::He, Species::He},
      sys.props().temperature()});
  Procs procs{};
  procs.add(std::move(spawn));
  procs.add(std::move(jumpWrap));
  procs.add(std::move(emission));
  return procs;
}

template <class Sys, class Cap> auto processesVirtualTypeVirtual(Sys& sys, Cap &cap) {
  using T1 = PbcParticle<SimpleParticleVirtual<Props::Coords>>;
  using T2 = SizeParticle<PbcParticle<SimpleParticleVirtual<Props::Coords>>>;
  using Spawn = SpawnSingle<T1, Species>;
  using Jump = Jump3DRandom<Sys, T1>;
  using Emission = Emission<Sys, T2, T1, typename Cap::RadiusCalc &>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using EmissionWrap = HeEmissionWrapper<Emission>;
  using Procs = VectorProcessCollection<Sys, std::unique_ptr<IProcess<Sys>>>;
  using namespace probInvars;
  auto spawn = std::make_unique<Virtualize<Spawn, Sys>>(Spawn{spawnRate, Species::He});
  auto jumpWrap = std::make_unique<Virtualize<JumpWrap>>(
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}});
  auto emission = std::make_unique<Virtualize<EmissionWrap>>(EmissionWrap{
      Emission{cap.radiusCalc(), emitDist, 1.0, Species::HeC, Species::He, Species::He},
      sys.props().temperature()});
  Procs procs{};
  procs.add(std::move(spawn));
  procs.add(std::move(jumpWrap));
  procs.add(std::move(emission));
  return procs;
}

auto neighBasic() { return NeighBasic{true}; }

template <class Props> auto neighCell(const Props &p, double cellLen) {
  return NeighCell<Props>{p.boxDimensions(), cellLen, true};
}

template <class Props, class P1, class P2, class Reaction>
auto dirty(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, P1, P2>;
  auto dirt = Dirt{reaction};
  return dirt;
}

template <class Props, class P, class Reaction>
auto dirtySingleType(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, P>;
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

template <class Hook>
auto system(Props &props, Hook &hook) {
  using Sys = SystemTuple<Props, Hook &, P1, P2>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<P1>(Species::He, "He");
  sys.template addSpecies<P2>(Species::HeC, "HeC");
  return sys;
}

template <class Hook>
auto systemSingleType(Props &props, Hook &hook) {
  using Sys = SystemSingle<Props, Hook &, P2>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<P2>(Species::He, "He");
  sys.template addSpecies<P2>(Species::HeC, "HeC");
  return sys;
}

template <class P, class Hook>
auto systemUniquePtr(Props &props, Hook &hook) {
  using Sys = SystemPtr<Props, Hook &, P>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<P>(Species::He, "He");
  sys.template addSpecies<P>(Species::HeC, "HeC");
  return sys;
}

template <class Sys, class P1, class P2> struct Logger {
  auto operator () (Sys &p, double t, size_t steps) {
    if ((int)t <= lastSecPrinted) return;
    using std::chrono::system_clock;
    std::chrono::duration<double> elapsed = system_clock::now() - start;
    lastSecPrinted = (int)t;
    pCounts[1] = p.template particles<P1>(Species::He).size();
    auto it = p.template begin<P2>(Species::HeC);
    auto last = p.template end<P2>(Species::HeC);
    for (;it != last; ++it) pCounts[it->size()]++;
    for (const auto &it : pCounts)
      std::cout << elapsed.count() << '\t' << t << '\t' << steps << '\t'
                << it.first << '\t' << it.second << '\n';
    pCounts.clear();
    // printSystemTuple(p);
    start = system_clock::now();
  }
private:
  std::chrono::system_clock::time_point start =
      std::chrono::system_clock::now();
  double lastSecPrinted = 0;
  std::unordered_map<int, int> pCounts;
};

template <class Sys, class P1, class P2>
auto predicate(double maxTime) {
  using Log = Logger<Sys, P1, P2>;
  return TimePredicate<Sys, Log>{maxTime, Log{}};
}

template <class Props, class Particle>
struct SysHook {
  template <class Hook>
  using Sys = SystemSingle<Props, Hook, Particle>;
};

template <class Props, class Particle>
struct SysHookVirtual {
  template <class Hook>
  using Sys = SystemPtr<Props, Hook, Particle>;
};

auto runBasic(double maxTime) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist);
  auto hook = dirty<Props, P1, P2>(reacts);
  auto sys = system(props, hook);
  auto procs = processes(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runNoBufferRadius(double maxTime) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = DistNoBuffer{};
  auto reacts = reactions(neigh, dist);
  auto hook = dirty<Props, P1, P2>(reacts);
  auto sys = system(props, hook);
  auto procs = processes(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runCell(double maxTime, double cellLen) {
  auto props = sysProps();
  auto dist = Dist{};
  auto neigh = neighCell(props, cellLen);
  auto reacts = reactions(neigh, dist);
  auto dirt = dirty<Props, P1, P2>(reacts);
  auto hook = both(dirt, neigh);
  auto sys = system(props, hook);
  auto procs = processes(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runEager(double maxTime) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist);
  auto hook = eager<Props>(reacts);
  auto sys = system(props, hook);
  auto procs = processes(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runEagerCell(double maxTime, double cellLen) {
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = system(props, hook);
  auto procs = processes(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runVirtualProcess(double maxTime) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist);
  auto hook = dirty<Props, P1, P2>(reacts);
  auto sys = system(props, hook);
  auto procs = virtualProcesses(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runVirtualProcessEagerCell(double maxTime, double cellLen) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = system(props, hook);
  auto procs = virtualProcesses(sys, dist);
  auto pred = predicate<decltype(sys), P1, P2>(maxTime);
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

auto runSingleType(double maxTime) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactionsSingleType(neigh, dist);
  auto hook = dirtySingleType<Props, P2>(reacts);
  auto sys = systemSingleType(props, hook);
  auto procs = processesSingleType(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runSingleTypeEagerCell(double maxTime, double cellLen) {
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto dist = Dist{};
  auto reacts = reactionsSingleType(neigh, dist);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = systemSingleType(props, hook);
  auto procs = processesSingleType(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runVirtualReaction(double maxTime) {
  auto props = sysProps();
  auto neigh = neighBasic();
  auto dist = Dist{};
  using CurSysHook = SysHook<Props, P2>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, P2>;
  auto hook = Hook{};
  using Recomb = All<AbsorbNew<P2, Species>, NeighBasic &, Dist &, P2>;
  auto recomb = Recomb{AbsorbNew<P2, Species>{Species::HeC, 2}, neigh, dist};
  hook.addReaction(PBC{}, Species::He);
  hook.addReaction(recomb);
  auto sys = systemSingleType(props, hook);
  auto procs = processesSingleType(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runVirtualReactionCell(double maxTime, double cellLen) {
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto dist = Dist{};
  using CurSysHook = SysHook<Props, P2>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, P2, NeighCell<Props> &>;
  auto hook = Hook{neigh};
  using Recomb = All<AbsorbNew<P2, Species>, NeighCell<Props> &, Dist &, P2>;
  auto recomb = Recomb{AbsorbNew<P2, Species>{Species::HeC, 2}, neigh, dist};
  hook.addReaction(PBC{}, Species::He);
  hook.addReaction(recomb);
  auto sys = systemSingleType(props, hook);
  auto procs = processesSingleType(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runVirtualReactionProcess(double maxTime) {
  auto props = sysProps();
  auto neigh = neighBasic();
  auto dist = Dist{};
  using CurSysHook = SysHook<Props, P2>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, P2>;
  auto hook = Hook{};
  using Recomb = All<AbsorbNew<P2, Species>, NeighBasic &, Dist &, P2>;
  auto recomb = Recomb{AbsorbNew<P2, Species>{Species::HeC, 2}, neigh, dist};
  hook.addReaction(PBC{}, Species::He);
  hook.addReaction(recomb);
  auto sys = systemSingleType(props, hook);
  auto procs = processesSingleTypeVirtual(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runUniquePtr(double maxTime) {
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactionsSingleType(neigh, dist);
  auto hook = dirtySingleType<Props, P2>(reacts);
  auto sys = systemUniquePtr<P2>(props, hook);
  auto procs = processesSingleType(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runUniquePtrEagerCell(double maxTime, double cellLen) {
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto dist = Dist{};
  auto reacts = reactionsSingleType(neigh, dist);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = systemUniquePtr<P2>(props, hook);
  auto procs = processesSingleType(sys, dist);
  auto pred = predicate<decltype(sys), P2, P2>(maxTime);
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

auto runUniquePtrVirtual(double maxTime) {
  // using Particle = SizeParticle<PbcParticle<SimpleParticleVirtual<typename Props::Coords>>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactionsVirtual(neigh, dist);
  auto hook = dirtySingleType<Props, VirtualParticle>(reacts);
  auto sys = systemUniquePtr<VirtualParticle>(props, hook);
  auto procs = processesVirtualType(sys, dist);
  auto pred = predicate<decltype(sys), VirtualParticle, VirtualParticle>(maxTime);
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

auto runUniquePtrVirtualEagerCell(double maxTime, double cellLen) {
  // using Particle = SizeParticle<PbcParticle<SimpleParticleVirtual<typename Props::Coords>>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto dist = Dist{};
  auto reacts = reactionsVirtual(neigh, dist);
  auto eag = eager<Props>(reacts);
  auto hook = both(neigh, eag);
  auto sys = systemUniquePtr<VirtualParticle>(props, hook);
  auto procs = processesVirtualType(sys, dist);
  auto pred = predicate<decltype(sys), VirtualParticle, VirtualParticle>(maxTime);
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

auto runAllVirtual(double maxTime) {
  using Particle = SizeParticle<PbcParticle<SimpleParticleVirtual<typename Props::Coords>>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords>;
  auto props = sysProps();
  NeighBasic neigh = neighBasic();
  auto dist = Dist{};
  using CurSysHook = SysHookVirtual<Props, VirtualParticle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, VirtualParticle>;
  auto hook = Hook{};
  using Recomb = All<AbsorbNew<Particle, Species>, NeighBasic &, Dist &, VirtualParticle>;
  auto recomb = Recomb{AbsorbNew<Particle, Species>{Species::HeC, 2}, neigh, dist};
  hook.addReaction(PBC{}, Species::He);
  hook.addReaction(recomb);
  auto sys = systemUniquePtr<VirtualParticle>(props, hook);
  auto procs = processesVirtualTypeVirtual(sys, dist);
  auto pred = predicate<decltype(sys), VirtualParticle, VirtualParticle>(maxTime);
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

auto runAllVirtualCell(double maxTime, double cellLen) {
  using Particle = SizeParticle<PbcParticle<SimpleParticleVirtual<typename Props::Coords>>>;
  using VirtualParticle = SimpleParticleVirtual<typename Props::Coords>;
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  auto dist = Dist{};
  using CurSysHook = SysHookVirtual<Props, VirtualParticle>;
  using Hook = ReactionCollection<CurSysHook::Sys, Props, VirtualParticle, NeighCell<Props>&>;
  auto hook = Hook{neigh};
  using Recomb = All<AbsorbNew<Particle, Species>, NeighCell<Props> &, Dist &, VirtualParticle>;
  auto recomb = Recomb{AbsorbNew<Particle, Species>{Species::HeC, 2}, neigh, dist};
  hook.addReaction(PBC{}, Species::He);
  hook.addReaction(recomb);
  auto sys = systemUniquePtr<VirtualParticle>(props, hook);
  auto procs = processesVirtualTypeVirtual(sys, dist);
  auto pred = predicate<decltype(sys), VirtualParticle, VirtualParticle>(maxTime);
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

auto timeAll(double maxTime, double cellLen) {
  std::cout << "\n==============\n";
  auto t = runAllVirtual(maxTime);
  std::cout << "total time takes in AllVirtual: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runNoBufferRadius(maxTime);
  std::cout << "total time takes in NoBufferRadius: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runAllVirtualCell(maxTime, cellLen);
  std::cout << "total time takes in AllVirtualCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualReactionCell(maxTime, cellLen);
  std::cout << "total time takes in VirtualReactionCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtr(maxTime);
  std::cout << "total time takes in UniquePtr: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtrEagerCell(maxTime, cellLen);
  std::cout << "total time takes in UniquePtrCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtrVirtual(maxTime);
  std::cout << "total time takes in UniquePtrVirtual: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runUniquePtrVirtualEagerCell(maxTime, cellLen);
  std::cout << "total time takes in UniquePtrVirtualEagerCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualReactionProcess(maxTime);
  std::cout << "total time takes in VirtualReactionProcess: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runSingleType(maxTime);
  std::cout << "total time takes in SingleType: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runSingleTypeEagerCell(maxTime, cellLen);
  std::cout << "total time takes in SingleTypeEagerCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runBasic(maxTime);
  std::cout << "total time takes in Basic: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runCell(maxTime, cellLen);
  std::cout << "total time takes in Cell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runEager(maxTime);
  std::cout << "total time takes in Eager: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runEagerCell(maxTime, cellLen);
  std::cout << "total time takes in EagerCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualProcess(maxTime);
  std::cout << "total time takes in VirtualProcess: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualProcessEagerCell(maxTime, cellLen);
  std::cout << "total time takes in VirtualProcessEagerCell: " << t;
  std::cout << "\n==============\n";

  std::cout << "\n==============\n";
  t = runVirtualReaction(maxTime);
  std::cout << "total time takes in VirtualReaction: " << t;
  std::cout << "\n==============\n";
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Please provide cell length and maxTime as command line argument.\n";
    return 1;
  }
  auto cellLen = 50.0;
  auto maxTime = 120.0;
  try {
    cellLen = std::stod(argv[1]);
    maxTime = std::stod(argv[2]);
  } catch (...) {
    std::cerr << "couldn't convert first argument to double.\n";
    return 1;
  }
  // runVirtualProcessEagerCell(maxTime);
  // runAllVirtualCell(maxTime, cellLen);
  // runBasic(maxTime);
  for (int i = 0; i < 3; ++i) timeAll(maxTime, cellLen);
  return 0;
}


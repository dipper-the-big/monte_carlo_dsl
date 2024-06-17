#include <chrono>
#include <iostream>
#include <random>
#include <string>

#include <pokmc.hpp>
#include <helper/RandReal.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/Eager.hpp>
#include <okmc/hooks/NeighCell2D.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/hooks/Both.hpp>
#include <pokmc/processcollections/OrProcess.hpp>
#include <pokmc/processcollections/AndProcess.hpp>
#include <pokmc/processes/ParticleProcessCounter.hpp>
#include <pokmc/processes/SystemProcessWrap.hpp>
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
constexpr real_type boxDim = 1000.0;
constexpr real_type wForAll = 1e13;
constexpr real_type spawnRateH = 1e24 * 1e-20 * boxDim * boxDim;
constexpr real_type emJumpH = 0.9;
constexpr real_type lenJumpH = 34.4;
constexpr real_type emDestroyH = 0.9;
constexpr real_type emDestroyH2 = 0.06;
constexpr size_t maxSteps = 1e6;
// constexpr auto maxSteps = 5000;
constexpr auto dumpInterval = 100;
}

auto sysProps(double temperature) {
  auto p = Props{};
  p.temperature(temperature);
  p.boxDimensions(
  Props::Coords{{probInvars::boxDim, probInvars::boxDim}});
  return p;
}

template <class Particle> struct RecombineH2 {
  template <class Sys>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, Particle &p, Particle &n) {
    sys.template remove<Particle>(pid);
    sys.template remove<Particle>(nid);
    sys.template add<Particle>(Species::H2, p);
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

auto neighBasic() { return NeighBasic{true}; }

template <class Props> auto neighCell(const Props &p, double cellLen) {
  return NeighCell2D<Props>{p.boxDimensions(), cellLen, true};
}

template <class Particle, class Sys> auto processes(Sys &sys) {
  using Spawn = SpawnSingle<Particle, Species>;
  using Jump = Jump2DRandom<Sys, Particle>;
  using Destroy = Remove<Sys, Particle>;
  using ProcsH = OrProcessVirtual<Sys, Particle>;
  using ProcsH2 = OrProcess<Destroy>;
  using Procs = AndProcessVirtual<Sys>;
  using namespace probInvars;
  auto spawnH = Spawn{spawnRateH, Species::H};
  auto jumpH = Jump{sys, Species::H, emJumpH, wForAll, lenJumpH};
  auto destroyH = Destroy{sys.props().temperature(), Species::H, emDestroyH, wForAll};
  auto destroyH2 = Destroy{sys.props().temperature(), Species::H2, emDestroyH2, wForAll};
  ProcsH procsH;
  procsH.add(std::move(jumpH));
  procsH.add(std::move(destroyH));
  ProcsH2 procsH2;
  procsH2.add(destroyH2);
  SystemProcessWrap<Spawn> s{spawnH};
  Procs procs;
  //procsH.isDecide(false);
  procsH2.isDecide(false);
  s.isDecide(false);
  procs.add(std::move(procsH));
  procs.add(std::move(procsH2));
  procs.add(s);
  return procs;
}

template <class Particle, class Neigh> auto reactions(Neigh &neigh) {
  using Add = Reflexive<RecombineH2<Particle>, Neigh &, Dist, Species, Particle>;
  using Reaction = Fallback<PBC, Add>;
  using ReactionSp = SpeciesSelect<Reaction, NoReact, Species>;
  auto add = Add{RecombineH2<Particle>{}, Species::H, neigh, Dist{}};
  return ReactionSp{Reaction{PBC{}, add}, NoReact{}, Species::H};
}

template <class Props, class Particle, class Reaction>
auto dirty(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, Particle>;
  auto dirt = Dirt{reaction};
  dirt.blackList(Species::H2);
  return dirt;
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

auto run(double temperature, size_t maxSteps, size_t dumpInterval, double cellLen) {
  using Particle = PbcParticle<SimpleParticle<typename Props::Coords>>;
  auto props = sysProps(temperature);
  //auto neigh = neighCell(props, cellLen);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions<Particle>(neigh);
  //auto dirt = dirty<Props, Particle>(reacts);
  auto hook = dirty<Props, Particle>(reacts);
  //auto hook = both(dirt, neigh);
  auto sys = system<Particle>(props, hook);
  auto procs =  processes<Particle>(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps/10, dumpInterval/10);
  using std::chrono::system_clock;
  auto start = system_clock::now();
  auto sysTime = pokmc(procs, pred, sys, 5.0);
  std::chrono::duration<double> elapsed = system_clock::now() - start;
  std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << pred.steps() << std::endl;
  return elapsed.count();
}

auto timeAll(double cellLen) {
  constexpr auto maxSteps = probInvars::maxSteps;
  constexpr auto dumpSteps = 500 * 100;
  for (auto temperature : {1200, 900, 600, 300}) {
    std::cout << "\n==============" << temperature << "=================\n";
    auto t = run(temperature, maxSteps, dumpSteps, cellLen);
    std::cout << "\n exec time: ====== " << t << " ========\n";
  }
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

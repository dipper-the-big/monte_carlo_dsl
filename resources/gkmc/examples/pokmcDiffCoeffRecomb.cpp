#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <okmc/particles/SystemSingle.hpp>
#include <okmc/SystemProps.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/hooks/Dummy.hpp>

#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <pokmc/processcollections/OrProcess.hpp>
#include <pokmc/processcollections/AndProcess.hpp>

#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/SpawnSingle.hpp>

#include <pokmc.hpp>
#include <helper/RandReal.hpp>
#include <helper/dist.hpp>

#include <okmc/loggers/Logger.hpp>
#include <okmc/predicates/benchPreds.hpp>

#include <okmc/reactions/Fallback.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/reactions/pbc.hpp>
#include <okmc/reactions/Absorb.hpp>


#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/Eager.hpp>
#include <okmc/hooks/NeighCell.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/hooks/Both.hpp>

using namespace gkmc;

enum class Species : int {
  He,
  HeC
};

using Props = SystemProps<3, Species>;
using P = SizeParticle<PbcParticle<SimpleParticle<Props::Coords>>>;

namespace probInvars {
constexpr double spawnRate = 22;
constexpr double jumpEm = 1.257;
constexpr double jumpW = 6.4e12;
constexpr double jumpLen = 9.0;
constexpr double emitDist = 10.0;
constexpr double minDist = 5.0;
constexpr double boxDim = 1e4;
constexpr double temperature = 800;
constexpr auto nParticles = 500;
constexpr auto maxSteps = 50000;
constexpr auto dumpInterval = 5000;
constexpr auto maxTime = 20.0;
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

template <class Sys> struct Logger {
  auto operator () (Sys &p, double t, size_t steps) {
    if ((int)t <= lastSecPrinted) return;
    using std::chrono::system_clock;
    std::chrono::duration<double> elapsed = system_clock::now() - start;
    lastSecPrinted = (int)t;
    pCounts[1] = p.template particles<P>(Species::He).size();
    auto it = p.template begin<P>(Species::HeC);
    auto last = p.template end<P>(Species::HeC);
    for (;it != last; ++it) pCounts[it->size()]++;
    for (const auto &it : pCounts)
      std::cout << elapsed.count() << '\t' << t << '\t' << steps << '\t'
                << it.first << '\t' << it.second << '\n';
    std::cout << std::endl;
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


auto sysProps() {
  auto p = Props{};
  p.temperature(probInvars::temperature);
  p.boxDimensions(Props::Coords{
      {probInvars::boxDim, probInvars::boxDim, probInvars::boxDim}});
  return p;
}

template <class Neigh, class Cap> auto reactions(Neigh &neigh, Cap &cap) {
  using AbsorbT = All<AbsorbNew<P, Species>, Neigh &, Cap &, P>;
  using SpeciesPBC = SpeciesSelect<PBC, NoReact, Species>;
  using Reaction = Fallback<SpeciesPBC, AbsorbT>;
  return Reaction{
      SpeciesPBC{PBC{}, NoReact{}, Species::He},
      AbsorbT{AbsorbNew<P, Species>{Species::HeC, 2}, neigh, cap}};
}

template <class Reaction>
auto dirty(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, P>;
  auto dirt = Dirt{reaction};
  return dirt;
}

auto neighBasic() { return NeighBasic{true}; }

template <class Props> auto neighCell(const Props &p, double cellLen) {
  return NeighCell<Props>{p.boxDimensions(), cellLen, true};
}

template <class Uno, class Dos> auto both(Uno &uno, Dos &dos) {
  return Both<Uno &, Dos &>{uno, dos};
}

template <class Sys>
auto fillSys(Sys& sys) {
  auto spawn = SpawnSingle<P, Species>{0.0, Species::He};
  for (auto i = 0; i < probInvars::nParticles; ++i) {
    spawn.execute(sys);
  }
}


template <class Hook>
auto system(Props &props, Hook &hook) {
  using Sys = SystemSingle<Props, Hook &, P>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<P>(Species::He, "He");
  sys.template addSpecies<P>(Species::HeC, "HeC");
  fillSys(sys);
  return sys;
}

template <class Sys>
auto processes(Sys& sys) {
  using Proc = Jump3DRandom<Sys, P>;
  using ProcWrap = OrProcess<Proc>;
  using Procs = AndProcess<ProcWrap>;
  using probInvars::jumpEm; using probInvars::jumpW; 
  using probInvars::jumpLen; using probInvars::temperature;
  auto proc1 = Proc{Species::He, jumpEm, jumpW, jumpLen, temperature};
  auto procWrap1 = ProcWrap{};
  procWrap1.add(proc1);
  auto procs = Procs{};
  procs.add(std::move(procWrap1));
  return procs;
}

void testDiffCoeff(double w) {
  auto props = sysProps();
  auto neigh = neighCell(props, 100.0);
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist); 
  auto dirt = dirty(reacts);
  auto hook = both(dirt, neigh);
  auto sys = system(props, hook);
  auto procs = processes(sys);
  auto oo = Logger<decltype(sys)>();
  auto pred = TimePredicate<decltype(sys), decltype(oo)>{probInvars::maxTime, oo};
  auto sysTime = pokmc(procs, pred, sys, w);
  std::cout << "Finished simulation time: " << sysTime << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr<<"pass temperature as first argument.\n";
    return 1;
  }
  auto w = std::stod(argv[1]);
  testDiffCoeff(w);
  return 0;
}

/*

  System sys;
  for (auto i = 0; i < 100; ++i) sys.add("W", 1, System::Coords{{0.0, 0.0, 0.0}});
  ProcessCollection pc;
  //pc.add(sp, std::make_unique<Jump3DBcc>(0.3954875, 1.7e12, 1100, 3.157));  //3.5
  pc.add(sp, std::make_unique<JumpRODBccProcess>(0.35454, 2.97e12, 1100, 3.160));  //3.5
  pc.add(sp, std::make_unique<Jump1DBccProcess>(0.29504, 2.97e12, 1100, 3.160));  //3.5
  bkl(pc, RandReal<double>(0.0, 1.0), Predicate(100000), Logger(1), sys);
  //bkl(pc, RandReal<double>(0.0, 1.0), Predicate(100000), Logger(1000), sys);

*/

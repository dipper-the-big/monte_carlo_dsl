#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>

#include <okmc/particles/SystemSingle.hpp>
#include <okmc/SystemProps.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/hooks/Dummy.hpp>

#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processcollections/TupleProcessCollection.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <pokmc/processcollections/OrProcess.hpp>
#include <pokmc/processcollections/AndProcess.hpp>

#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/SpawnSingle.hpp>

#include <bkl.hpp>
#include <helper/RandReal.hpp>
#include <helper/dist.hpp>

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
constexpr double jumpW = 5.390264868e13;//6.0e12;
constexpr double jumpLen = 3.0;
constexpr double emitDist = 10.0;
constexpr double minDist = 5.0;
constexpr double boxDim = 1e4;
constexpr double temperature = 800;
constexpr auto nParticles = 100;
constexpr auto maxSteps = 50000;
constexpr auto dumpInterval = 5000;
constexpr auto maxTime = 2000000.0;
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
private:
  std::ofstream _fout;
public:
  Logger(std::string fname) : _fout{fname, std::ios::out} {}
  auto operator () (Sys &p, double t, size_t steps) {
    if ((int)t <= lastSecPrinted) return;
    using std::chrono::system_clock;
    std::chrono::duration<double> elapsed = system_clock::now() - start;
    lastSecPrinted = (int)t;
    pCounts[1] = p.template particles<P>(Species::He).size();
    auto it = p.template begin<P>(Species::HeC);
    auto last = p.template end<P>(Species::HeC);
    for (;it != last; ++it) pCounts[it->size()]++;
    for (const auto &it : pCounts) {
      _fout << elapsed.count() << '\t' << t << '\t' << steps << '\t'
                << it.first << '\t' << it.second << '\n';
    }
    _fout << std::endl;
    pCounts.clear();
    // printSystemTuple(p);
    start = system_clock::now();
  }
  ~Logger() {
    _fout.close();
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
  // fillSys(sys);
  return sys;
}

template <class Sys>
auto processes(Sys& sys) {
  using Spawn = SpawnSingle<P, Species>;
  using Jump = Jump3DRandom<Sys, P>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap>;
  using namespace probInvars;
  auto spawn = Spawn{spawnRate, Species::He};
  auto jumpWrap =
      JumpWrap{Jump{Species::He, jumpEm, jumpW, jumpLen, temperature}};
  return Procs{spawn, jumpWrap};
}

void testDiffCoeff(double w, std::string fname) {
  auto props = sysProps();
  auto neigh = neighCell(props, 100.0);
  // auto neigh = neighBasic();
  auto dist = Dist{};
  auto reacts = reactions(neigh, dist); 
  auto dirt = dirty(reacts);
  auto hook = both(dirt, neigh);
  auto sys = system(props, hook);
  auto procs = processes(sys);
  Logger<decltype(sys)> oo(fname);
  auto pred = TimePredicate<decltype(sys), decltype(oo)&>{probInvars::maxTime, oo};
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  auto rnd = RandReal<>{0.0, 1.0};
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  if (status == bklStatus::NoProcesses) { std::cout<<"Error in rate calculation.\n"; }
  std::cout<<"Finished simulation time: "<<sysTime<<std::endl;
}

int main(int argc, char* argv[]) {
  assert(argc > 0);
  MPI_Init(NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  using std::to_string; using std::string;
  //for (auto w : {50.0, 10.0, 5.0, 3.0, 2.0, 1.0, 0.1, 0.05, 0.01}) {
  for (auto w : {1.0}) {
    string fn = string(argv[1]) + string("_") +
                to_string(world_rank) + string("_") + to_string(w);
    testDiffCoeff(w, fn);
  }
  MPI_Finalize();
  return 0;
}

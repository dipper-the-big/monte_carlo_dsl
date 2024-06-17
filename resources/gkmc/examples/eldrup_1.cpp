#include <chrono>
#include <iostream>
#include <random>
#include <string>

#include <bkl.hpp>
#include <helper/RandReal.hpp>
#include <okmc/loggers/Logger.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/Eager.hpp>
#include <okmc/hooks/NeighCell.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/hooks/Both.hpp>
#include <okmc/processcollections/TupleProcessCollection.hpp>
#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <okmc/processes/Jump3DBcc.hpp>
#include <okmc/processes/SpawnSingle.hpp>
#include <okmc/processes/VirtualProcess.hpp>
#include <okmc/reactions/Fallback.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/reactions/pbc.hpp>
#include <okmc/reactions/Absorb.hpp>
#include <okmc/reactions/Annihilation.hpp>
#include <okmc/particles/SystemTuple.hpp>
#include <okmc/particles/SystemSingle.hpp>
#include <okmc/particles/SystemPtr.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/predicates/benchPreds.hpp>
#include <okmc/SystemProps.hpp>

enum class Sp : int { I, IC, V, VC };

using namespace gkmc;

using Props = SystemProps<3, Sp>;

using Particle = PbcParticle<SimpleParticle<Props::Coords>>;

namespace probInvars {
constexpr auto temperature = 600.0;
constexpr auto latticeConstant = 2.87;
constexpr auto a0 = latticeConstant;
constexpr auto boxDim = 1e4;
constexpr auto biasI = 1.15;
constexpr auto biasV = 1.00;
const auto r0 = std::sqrt(3) * a0 / 4.0;
constexpr auto aCube = a0 * a0 * a0;
constexpr auto epsilon = 1e-6;
constexpr auto v0 = 5 * 1e12;
constexpr auto emJumpV = 0.53;
constexpr auto emJumpI = 0.30;
const auto jumpLen = std::sqrt(3) * a0 / 2.0;
const auto spawnRate = 49.0;
const auto nParticlesInitial = 5000;
}

auto sysProps() {
  using probInvars::boxDim; using probInvars::temperature;
  using probInvars::latticeConstant;
  auto p = Props{};
  p.temperature(temperature);
  p.latticeConstant(latticeConstant);
  p.boxDimensions(std::array<double, 3>{boxDim, boxDim, boxDim});
  return p;
}

auto radiusC(Particle p, double bias = 1.15) {
  using invars::pi;
  using probInvars::aCube; using probInvars::r0; using probInvars::epsilon;
  auto x = std::cbrt((3.0 * aCube ) / (4.0 * pi * 2.0));
  return bias * (r0 + epsilon + x * std::cbrt(p.size()) - x);
}

struct Dist {
  using P = Particle;
  template <class Sys>
  double operator()(Sys sys, Particle p, Particle n, typename Sys::Pid pid,
                    typename Sys::Pid nid) {
    auto sp1 = sys.template getSpecies<P>(pid);
    auto sp2 = sys.template getSpecies<P>(nid);
    return _calc(p, n, sp1) + _calc(n, p, sp2);
  }
private:
  double _calc(Particle& p, Particle& n, Sp& sp) {
    using probInvars::biasI; using probInvars::biasV;
    auto bias = (sp == Sp::I || sp == Sp::IC) ? biasI : biasV;
    return radiusC(p, bias);
  }
};

auto neighBasic() { return NeighBasic{true}; }

template <class Props> auto neighCell(const Props &p, double cellLen) {
  return NeighCell<Props>{p.boxDimensions(), cellLen, true};
}

template <class Neigh> auto reactions(Neigh &neigh) {
  using ClusterI = Reflexive<AbsorbEq<Sp>, Neigh &, Dist, Sp, Particle>;
  using ClusterV = Reflexive<AbsorbEq<Sp>, Neigh &, Dist, Sp, Particle>;
  using ClusterIC = Any<Absorb, Neigh &, Dist, Sp, Particle>;
  using ClusterVC = Any<Absorb, Neigh &, Dist, Sp, Particle>;
  auto fn = [](auto &sys, auto &pid, auto &nid, Particle &p, Particle &n) {
    auto newSize = abs(p.size() - n.size());
    if (newSize != 1) return false;
    auto& coords = (p.size() > n.size()) ? p.coords() : n.coords();
    auto& bigSp = (p.size() > n.size()) ? sys.template getSpecies<Particle>(pid) 
                                        : sys.template getSpecies<Particle>(nid);
    auto newSp = (bigSp == Sp::IC) ? Sp::I : Sp::V;
    sys.template add<Particle>(newSp, Particle{coords, newSize});
    return true;
  };
  using Anni = AnnihilationFn<decltype(fn)>;
  using AnnihilateIV = Symmetric<Anni, Neigh &, Dist, Sp, Particle, Particle>;
  return fallBack(
      PBC{}
      , ClusterI{AbsorbEq<Sp>{Sp::IC}, Sp::I, neigh, Dist{}}
      , ClusterV{AbsorbEq<Sp>{Sp::VC}, Sp::V, neigh, Dist{}}
      , ClusterIC{Absorb{}, std::vector<Sp>{Sp::IC, Sp::I}, neigh, Dist{}}
      , ClusterVC{Absorb{}, std::vector<Sp>{Sp::VC, Sp::V}, neigh, Dist{}}
      , AnnihilateIV{Anni{fn},
                     std::vector<Sp>{Sp::I, Sp::IC},
                     std::vector<Sp>{Sp::V, Sp::VC},
                     neigh,
                     Dist{}}
  );
}

template <class Sys> auto processes(Sys &sys) {
  using probInvars::emJumpI; using probInvars::v0; using probInvars::emJumpV; 
  using probInvars::spawnRate; using probInvars::nParticlesInitial;
  using Jump = Jump3DBcc<Sys, Particle>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  auto jumpI = JumpWrap{Jump{sys, Sp::I, emJumpI, v0}};
  auto jumpV = JumpWrap{Jump{sys, Sp::V, emJumpV, v0}};
  using Jumps = VectorProcessCollection<Sys, JumpWrap>;
  Jumps js;
  js.add(jumpI);
  js.add(jumpV);
  using Spawn = SpawnSingleOfEach<Sys, Particle>;
  auto spawn = Spawn{spawnRate, std::vector<Sp>{Sp::I, Sp::V}};
  // filling the system with some particles
  for (auto i = 0; i < nParticlesInitial; ++i) {
    spawn.execute(sys);
  }
  if (nParticlesInitial > 0) sys.update(0.0,0.0);
  using Procs = TupleProcessCollection<Sys, Jumps, Spawn>;
  return Procs{js, spawn};
}

template <class Props, class Reaction>
auto dirty(Reaction &reaction) {
  using Dirt = DirtyList<Props, Reaction&, Particle>;
  auto dirt = Dirt{reaction};
  return dirt;
}

template <class Uno, class Dos> auto both(Uno &uno, Dos &dos) {
  return Both<Uno &, Dos &>{uno, dos};
}

template <class Props, class Hook>
auto system(Props &props, Hook &hook) {
  using Sys = SystemTuple<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Sp::I, "I");
  sys.template addSpecies<Particle>(Sp::IC, "IC");
  sys.template addSpecies<Particle>(Sp::V, "V");
  sys.template addSpecies<Particle>(Sp::VC, "VC");
  return sys;
}

auto run(double maxTime, size_t dumpInterval, double cellLen) {
  auto props = sysProps();
  auto neigh = neighCell(props, cellLen);
  // NeighBasic neigh = neighBasic();
  auto reacts = reactions(neigh);
  //auto dirt = dirty<Props, Particle>(reacts);
  auto dirt = dirty<Props>(reacts);
  auto hook = both(dirt, neigh);
  auto sys = system(props, hook);
  auto procs =  processes(sys);
  auto oo = LogCounts<decltype(sys)>(dumpInterval);
  auto pred = TimePredicate<decltype(sys), decltype(oo)>{maxTime, oo};
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

auto timeAll(double cellLen) {
  constexpr auto maxTime = 60.0;
  constexpr size_t dumpSteps = 1000000;
  run(maxTime, dumpSteps, cellLen);
}

int main(int argc, char *argv[]) {
  auto cellLen = 100.0;
  timeAll(cellLen);
  return 0;
}

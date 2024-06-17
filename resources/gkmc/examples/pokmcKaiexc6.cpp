#include <array>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <pokmc.hpp>
#include <helper/RandReal.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <pokmc/processcollections/OrProcess.hpp>
#include <pokmc/processcollections/AndProcess.hpp>
#include <pokmc/processes/ParticleProcessCounter.hpp>
#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/ProcessCounter.hpp>
#include <okmc/reactions/Annihilation.hpp>
#include <okmc/reactions/ReactionInactiveCounter.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/particles/SystemTuple.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/SystemProps.hpp>
#include <helper/invars.hpp>

using coords = std::array<double, 3>;

template <class Gen> auto givePoint(Gen &gen) {
  auto _dist = gen();
  RandReal<> _rnd{0.0, 1.0};

  auto th = 2 * gkmc::invars::pi * _rnd();
  auto phi = acos(1. - 2. * _rnd());
  return coords{{_dist * std::sin(phi) * std::cos(th),
                 _dist * std::sin(phi) * std::sin(th),
                 _dist * std::cos(phi)}};
}

struct Result {
  std::size_t steps;
  double time;
  std::size_t SystemTuple;
  std::size_t iJumps;
  std::size_t vJumps;
  double ratio;
};

enum class Species : int { I, V };

using namespace gkmc;
auto run(double temperature, double w) {
  using std::array;
  using std::normal_distribution;
  using std::size_t;
  using std::vector;

  // problem constants
  constexpr double ws[] = {1.717e15, 0.001282e15};
  constexpr double ems[] = {1.37, 0.1};
  constexpr int counts[] = {100, 100};
  constexpr double initRadius[] = {60.0, 20.0};
  const auto rCalc = [](const auto &, const auto &, const auto &, const auto &, const auto &) { return 4.0; };
  constexpr double jump = 10.0;//2.35;

  using Props = SystemProps<3, Species>;
  using Particle = SimpleParticle<typename Props::Coords>;

  using Reaction = Reflexive<Annihilation, NeighBasic, decltype(rCalc), Species, Particle>;
  using Hook = DirtyList<Props, Reaction, Particle>;
  using Sys = SystemTuple<Props, Hook &, Particle>;
  using Proc = Jump3DRandom<Sys, Particle>;
  using ProcCounter = ParticleProcessCounter<Proc>;
  using ProcWrap = OrProcess<ProcCounter>;
  using Procs = AndProcess<ProcWrap>;

  auto props = Props{};
  props.temperature(temperature);
  auto reaction = Reaction{Annihilation{}, Species::I, Species::V, NeighBasic{false}, rCalc};
  auto hook = Hook{reaction};
  auto sys = Sys{props, hook};
  // adding process with wrapper for counting
  ProcWrap procI;
  procI.add(ProcCounter{Proc{Species::I, ems[0], ws[0], jump, temperature}});
  ProcWrap procV;
  procV.add(ProcCounter{Proc{Species::V, ems[1], ws[1], jump, temperature}});
  Procs procand;
  procand.add(procI);
  procand.add(procV);
  // adding particles
  for (auto _ = 0; _ <= counts[0]; ++_) {
    RandReal<normal_distribution<>> rnd{0,
                                        initRadius[0]}; // TODO: outside
    auto particle = Particle{givePoint(rnd)};
    sys.add(Species::I, particle);
  }
  for (auto _ = 0; _ <= counts[0]; ++_) {
    RandReal<normal_distribution<>> rnd{-initRadius[1],
                                        initRadius[1]}; // TODO: outside
    auto particle = Particle{givePoint(rnd)};
    sys.add(Species::V, particle);
  }
  sys.update(0.0, 0.0);

  size_t steps = 0;
  constexpr auto dumpInterval = 5000000;
  auto printSystemTupleCount = [&steps, dumpInterval](Sys &p, double time) {
    if (!steps || steps % dumpInterval)
      return;
    auto ps = p.template particles<Particle>(Species::I); // TODO: for V
    for (auto it : ps) {
      std::cout << time << "\t ";
      std::for_each(begin(it.coords()), end(it.coords()),
                    [](const double &c) { std::cout << c << "\t "; });
      std::cout << '\n';
    }
  };
  auto pred = [&](Sys &sys, double time) {
    //printSystemTuple(sys, time);
    ++steps;
    return time < 1e-7;
  };

  auto sysTime = 0.0;
  sysTime = pokmc(procand, pred, sys, w);
  std::cout << "Finished simulation time: " << sysTime << std::endl;
  return Result{steps, sysTime,
                sys.template particles<Particle>(Species::I).size(),
                procand.procs()[0].procs()[0].count(), procand.procs()[1].procs()[0].count(),
                (double)procand.procs()[0].procs()[0].count() / procand.procs()[1].procs()[0].count()};
  /*
  std::tie(sysTime, status) = bkl(procs, rnd, pred, SystemTuple);
  return Result{steps, sysTime,
                SystemTuple.template particles<Particle>(Species::I).size(),
                procs.procs()[0].count(), procs.procs()[1].count()};
                */
}

int main(int argc, char *argv[]) {
  // run for different temperatures
  /*
  if (argc < 2) {
    std::cerr << "Please provide factor (w) as argument.\n";
    return 1;
  }
  */
  constexpr auto nSamples = 10;
  using sample = std::array<Result, nSamples>;
  using temperature = int;
  std::map<temperature, sample> results;
  for (auto it : {0.01,0.05,0.1,0.5,1.0,2.0,3.0,5.0,7.0,10.0}) {
    for (auto i = 0; i < nSamples; ++i) {
      results[int(it * 1000)][i] = run(500, it);
    }
  }
  std::cout << "temperature\tsteps\ttime\tnInterstitial\tiJumps\tvJumps\n";
  for (auto it : results) {
    for (auto jt : it.second) {
      std::cout << it.first / 1000.0 << '\t' << jt.steps << '\t' << jt.time << '\t'
                << jt.SystemTuple << '\t' << jt.iJumps << '\t' << jt.vJumps
                << '\t' << jt.ratio << '\n';
    }
  }
  return 0;
}

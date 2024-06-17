#include <array>
#include <iostream>
#include <random>
#include <string>
#include <vector>

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

using coords = std::array<double, 3>;

template <class Gen> auto givePoint(Gen &gen) {
  coords c;
  std::generate(begin(c), end(c), gen);
  return c;
}

struct Result {
  std::size_t steps;
  double time;
  std::size_t SystemTuple;
  std::size_t iJumps;
  std::size_t vJumps;
};

enum class Species : int { I, V };

using namespace gkmc;
auto run(double temperature) {
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
  constexpr double jump = 2.35;
  constexpr size_t inactiveLimit = 1e6;

  using Props = SystemProps<3, Species>;
  using Particle = SimpleParticle<typename Props::Coords>;

  using Reaction = Reflexive<Annihilation, NeighBasic, decltype(rCalc), Species, Particle>;
  using ReactionCounter = ReactionInactiveCounter<Reaction &>;
  using Hook = DirtyList<Props, ReactionCounter &, Particle>;
  using Sys = SystemTuple<Props, Hook &, Particle>;
  using Proc = Jump3DRandom<Sys, Particle>;
  using ProcWrap = nParticleProcessWrap<Proc>;
  using ProcCounter = ProcessCounter<ProcWrap>;
  using Procs = VectorProcessCollection<Sys, ProcCounter>;

  auto props = Props{};
  props.temperature(temperature);
  size_t curInactiveCount = 0;
  auto reaction = Reaction{Annihilation{}, Species::I, Species::V, NeighBasic{false}, rCalc};
  auto reactionCounter = ReactionCounter{reaction, curInactiveCount};
  auto hook = Hook{reactionCounter};
  auto sys = Sys{props, hook};
  // adding process with wrapper for counting
  Procs procs;
  procs.add(ProcCounter{ProcWrap{Proc{Species::I, ems[0], ws[0], jump, temperature}}});
  procs.add(ProcCounter{ProcWrap{Proc{Species::V, ems[1], ws[1], jump, temperature}}});
  // adding particles
  for (auto _ = 0; _ <= counts[0]; ++_) {
    RandReal<normal_distribution<>> rnd{-initRadius[0],
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
    return inactiveLimit > ++curInactiveCount;
  };

  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  auto rnd = RandReal<>{0.0, 1.0};
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  }
  std::cout << "Finished simulation time: " << sysTime << std::endl;
  return Result{steps, sysTime,
                sys.template particles<Particle>(Species::I).size(),
                procs.procs()[0].count(), procs.procs()[1].count()};
}

int main(int argc, char *argv[]) {
  // run for different temperatures
/*
  constexpr auto nSamples = 2;
  using sample = std::array<Result, nSamples>;
  using temperature = int;
  std::map<temperature, sample> results;
  for (auto it : {500, 1500, 2500}) {
    for (auto i = 0; i < nSamples; ++i) {
      results[it][i] = run(it);
    }
  }
  std::cout << "temperature\tsteps\ttime\tnInterstitial\tiJumps\tvJumps\n";
  for (auto it : results) {
    for (auto jt : it.second) {
      std::cout << it.first << '\t' << jt.steps << '\t' << jt.time << '\t'
                << jt.SystemTuple << '\t' << jt.iJumps << '\t' << jt.vJumps
                << '\n';
    }
  }
  */
  return 0;
}

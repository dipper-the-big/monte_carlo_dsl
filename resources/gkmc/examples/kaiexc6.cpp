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
  std::size_t particles;
  std::size_t iJumps;
  std::size_t vJumps;
  double ratio;
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
  constexpr double jump = 10.0;//2.35;

  using Props = SystemProps<3, Species>;
  using Particle = SimpleParticle<typename Props::Coords>;

  using Reaction = Reflexive<Annihilation, NeighBasic, decltype(rCalc), Species, Particle>;
  using Hook = DirtyList<Props, Reaction, Particle>;
  using Sys = SystemTuple<Props, Hook &, Particle>;
  using Proc = Jump3DRandom<Sys, Particle>;
  using ProcWrap = nParticleProcessWrap<Proc>;
  using ProcCounter = ProcessCounter<ProcWrap>;
  using Procs = VectorProcessCollection<Sys, ProcCounter>;

  auto props = Props{};
  props.temperature(temperature);
  auto reaction = Reaction{Annihilation{}, std::vector<Species>{Species::I}, std::vector<Species>{Species::V}, NeighBasic{false}, rCalc};
  auto hook = Hook{reaction};
  auto sys = Sys{props, hook};
  // adding process with wrapper for counting
  Procs procs;
  procs.add(ProcCounter{ProcWrap{Proc{Species::I, ems[0], ws[0], jump, temperature}}});
  procs.add(ProcCounter{ProcWrap{Proc{Species::V, ems[1], ws[1], jump, temperature}}});
  // adding SystemTuple
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
  auto status = bklStatus::NoProcesses;
  auto rnd = RandReal<>{0.0, 1.0};
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  }
  std::cout << "Finished simulation time: " << sysTime << std::endl;
  return Result{steps, sysTime,
                sys.template particles<Particle>(Species::I).size(),
                procs.procs()[0].count(), procs.procs()[1].count(),
                (double)procs.procs()[0].count()/ procs.procs()[1].count()};
}

int main(int argc, char *argv[]) {
  // run for different temperatures
  constexpr auto nSamples = 10;
  using sample = std::array<Result, nSamples>;
  using temperature = int;
  std::map<temperature, sample> results;
  for (auto it : {500}/*{500, 1500, 2500}*/) {
    for (auto i = 0; i < nSamples; ++i) {
      results[it][i] = run(it);
    }
  }
  std::cout << "temperature\tsteps\ttime\tnInterstitial\tiJumps\tvJumps\n";
  for (auto it : results) {
    for (auto jt : it.second) {
      std::cout << it.first << '\t' << jt.steps << '\t' << jt.time << '\t'
                << jt.particles << '\t' << jt.iJumps << '\t' << jt.vJumps
                << '\t' << jt.ratio << '\n';
    }
  }
  return 0;
}

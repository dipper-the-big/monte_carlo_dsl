// HeMCD.cpp

#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <bkl.hpp>
#include <helper/RandReal.hpp>
#include <okmc/hooks/NeighBasic.hpp>
#include <okmc/hooks/NeighCell.hpp>
#include <okmc/hooks/DirtyList.hpp>
#include <okmc/hooks/Dummy.hpp>
#include <okmc/hooks/Both.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>
#include <okmc/processcollections/TupleProcessCollection.hpp>
#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/SpawnSingle.hpp>
#include <okmc/processes/Emission.hpp>
#include <okmc/reactions/pbc.hpp>
#include <okmc/reactions/Fallback.hpp>
#include <okmc/reactions/Binary.hpp>
#include <okmc/particles/SystemSingle.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/SystemProps.hpp>

using namespace gkmc;

using Species = int;
using Props = SystemProps<3, Species>;
using P = PbcParticle<SimpleParticle<Props::Coords>>;

struct AbsorbHe {
  template <class Sys, class Particle>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, Particle &p, Particle &n) {
    sys.template remove<Particle>(pid);
    sys.template remove<Particle>(nid);
    sys.template add<Particle>(pid.first + nid.first, Particle{p.coords()});
    return true;
  }
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

template <class Sz> class HeEmission : public Dummy {
private:
public:
  using Pid = typename Props::Pid;
  HeEmission(double temperature, double dist, Sz sz)
      : _dist{dist}, _calcRate{HeEmissionRate{temperature}}, _rates{0.0, 4}, _sizeCalc{sz} {}

  template <class Sys>
  double rate(Sys &sys) {
    return _rates.last();
  }

  template <class Sys>
  const std::vector<double> &cumulativeRates(Sys &sys) {
    return _rates.buffer();
  }

  template <class Sys>
  auto executeSelected (size_t sp, Sys &sys, double dice) {
    auto prevCumulativeRate = 0.0;
    if (sp != 0) prevCumulativeRate = _rates.buffer()[sp - 1];
    auto subDice = dice - prevCumulativeRate;
    auto temp = (subDice / _calcRate(sp));
    auto particleIndex = (size_t)(temp);
    return emit(sp, particleIndex, sys, temp - particleIndex);
  }

  template <class Sys>
  void emit (Species sp, const size_t &pIndex, Sys &sys,
                       const double subDice) {
    // select particle using subDice.
    auto pid = sys.cookPid(sp, pIndex);
    auto& particle = sys.template particle<P>(pid);
    auto away = _sizeCalc(pid.first) + _dist; // TODO
    auto coords = rndOnSphere(subDice, _rnd(), away);

    std::transform(std::begin(coords), std::end(coords),
                   std::begin(particle.coords()), std::begin(coords),
                   std::plus<double>());

    sys.template remove<P>(pid);
    sys.template add<P>(pid.first - 1, P{particle.coords()});
    sys.template add<P>(1, P{coords});
  }

  template <class Sys>
  auto execute(Sys &sys, const double &dice) {
    auto pickedCumulative =
        std::lower_bound(begin(_rates.buffer()), end(_rates.buffer()), dice);
    auto indexSpecies = std::distance(begin(_rates.buffer()), pickedCumulative);
    return executeSelected(indexSpecies, sys, dice);
  }

  void _ratesIncrement(size_t index) {
    _rates.incrementAt(index, _calcRate(index), 0.0);
  }

  void _ratesDecrement(size_t index) {
   assert(_rates.buffer().size() > index);
    _rates.decrementAt(index, _calcRate(index), 0.0);
  }

  template <class Sys>
  void added(const Sys &, P &p, const Pid &pid) {
    if (pid.first > 1) _ratesIncrement(pid.first);
  }
  template <class Sys>
  void removed(const Sys &, P &p, const Pid &pid) {
    if (pid.first > 1) _ratesDecrement(pid.first);
  }
private:
  double _dist;
  std::vector<double> _ratesCurrent;
  SizeBuffer<HeEmissionRate> _calcRate;
  BufferEdit<double> _rates;
  RandReal<> _rnd{0.0, 1.0};
  Sz _sizeCalc;
};

auto run(double temperature) {
  using std::array;
  using std::normal_distribution;
  using std::size_t;
  using std::vector;
  // problem constants
  constexpr auto spawnRate = 22;
  constexpr auto jumpEm = 1.257;
  constexpr auto jumpW = 6.3e12;
  constexpr auto jumpLen = 9.0;
  constexpr auto emitDist = 10.0;
  constexpr auto maxTime = 360.0;
  constexpr auto minDist = 5.0;
  constexpr auto boxDim = 1e4;
  auto radius = [](std::size_t s) { return pow(8.37 * pow(s, 1.02), 0.33333); };
  auto sizeBuffer = SizeBuffer<decltype(radius)>{radius, 80};
  auto cap = [&sizeBuffer, &minDist](auto &s, auto &p, auto &n, auto &pid, auto nid) {
    return sizeBuffer(p.size()) + sizeBuffer(n.size()) + minDist;
  };
  //using Neigh = NeighBasic;
  using Neigh = NeighCell<Props>;
  using Recomb = All<AbsorbHe, Neigh&, decltype(cap) &, P>;
  using TypePBC = SpeciesSelect<PBC, NoReact, Species>;

  using Reaction = Fallback<TypePBC, Recomb>;

  using Emit = HeEmission<decltype(sizeBuffer) &>;

  using Dirt = DirtyList<Props, Reaction &, P>;
  //using Hook = Both<Emit &, Dirt>;
  using Hook1 = Both<Emit &, Dirt>;
  using Hook = Both<Both<Emit &, Dirt>, Neigh&>;
  using Sys = SystemSingle<Props, Hook &, P>;

  using Spawn = SpawnSingle<P, Species>;
  using Jump = Jump3DRandom<Sys, P>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  //using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap>;
  using Procs = TupleProcessCollection<Sys, Spawn, JumpWrap, Emit&>;
  
  auto props = Props{};
  props.temperature(temperature);
  props.boxDimensions(Props::Coords{{boxDim, boxDim, boxDim}});

  //auto neigh = Neigh{true};
  auto neigh = Neigh{props.boxDimensions(), 100.0, true};
  auto reaction = Reaction{TypePBC{PBC{}, NoReact{}, Species{1}},
                           Recomb{AbsorbHe{}, neigh, cap}};
  auto emit = Emit{temperature, emitDist, sizeBuffer};
  auto hook = Hook{Hook1{emit, Dirt{reaction}}, neigh};
  auto ps = Sys{props, hook};
  auto spawn = Spawn{spawnRate, Species{1}};
  auto jumpWrap =
      JumpWrap{Jump{Species{1}, jumpEm, jumpW, jumpLen, temperature}};
  // Procs procs{spawn, jumpWrap};
  Procs procs{spawn, jumpWrap, emit};

  size_t steps = 0;
  using std::chrono::system_clock;
  auto start = system_clock::now();
  auto lastSecPrinted = 0;
  auto printSystemCount =
      [&steps, &start, &lastSecPrinted](Sys &p, double t) {
        if ((int)t <= lastSecPrinted) return;
        std::chrono::duration<double> elapsed = system_clock::now() - start;
        lastSecPrinted = (int)t;
        for (int i = 1; p.template isSpecies<P>(i); ++i) {
          std::cout << elapsed.count() << '\t' << t << '\t' << steps << '\t'
                    << i << '\t' << p.template particles<P>(i).size() << '\n';
        }
        std::cout<<std::endl;
        // printSystemTuple(p);
        start = system_clock::now();
      };
  auto pred = [&](Sys &sys, double time) {
    printSystemCount(sys, time);
    ++steps;
    return maxTime > time;
  };
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  auto rnd = RandReal<>{0.0, 1.0};

  std::tie(sysTime, status) = bkl(procs, rnd, pred, ps);
  if (status == bklStatus::NoProcesses) {
    std::cout << "Error in rate calculation.\n";
  } else {
    std::cout << "Finished simulation.\ttime: " << sysTime
              << "\tsteps: " << steps << std::endl;
    //printSystemTuple(SystemTuple);
  }
}

int main(int argc, char *argv[]) {
  run(800);
  return 0;
}
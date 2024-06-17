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
#include <okmc/processes/Jump3DBcc.hpp>
#include <okmc/processes/Emission.hpp>
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

enum class Sp : int { I, IC, V, VCS, S};

using namespace gkmc;

using Props = SystemProps<3, Sp>;

struct Particle : public SimpleParticle<Props::Coords> {
  Particle(const Coords &c, int size = 1, int nSoulte = 0)
      : SimpleParticle{c}, _size{size}, _nSolute{nSoulte}, _prevCoords{c} {}
  void size(int val) { _size = val; }
  int size() const { return _size; }
  void nSolute(int val) { _nSolute = val; }
  int nSolute() const { return _nSolute; }

  const auto &coords() const { return SimpleParticle::coords(); }
  const auto &prevCoords() const { return _prevCoords; }
  void coords(const Coords &c) {
    _prevCoords = SimpleParticle::coords();
    SimpleParticle::coords(c);
  }
  void coords(value_type val, size_t i) {
    _prevCoords[i] = SimpleParticle::coords()[i];
    SimpleParticle::coords(val, i);
  }
private:
  int _size;
  int _nSolute;
  Coords _prevCoords;
};

namespace probInvars {
constexpr auto temperature = 600.0;
constexpr auto latticeConstant = 2.87;
constexpr auto a0 = latticeConstant;
constexpr auto boxDim = a0 * 100.0;
// 20keV, Fe-0.2%Cu, absorbing boundaries, set B
// set B according to atomistic simulations, all are mobile
// attempt frequency, v0 = 6 * 1e12 / s
constexpr auto biasI = 1.15;
constexpr auto biasV = 1.00;
const auto r0 = std::sqrt(3) * a0 / 4.0;
constexpr auto aCube = a0 * a0 * a0;
constexpr auto epsilon = 1e-6;
constexpr auto v0 = 6 * 1e12;
constexpr auto emJumpV = 0.69;
constexpr auto emJumpI = 0.04;
const auto jumpLen = std::sqrt(3) * a0 / 2.0;
constexpr auto emitEbV = 0.2;
constexpr auto emitEforV = 1.6;
constexpr auto emitEbVS = 0.1;
constexpr auto emitEforVS = 1.88;
constexpr auto emitEbI = 1.0;
constexpr auto emitEforI = 4.0;
}

auto sysProps() {
  auto p = Props{};
  p.temperature(probInvars::temperature);
  p.box(100, probInvars::latticeConstant);
  return p;
}

auto radiusC(Particle p, double bias = 1.15) {
  using invars::pi;
  using probInvars::aCube; using probInvars::r0; using probInvars::epsilon;
  auto x = std::cbrt((3.0 * aCube ) / (4.0 * pi * 2.0));
  return bias * (r0 + epsilon + x * std::cbrt(p.size() + p.nSolute()) - x);
}

auto nSVCS(Particle p) {
  using invars::pi;
  using probInvars::a0; using probInvars::aCube; using probInvars::epsilon;
  auto x = std::cbrt((3.0 * aCube ) / (4.0 * pi * 2.0));
  return (a0 / 2.0 + epsilon + x * std::cbrt(p.nSolute()) - x);
}

auto nSpS(Particle p) {
  using invars::pi;
  using probInvars::aCube; using probInvars::r0; using probInvars::epsilon;
  auto x = std::cbrt((3.0 * aCube ) / (4.0 * pi * 2.0));
  return (r0 / 2.0 + epsilon + x * std::cbrt(p.nSolute()) - x);
}

struct Dist {
  using P = Particle;
  template <class Sys>
  double operator()(Sys sys, Particle p, Particle n, typename Sys::Pid pid,
                    typename Sys::Pid nid) {
    auto sp1 = sys.template getSpecies<P>(pid);
    auto sp2 = sys.template getSpecies<P>(nid);
    return calc(p, n, sp1, sp2) + calc(n, p, sp1, sp2);
  }
private:
  double calc(Particle& p, Particle& n, Sp& sp1, Sp& sp2) {
    if (sp1 == Sp::S || (sp1 == Sp::VCS && p.size() == 0)) {
      if (sp2 == Sp::S || (sp2 == Sp::VCS && n.size() == 0)) {
        return nSpS(p);
      } else {
        return nSVCS(p);
      }
    } else {
      using probInvars::biasI; using probInvars::biasV;
      auto bias = (sp1 == Sp::I || sp1 == Sp::IC) ? biasI : biasV;
      return radiusC(p, bias);
    }
  }
};

auto neighBasic() { return NeighBasic{true}; }

template <class Props> auto neighCell(const Props &p, double cellLen) {
  return NeighCell<Props>{p.boxDimensions(), cellLen, true};
}

struct AnniCustom {
  using P = Particle;
  template <class Sys>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, P &n) {
    // std::cout << "============ AnniCustom ============\n";
    if (p.size() == n.size()) {
      sys.template remove<P>(pid);
      sys.template remove<P>(nid);
      auto nSolute = abs(p.nSolute() - n.nSolute());
      if (nSolute != 0) {
        sys.template add<P>(Sp::S, P{p.coords(), 0, nSolute});
      }
    } else {
      auto sp = sys.template getSpecies<P>(pid);
      if (p.size() > n.size()) {
        if (sp == Sp::I) {
          addInterstitial(sys, pid, nid, p, n);
        } else {
          addVacancy(sys, pid, nid, p, n);
        }
      } else {
        if (sp == Sp::I) {
          addVacancy(sys, nid, pid, n, p);
        } else {
          addInterstitial(sys, nid, pid, n, p);
        }
      }
    }
    return true;
  }
private:
  template <class Sys>
  auto addInterstitial(Sys &sys, const typename Sys::Pid &pid,
                       const typename Sys::Pid &nid, P &p, P &n) {
    sys.template remove<P>(nid);
    auto sz = p.size() - n.size();
    if (sz > 1) p.size(p.size() - n.size());
    else {
      sys.template remove<P>(pid);
      sys.template add<P>(Sp::I, P{p.coords(), sz});
    }
    if (n.nSolute() != 0) {
      sys.template add<P>(Sp::S, P{n.coords(), 0, n.nSolute()});
    }
  }
  template <class Sys>
  auto addVacancy(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, P &n) {
    sys.template remove<P>(nid);
    auto sz = p.size() - n.size();
    if (sz > 1 || p.nSolute() > 0) p.size(p.size() - n.size());
    else {
      sys.template remove<P>(pid);
      sys.template add<P>(Sp::V, P{p.coords(), sz});
    }
  }
};

class AbsorbS {
public:
  template <class Sys, class P, class N>
  auto operator()(Sys &sys, const typename Sys::Pid &pid,
                  const typename Sys::Pid &nid, P &p, N &n) {
    if (p.size() > n.size()) {
      p.size(p.size() + n.size());
      p.nSolute(p.nSolute() + n.nSolute());
      sys.template edited<P>(pid);
      sys.template remove<N>(nid);
    } else {
      n.size(n.size() + p.size());
      n.nSolute(n.nSolute() + p.nSolute());
      sys.template edited<N>(nid);
      sys.template remove<P>(pid);
    }
    return true;
  }
};

template <class Neigh> auto reactions(Neigh &neigh) {
  using ClusterI = Reflexive<AbsorbEq<Sp>, Neigh &, Dist, Sp, Particle>;
  using ClusterV = Reflexive<AbsorbEq<Sp>, Neigh &, Dist, Sp, Particle>;
  using ClusterIC = Any<Absorb, Neigh &, Dist, Sp, Particle>;
  using ClusterVC = Any<AbsorbS, Neigh &, Dist, Sp, Particle>;
  using AnnihilateIV = Symmetric<AnniCustom, Neigh &, Dist, Sp, Particle, Particle>;
  return fallBack(
      ABC{}
      , ClusterI{AbsorbEq<Sp>{Sp::IC}, Sp::I, neigh, Dist{}}
      , ClusterV{AbsorbEq<Sp>{Sp::VCS}, Sp::V, neigh, Dist{}}
      , ClusterIC{Absorb{}, std::vector<Sp>{Sp::IC, Sp::I}, neigh, Dist{}}
      , ClusterVC{AbsorbS{}, std::vector<Sp>{Sp::VCS, Sp::V, Sp::S}, neigh, Dist{}}
      , AnnihilateIV{AnniCustom{},
                     std::vector<Sp>{Sp::I, Sp::IC},
                     std::vector<Sp>{Sp::V, Sp::VCS},
                     neigh,
                     Dist{}}
  );
  //return NoReact{};
}

template <class Sys> class EmissionVS {
public:
  using P = Particle;
  using System = Sys;
  using Species = Sp;
  using ParticleType = Particle;
  using rate_type = double;

private:
  double _dist = 2.0;
  double _rate;
  RandReal<> _rnd{0.0, 1.0};
  Sp _spOn = Sp::VCS;
  Sp _spEmit = Sp::S;
  Sp _spShift = Sp::V;
  int _sizeShift = 1;

public:
  EmissionVS(rate_type rate) : _rate{rate} {}

  auto executeSelected(const size_t &pIndex, System &sys,
                       const rate_type subDice) {
    auto pid = sys.cookPid(_spOn, pIndex);
    auto& particle = sys.template particle<P>(pid);
    auto away = radiusC(particle) + _dist;
    //  http://corysimon.github.io/articles/uniformdistn-on-sphere/
    auto th = 2 * invars::pi * _rnd(); // subDice; TODO
    auto phi = acos(1. - 2. * _rnd());
    typename System::Coords coords{{away * std::sin(phi) * std::cos(th),
                                    away * std::sin(phi) * std::sin(th),
                                    away * std::cos(phi)}};
    std::transform(std::begin(coords), std::end(coords),
                   std::begin(particle.coords()), std::begin(coords),
                   std::plus<rate_type>());
    if (particle.size() == 1 && particle.nSolute() == 1) {
      return;
    } else if (particle.size() == 1 && particle.nSolute() > 1) {
      sys.template remove<P>(pid);
      sys.template add<P>(Sp::S, P{particle.coords(), 0, particle.nSolute() - 1});
    } else if (particle.size() == 2 && particle.nSolute() == 1) {
      sys.template remove<P>(pid);
      sys.template add<P>(Sp::V, P{particle.coords(), 1, 0});
    } else {
      particle.size(particle.size() - 1);
      particle.nSolute(particle.nSolute() - 1);
      sys.template edited<P>(pid);
    }
    sys.template add<P>(_spEmit, P{coords, 1, 1});
  }

  const rate_type &ratePerParticle(const System &s) const { return _rate; }

  const Sp &species() const { return _spOn; }
};

struct RateEmit {
  RateEmit(double emB, double emFor) : _emB{emB}, _emFor{emFor} {}
  auto operator () (Particle p) {
    using probInvars::v0; using probInvars::temperature;
    auto x = cbrt(p.size() * p.size()) - cbrt((p.size() - 1.0) * (p.size() - 1.0));
    x /= (cbrt(4.0) - 1.0);
    x *= (_emB - _emFor);
    x += _emFor;
    return v0 * exp(-x / (temperature * invars::kB));
  }
private:
  double _emB;
  double _emFor;
};

struct RateEmitMod {
  RateEmitMod() : _r{probInvars::emitEbVS, probInvars::emitEforVS} {}
  auto operator () (Particle p) {
    if (p.nSolute() == 0) { return 0.0; }
    return _r(p);
  }
private:
  RateEmit _r;
};

template <class Sys> auto processes(Sys &sys) {
  using probInvars::emJumpI; using probInvars::v0; using probInvars::emJumpV; 
  using probInvars::emitEbV; using probInvars::emitEforV;
  using probInvars::emitEbI; using probInvars::emitEforI;
  using Jump = Jump3DBcc<Sys, Particle>;
  using JumpWrap = nParticleProcessWrap<Jump>;
  auto jumpI = JumpWrap{Jump{sys, Sp::I, emJumpI, v0}};
  auto jumpV = JumpWrap{Jump{sys, Sp::V, emJumpV, v0}};
  constexpr auto pinv = 1./100.; 
  constexpr auto qinv = 1./1000.; 
  constexpr auto s = 0.51; 
  auto v0VCS = [pinv, qinv](const Particle& p) {
    auto res = 1.0; // v0;
    if (p.size() > 2) res *= pow(pinv, p.size() - 2);
    if (p.nSolute() > 1) res *= pow(qinv, p.nSolute() - 1);
    return res;
  };
  auto v0IC = [s](Particle p) {
    return pow(p.size(), -s);
  };
  using JumpVCS = ParticleRerateWrapper<Jump, decltype(v0VCS)>;
  auto jumpVCS = JumpVCS{Jump{sys, Sp::VCS, emJumpV, v0}, v0VCS};
  using Jump1D = Jump1DBcc<Sys, Particle>;
  using JumpIC = ParticleRerateWrapper<Jump1D, decltype(v0IC)>;
  auto jumpIC = JumpIC{Jump1D{sys, Sp::IC, emJumpI, v0}, v0IC};
  using Emit = Emission<Sys, Particle, Particle, decltype(radiusC)*>;
  using EmitC = ParticleRerateWrapper<Emit, RateEmit>;
  auto emitV = EmitC{Emit{radiusC, 2.0, 1.0, Sp::VCS, Sp::V, Sp::V},
                     RateEmit{emitEbV, emitEforV}};
  auto emitI = EmitC{Emit{radiusC, 2.0, 1.0, Sp::IC, Sp::I, Sp::I},
                     RateEmit{emitEbI, emitEforI}};
  using EmitVS = ParticleRerateWrapper<EmissionVS<Sys>, RateEmitMod>;
  auto emitVS = EmitVS{EmissionVS<Sys>{1.0}, RateEmitMod{}};
  /*
  std::cout << jumpIC.rate(sys) << " rate of jumpIC\n";
  std::cout << emitI.rate(sys) << " rate of emitI\n";
  std::cout << jumpV.rate(sys) << " rate of jumpV\n";
  std::cout << jumpVCS.rate(sys) << " rate of jumpVCS\n";
  std::cout << emitV.rate(sys) << " rate of emitV\n";
  std::cout << emitVS.rate(sys) << " rate of emitVS\n";
  */
  using JumpSingles = VectorProcessCollection<Sys, JumpWrap>;
  JumpSingles js;
  js.add(jumpI);
  js.add(jumpV);
  // using EmitSingles = VectorProcessCollection<Sys, EmitC>;
  //EmitSingles es;
  // es.add(emitI);
  // es.add(emitV);
  using Procs = TupleProcessCollection<Sys, JumpSingles, JumpVCS, JumpIC>;
  return Procs{js, jumpVCS, jumpIC};
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

template <class Gen> auto givePoint(Gen &gen) {
  auto _dist = gen();
  auto mid = probInvars::boxDim / 2.0;
  RandReal<> _rnd{0.0, 1.0};
  auto th = 2 * gkmc::invars::pi * _rnd();
  auto phi = acos(1. - 2. * _rnd());
  return Props::Coords{{mid + _dist * std::sin(phi) * std::cos(th),
                 mid + _dist * std::sin(phi) * std::sin(th),
                 mid + _dist * std::cos(phi)}};
}

template <class Props, class Hook>
auto system(Props &props, Hook &hook) {
  using Sys = SystemTuple<Props, Hook &, Particle>;
  auto sys = Sys{props, hook};
  sys.template addSpecies<Particle>(Sp::I, "I");
  sys.template addSpecies<Particle>(Sp::IC, "IC");
  sys.template addSpecies<Particle>(Sp::V, "V");
  sys.template addSpecies<Particle>(Sp::VCS, "VCS");
  sys.template addSpecies<Particle>(Sp::S, "S");
  return sys;
}

template <class Sys>
auto fill(Sys& sys) {
  for (auto _ = 0; _ < 24; ++_) {
    RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 1, 0};
    sys.add(Sp::I, particle);
  }
  for (auto _ = 0; _ < 3; ++_) {
    RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 2, 0};
    sys.add(Sp::IC, particle);
  }
  for (auto _ = 0; _ < 2; ++_) {
    RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 3, 0};
    sys.add(Sp::IC, particle);
  }
  for (auto _ = 0; _ < 1; ++_) {
    RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 4, 0};
    sys.add(Sp::IC, particle);
  }
  for (auto _ = 0; _ < 1; ++_) {
    RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 6, 0};
    sys.add(Sp::IC, particle);
  }
  for (auto _ = 0; _ < 1; ++_) {
    RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 10, 0};
    sys.add(Sp::IC, particle);
  }
  // Vacancy
  for (auto _ = 0; _ < 36; ++_) {
    RandReal<std::normal_distribution<>> rnd{40.0, 40.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 1, 0};
    sys.add(Sp::V, particle);
  }
  for (auto _ = 0; _ < 5; ++_) {
    RandReal<std::normal_distribution<>> rnd{40.0, 40.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 2, 1};
    sys.add(Sp::VCS, particle);
  }
  for (auto _ = 0; _ < 1; ++_) {
    RandReal<std::normal_distribution<>> rnd{40.0, 40.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 4, 0};
    sys.add(Sp::VCS, particle);
  }
  for (auto _ = 0; _ < 1; ++_) {
    RandReal<std::normal_distribution<>> rnd{40.0, 40.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 6, 2};
    sys.add(Sp::VCS, particle);
  }
  for (auto _ = 0; _ < 5; ++_) {
    RandReal<std::normal_distribution<>> rnd{60.0, 40.0}; // TODO: outside
    auto particle = Particle{givePoint(rnd), 0, 1};
    sys.add(Sp::S, particle);
  }
  sys.update(0.0, 0.0);
}

template <class Sys>
auto printSys (Sys& sys, double time) { 
  auto SystemTuple = sys.template particles<Particle>(Sp::IC);
  for (auto it : SystemTuple) {
    std::cout<<time<<"\t ";
    std::for_each(begin(it.coords()), end(it.coords()), [](const double& c) {
      std::cout<<c<<"\t ";
    });
    std::cout<<'\n';
  }
};

template <class Sys, class Particle> struct Logger {
  Logger(size_t dumpInterval, bool isPrintPos)
    : _dumpInterval{dumpInterval}, _isPrintPos{isPrintPos} {}
  auto operator()(Sys &p, double time, size_t steps) {
    using std::chrono::system_clock;
    if (!steps || steps % _dumpInterval) return;
    if (_isPrintPos) { printSys(p, time); return; }
    std::chrono::duration<double> elapsed = system_clock::now() - start;
    auto nI = p.template particles<Particle>(Sp::I).size();
    auto& IC = p.template particles<Particle>(Sp::IC);
    auto nIC = 0;
    for (auto &x : IC) {
      nIC += x.size();
    }
    auto nV = p.template particles<Particle>(Sp::V).size();
    auto& VCS = p.template particles<Particle>(Sp::VCS);
    auto nVCS = 0;
    for (auto &x : VCS) {
      nVCS += x.size();
    }
    auto nS = p.template particles<Particle>(Sp::S).size();
    std::cout << elapsed.count() << '\t' << time << '\t' << steps << '\t'
         << nI << '\t'<< nIC << '\t'  << nV << '\t' << nVCS << '\t' << nS << '\n';
    start = system_clock::now();
  }
private:
  std::chrono::system_clock::time_point start =
      std::chrono::system_clock::now();
  size_t _dumpInterval;
  bool _isPrintPos;
};

template <class Sys, class Particle>
auto predicate(size_t maxSteps, size_t dumpInterval, bool isPrintPos = false) {
  using Log = Logger<Sys, Particle>;
  return StepsPredicate<Sys, Log>{maxSteps, Log{dumpInterval, isPrintPos}};
}

auto run(size_t maxSteps, size_t dumpInterval, double cellLen) {
  auto props = sysProps();
  //auto neigh = neighCell(props, cellLen);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions(neigh);
  //auto dirt = dirty<Props, Particle>(reacts);
  auto hook = dirty<Props>(reacts);
  //auto hook = both(dirt, neigh);
  auto sys = system(props, hook);
  fill(sys);
  auto procs =  processes(sys);
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

template <class Sys>
auto fillTest1DJump(Sys& sys) {
  RandReal<std::normal_distribution<>> rnd{80.0, 30.0}; // TODO: outside
  auto particle = Particle{givePoint(rnd), 3, 0};
  sys.add(Sp::IC, particle);
  sys.add(Sp::V, Particle{givePoint(rnd), 1, 0});
  sys.add(Sp::VCS, Particle{givePoint(rnd), 3, 1});
  sys.update(0.0, 0.0);
}

auto test1Djump(size_t maxSteps, size_t dumpInterval, double cellLen) {
  maxSteps = 100000;
  dumpInterval = 10;
  auto props = sysProps();
  //auto neigh = neighCell(props, cellLen);
  NeighBasic neigh = neighBasic();
  auto reacts = reactions(neigh);
  //auto dirt = dirty<Props, Particle>(reacts);
  auto hook = dirty<Props>(reacts);
  //auto hook = both(dirt, neigh);
  auto sys = system(props, hook);
  fillTest1DJump(sys);
  auto procs =  processes(sys);
  auto pred = predicate<decltype(sys), Particle>(maxSteps, dumpInterval, true);
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
  constexpr auto maxSteps = 1000000;
  constexpr auto dumpSteps = 1;
  //test1Djump(maxSteps, dumpSteps, cellLen);
  run(maxSteps, dumpSteps, cellLen);
}

int main(int argc, char *argv[]) {
  auto cellLen = 50.0;
  timeAll(cellLen);
  return 0;
}

#include <iostream>
#include <string>
#include <vector>

#include <okmc/particles/SystemTuple.hpp>
#include <okmc/SystemProps.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/hooks/Dummy.hpp>

#include <okmc/processcollections/VectorProcessCollection.hpp>
#include <okmc/processcollections/nParticleProcessWrap.hpp>

#include <okmc/processes/Dummy.hpp>

#include <okmc/processes/Jump3DRandom.hpp>
#include <okmc/processes/Jump3DBcc.hpp>

#include <bkl.hpp>
#include <helper/RandReal.hpp>

using namespace gkmc;

enum class Species : int {
  only
};

void testDiffCoeff(double temperature) {
  using Props = SystemProps<3, Species>;
  using Particle = gkmc::SimpleParticle<typename Props::Coords>;
  using Hook = Dummy;
  using Sys = SystemTuple<Props, Hook, Particle>;
  using Proc = Jump3DRandom<Sys, Particle>;
  //using Proc = Jump3DBcc<Sys, Particle>;
  using ProcWrap = nParticleProcessWrap<Proc>;
  using Procs = VectorProcessCollection<Sys, ProcWrap>;
  auto props = Props{};
  props.temperature(temperature);
  props.latticeConstant(2.855);
  auto sys = Sys{props, Hook{}};
  auto sp = Species::only;
  //auto proc1 = Proc{sys, sp, 0.21599278, 3.4124463302909754e13/*6e12*/};
  auto em = 1.257;// 0.21599278 
  auto w = 6e12;//5.390264868e13;//3.4124463302909754e13;// 6.0e12;//6.4e12;// 3.4124463302909754e13
  auto len = 9.0;// 2.472
  auto proc1 = Proc{sp, em, w, len, temperature};
  auto procWrap1 = ProcWrap{proc1};
  auto procs = Procs{};
  procs.add(std::move(procWrap1));
  auto rnd = RandReal<>{0.0, 1.0};
  auto steps = 0;
  auto sysTime = 0.0;
  auto status = bklStatus::NoProcesses;
  constexpr auto maxSteps = 2500000;
  constexpr auto dumpInterval = 5000;
  constexpr auto nSystemTuple = 500;
  auto pos = typename Sys::Coords{{0., 0., 0.}};
  for (auto i = 0; i < nSystemTuple; ++i) {
    sys.add(sp, Particle{pos});
  }
  auto printSys = [&](Sys& sys, double time) { 
    if (!steps || steps%dumpInterval) return;
    auto SystemTuple = sys.template particles<Particle>(sp);
    for (auto it : SystemTuple) {
      std::cout<<time<<"\t ";
      std::for_each(begin(it.coords()), end(it.coords()), [](const double& c) {
        std::cout<<c<<"\t ";
      });
      std::cout<<'\n';
    }
  };
  auto pred = [&](Sys& sys, double time) { 
    printSys(sys, time); 
    return (++steps) < maxSteps;  
  };
  std::tie(sysTime, status) = bkl(procs, rnd, pred, sys);
  if (status == bklStatus::NoProcesses) { std::cout<<"Error in rate calculation.\n"; }
  std::cout<<"Finished simulation time: "<<sysTime<<std::endl;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr<<"pass temperature as first argument.\n";
    return 1;
  }
  auto temperature = std::stod(argv[1]);
  testDiffCoeff(temperature);
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

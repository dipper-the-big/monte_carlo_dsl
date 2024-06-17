#include <iostream>
#include <string>
#include <vector>

#include <mcd.hpp>
#include <mcd/diffuse.hpp>
#include <okmc/particles/SystemSingle.hpp>
#include <okmc/SystemProps.hpp>
#include <okmc/particles/SimpleParticle.hpp>
#include <okmc/hooks/Dummy.hpp>

#include <helper/RandReal.hpp>
#include <helper/invars.hpp>

using namespace gkmc;

enum class Species : int { only };

using Props = SystemProps<3, Species>;
using Particle = gkmc::SimpleParticle<typename Props::Coords>;

void testDiffCoeff(double temperature) {
  using Hook = Dummy;
  using Sys = SystemSingle<Props, Hook, Particle>;
  auto props = Props{};
  props.temperature(temperature);
  props.latticeConstant(2.855);
  auto sys = Sys {props, Hook{}};
  auto sp = Species::only;

  auto steps = 0;
  constexpr auto maxSteps = 5000;
  constexpr auto dumpInterval = 10;
  constexpr auto nParticles = 500;
  auto pos = typename Sys::Coords{{0., 0., 0.}};
  for (auto i = 0; i < nParticles; ++i) {
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
  constexpr auto delTFromPokmc = 1.5e-6;//1.26334e-5;//1.30e-5;// 1.35e-5;// 1.26334e-5;
  constexpr auto w = 1.00;
  constexpr auto jumpLen = 3.0;
  auto sysTime = mcd(Diffuse3D<Particle, Species>{sp, jumpLen / w}, pred, sys, delTFromPokmc / (w * w));
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
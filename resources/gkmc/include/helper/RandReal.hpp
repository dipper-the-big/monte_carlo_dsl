// RandReal.hpp

#ifndef GKMC_RANDREAL_H
#define GKMC_RANDREAL_H

#include <algorithm>
#include <functional>
#include <random>

#include <helper/invars.hpp>

// http://codereview.stackexchange.com/questions/109260/seed-stdmt19937-from-stdrandom-device
template <class T = std::mt19937, std::size_t N = T::state_size>
auto ProperlySeededRandomEngine() -> typename std::enable_if<!!N, T>::type {
  typename T::result_type random_data[N];
  std::random_device source;
  std::generate(std::begin(random_data), std::end(random_data),
                std::ref(source));
  std::seed_seq seeds(std::begin(random_data), std::end(random_data));
  T seededEngine(seeds);
  return seededEngine;
}

template <class U = std::uniform_real_distribution<>> struct RandReal {
public:
  using T = typename U::result_type;
  RandReal(T min, T max, std::seed_seq customSeed) : dis{min, max} {
    gen.seed(customSeed);
  }

  RandReal(T min, T max) : dis{min, max} { gen = ProperlySeededRandomEngine(); }

  T operator()() { return dis(gen); }

private:
  std::mt19937 gen;
  U dis;
  //std::uniform_real_distribution<T> dis;
};

//  http://corysimon.github.io/articles/uniformdistn-on-sphere/
auto rndOnSphere(double rnd1, double rnd2, double radius) {
    auto th = 2 * gkmc::invars::pi * rnd2; // subDice; TODO
    auto phi = acos(1. - 2. * rnd1);
    return std::array<double, 3> {{radius * std::sin(phi) * std::cos(th),
                                    radius * std::sin(phi) * std::sin(th),
                                    radius * std::cos(phi)}};
}

#endif

// Logger.hpp

#ifndef GKMC_LOGGER_H
#define GKMC_LOGGER_H

#include <iostream>
#include <tuple>

namespace gkmc {

/*!
 * Logs on standard output according to particle print
 * */
 template <class System>
class LogParticles {
public:
  using Tup = typename System::TypesTuple;

  LogParticles (std::size_t interval) : _n{0}, _interval{interval} { }

  template <std::size_t I>
  bool printEm(System& sys) {
    using T = std::tuple_element_t<I, Tup>;
    const auto& particles = sys.template particles<T>();
    for (const auto& it : particles) {
      std::cout << sys.template speciesName<T>(it.first) << ":\n";
      size_t i = 0;
      for (const auto& jt : it.second) {
        std::cout << i++ << '\t';
        printParticle(jt);
        std::cout << '\n';
      }
      std::cout << '\n';
    }
    return false;
  }

  template<std::size_t... Is>
  void printEmHelper(System& sys, const std::index_sequence<Is...>&) {
    bool x[] = {printEm<Is>(sys)...};
  }

  void operator() (System& sys, double time, size_t steps) {
    if (steps == 1 || steps % _interval) return;
    std::cout << "\nIteration: " << steps << '\t' << "Time: " << time << '\n';
    _time = time;
    printEmHelper(sys, std::make_index_sequence<std::tuple_size<Tup>::value>{});
  }

  void operator() (System&, double time, std::string msg) {
    std::cout<<"Info- Iteration: "<<(_n)<<'\t'<<"Time: "<<time<<'\t'<<msg<<'\n';
  }

  const auto& steps() const { return _n; }
  const auto& time() const { return _time; }

private:
  size_t _n;
  double _time;
  size_t _interval;
};

/*!
 * Logs particle counts on standard output for each species
 * */
 template <class System>
class LogCounts {
public:
  using Tup = typename System::TypesTuple;

  LogCounts (std::size_t interval) : _n{0}, _interval{interval} { }

  template <std::size_t I>
  bool printEm(System& sys) {
    using T = std::tuple_element_t<I, Tup>;
    const auto& particles = sys.template particles<T>();
    for (const auto& it : particles) {
      std::cout << sys.template speciesName<T>(it.first) << ": " 
                << it.second.size() << '\n';
    }
    return false;
  }

  template<std::size_t... Is>
  void printEmHelper(System& sys, const std::index_sequence<Is...>&) {
    bool x[] = {printEm<Is>(sys)...};
  }

  void operator() (System& sys, double time, size_t steps) {
    if (steps == 1 || steps % _interval) return;
    std::cout << "\nIteration: " << steps << '\t' << "Time: " << time << '\n';
    _time = time;
    printEmHelper(sys, std::make_index_sequence<std::tuple_size<Tup>::value>{});
  }

  void operator() (System&, double time, std::string msg) {
    std::cout<<"Info- Iteration: "<<(_n)<<'\t'<<"Time: "<<time<<'\t'<<msg<<'\n';
  }

  const auto& steps() const { return _n; }
  const auto& time() const { return _time; }

private:
  size_t _n;
  double _time;
  size_t _interval;
};

}

#endif
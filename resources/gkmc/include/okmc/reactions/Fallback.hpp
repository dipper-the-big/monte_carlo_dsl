// Fallback.hpp

#ifndef GKMC_Fallback_H
#define GKMC_Fallback_H

#include <iostream>
#include <tuple>
#include <type_traits>

namespace gkmc {

/** 
 * Tries the next reaction unless one of them succeedes.
 * Returns false if none does.
 * */
template <class... Rs> class Fallback {
public:
  Fallback(Rs... rs) : _rs{std::forward<Rs>(rs)...} {}

  template <int sz, class T, class U>
  bool _helper(T &a, const typename T::Pid &b, U &c, std::false_type) {
    if ((std::get<sz>(_rs))(a, b, c)) return true;
    return _helper<sz + 1>(a, b, c, std::integral_constant<bool, sz + 1 == lastIndex>{});
  }

  template <int sz, class T, class U>
  bool _helper(T &a, const typename T::Pid &b, U &c, std::true_type) {
    return false;
  }

  template <class T, class U>
  auto operator()(T &a, const typename T::Pid &b, U &c) {
    return _helper<0>(a, b, c, std::integral_constant<bool, 0 == lastIndex>{});
  }

private:
  std::tuple<Rs...> _rs;
  static constexpr int lastIndex = sizeof...(Rs);
};

template <class... Rs>
auto fallBack(Rs&&... rs) {
  return Fallback<Rs...>{std::forward<Rs>(rs)...};
}

/** 
 * Tries the first reaction if the type is same as the input particle else 
 * tries second.
 * */
template <class R1, class R2, class Particle = typename R1::Particle> 
class TypeSelect {
public:
  TypeSelect(R1 r1, R2 r2) : _r1{r1}, _r2{r2} {}

  template <class T>
  auto operator()(T &a, const typename T::Pid &b, Particle &c) {
    return _r1(a, b, c);
  }

  template <class T, class U>
  auto operator()(T &a, const typename T::Pid &b, U &c) {
    return _r2(a, b, c);
  }

private:
  R1 _r1;
  R2 _r2;
};

/** 
 * Tries the first reaction if the input predicate succeedes else
 * tries second.
 * */
template <class R1, class R2, class Pred> class PredicateWrap {
public:
  PredicateWrap(R1 r1, R2 r2, Pred p) : _r1{r1}, _r2{r2}, _pred{p} {}

  template <class T, class U>
  auto operator()(T &a, const typename T::Pid &b, U &c) {
    if (_pred(a, b, c)) {
      return _r1(a, b, c);
    }
    return _r2(a, b, c);
  }
private:
  R1 _r1;
  R2 _r2;
  Pred _pred;
};

/** 
 * Tries the first reaction if the input species is same as input else
 * tries second.
 * */
template <class R1, class R2, class Species> class SpeciesSelect {
public:
  SpeciesSelect(R1 r1, R2 r2, Species sp) : _r1{r1}, _r2{r2}, _sp{sp} {}

  template <class T, class U>
  auto operator()(T &a, const typename T::Pid &b, U &c) {
    if (b.first == _sp) {
      return _r1(a, b, c);
    }
    return _r2(a, b, c);
  }
private:
  R1 _r1;
  R2 _r2;
  Species _sp;
};

/**
 * returns false always
 * */
class NoReact {
public:
  template <class T, class U>
  auto operator()(T &a, const typename T::Pid &b, U &c) {
    return false;
  }
};

/*
template <class R1, class R2> class Type2Fallback {
public:
  using Particle = typename R1::Particle;
  using Particle2 = typename R2::Particle;
  Type2Fallback(R1 r1, R2 r2) : _r1{r1}, _r2{r2} {}

  template <class T>
  auto operator()(T &a, const typename T::Pid &b, Particle &c) {
    return _r1(a, b, c);
  }

  template <class T>
  auto operator()(T &a, const typename T::Pid &b, Particle2 &c) {
    return _r1(a, b, c);
  }

  template <class T, class U>
  auto operator()(T &a, const typename T::Pid &b, U &c) {
    return _r2(a, b, c);
  }

private:
  R1 _r1;
  R2 _r2;
};
*/
}
#endif

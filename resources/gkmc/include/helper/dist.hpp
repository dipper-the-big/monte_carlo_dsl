// Dist.hpp

#ifndef GKMC_Dist_H
#define GKMC_Dist_H

#include <array>
#include <cmath>
#include <vector>

#include <helper/invars.hpp>
#include <helper/RandReal.hpp>

template <class Coords> auto calcDist(const Coords &a, const Coords &b) {
  auto dist = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    dist += (a[i] - b[i]) * (a[i] - b[i]);
  }
  return std::sqrt(dist);
}

/*
 * Calculate distance given the periodic boundary conditions
*/
template <class Coords>
auto calcDistMirror(const Coords &a, const Coords &b, const Coords &box) {
  auto dist = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    auto temp = std::abs(a[i] - b[i]);              // non-mirror distance
    temp = std::min(temp, std::abs(box[i] - temp)); // min of non-mirror, mirror
    dist += (temp * temp);
  }
  return std::sqrt(dist);
}

/*
 * Gives a predicate - if two particles are closer than a certain distance
 * with choice to use mirrored or non-mirrored distance calculation.
 */
template <class Coords> class AreClose {
  using T = typename Coords::value_type;

public:
  AreClose() : _isMirror{false}, _box{} {}
  AreClose(const Coords &box) : _isMirror{true}, _box{box} {}

  bool operator()(const Coords &a, const Coords &b, const T &threshold) const {
    return _isMirror ? (calcDistMirror(a, b, _box) < threshold)
                     : (calcDist(a, b) < threshold);
  }

private:
  bool _isMirror;
  Coords _box;
};

/*
 * capture radius calculation for clusters as given in LakiMoca paper for FeCu.
 * link: https://www.sciencedirect.com/science/article/pii/S002231151000231X 
 */
template <class T> // requires T double/float/... num
class CaptureRadiusLakimoca {
public:
  CaptureRadiusLakimoca(T latConst, T biasi = 1.0, T biasj = 1.0)
      : _bias{{biasi, biasj}}, _a0{latConst} {}

  T operator()(size_t sizei, size_t sizej) {
    using gkmc::invars::pi;
    using std::cbrt;
    auto capi =
        (_r0 + cbrt(3. * sizei / pi) - cbrt(3. / pi)) * _a0 * _bias[0] / 2.;
    auto capj =
        (_r0 + cbrt(3. * sizej / pi) - cbrt(3. / pi)) * _a0 * _bias[1] / 2.;
    return capi + capj;
  }

private:
  const std::array<T, 2> _bias;
  const T _a0;
  constexpr static T _r0 = 1.732 / 2. + .01; // .01 for correction
};

// Capture radius buffer wrapper
// Fn requires to be of type: size_t -> auto
// which will have C++ fn. sign. as auto fn(size_t);
template <class Fn, class DomainFn = int> class SizeBuffer {
  using CoDomainFn = std::result_of_t<Fn(DomainFn)>;
public:
  SizeBuffer(Fn fn, int initialBufferSize = 32) : _fn{fn} {
    _preCalc.reserve(initialBufferSize);
    for (DomainFn i = 0; i < initialBufferSize; ++i) {
      _preCalc[i] = _fn(i);
    }
  }
  const auto &operator()(int sizei) {
    if (sizei >= (int)_preCalc.size()) {
      _preCalc.reserve(sizei + 1);
      for (DomainFn i = _preCalc.size(); i <= sizei; ++i) {
        _preCalc.push_back(_fn(i));
      }
    }
    return _preCalc[sizei];
  }

  template<class Sys, class P>
  const auto &operator()(const Sys &sys, const P &p, const typename Sys::Pid &pid) {
    auto sizei = p.size();
    if (sizei >= (int)_preCalc.size()) {
      _preCalc.reserve(sizei + 1);
      for (DomainFn i = _preCalc.size(); i <= sizei; ++i) {
        _preCalc.push_back(_fn(i));
      }
    }
    return _preCalc[sizei];
  }

  auto &preCalc() { return _preCalc; }

private:
  Fn _fn;
  std::vector<CoDomainFn> _preCalc;
};

// capture radius buffer wrapper
// Fn requires to be of type: size_t -> size_t -> auto
// which will have C++ fn. sign. as auto fn(size_t, size_t);
template <class Fn> class CaptureRadiusBuffer {
  using DomainFn = std::result_of_t<Fn(size_t, size_t)>;

public:
  CaptureRadiusBuffer(Fn fn, size_t initialBufferSize = 32) : _fn{fn} {
    _preCalc.reserve(initialBufferSize);
    for (size_t i = 0; i < initialBufferSize; ++i) {
      _preCalc.emplace_back(initialBufferSize);
      for (size_t j = 0; j < initialBufferSize; ++j) {
        _preCalc[i][j] = _fn(i, j);
      }
    }
  }

  const auto &operator()(size_t sizei, size_t sizej) {
    if (sizei >= _preCalc.size() || sizej >= _preCalc.size()) {
      auto mx = std::max(sizei, sizej);
      for (size_t i = 0; i < _preCalc.size(); ++i) {
        _preCalc[i].reserve(mx + 1);
        for (size_t j = _preCalc.size(); j <= mx; ++j) {
          _preCalc[i].push_back(_fn(i, j));
        }
      }
      _preCalc.reserve(mx + 1);
      for (size_t i = _preCalc.size(); i <= mx; ++i) {
        _preCalc.emplace_back();
        _preCalc[i].reserve(mx + 1);
        for (size_t j = 0; j <= mx; ++j) {
          _preCalc[i].push_back(_fn(i, j));
        }
      }
    }
    return _preCalc[sizei][sizej];
  }

private:
  Fn _fn;
  std::vector<std::vector<DomainFn>> _preCalc;
};

/** 
 * Buffer that allows edit, increment or decrement stored value
 * by input value.
 * Right now specialized for only a single use. Can be easily
 * generalized for use with others
 **/
template <class CoDomain> class BufferEdit {
public:
  BufferEdit(CoDomain val, size_t initialBufferSize) {
    _preCalc.reserve(initialBufferSize);
    for (size_t i = 0; i < initialBufferSize; ++i) {
      _preCalc.push_back(val);
    }
  }

  const auto incrementAt(size_t index, CoDomain val, CoDomain defaultVal) {
    if (_preCalc.size() == 0) { _preCalc.push_back(defaultVal); }
    if (index >= _preCalc.size()) {
      _preCalc.reserve(index + 1);
      for (size_t i = _preCalc.size(); i < index; ++i) {
        _preCalc.push_back(_preCalc[i-1] + defaultVal);
      }
      _preCalc.push_back(val + _preCalc[_preCalc.size() - 1]);
    } else {
      for (size_t i = index; i < _preCalc.size(); ++i) {
        _preCalc[i] += val;
      }
    }
  }

  const auto decrementAt(size_t index, CoDomain val, CoDomain defaultVal) {
    for (size_t i = index; i < _preCalc.size(); ++i) {
      _preCalc[i] -= val;
    }
  }

  auto &buffer() { return _preCalc; }
  auto &last() { return _preCalc[_preCalc.size() - 1]; }

private:
  std::vector<CoDomain> _preCalc;
};

#endif

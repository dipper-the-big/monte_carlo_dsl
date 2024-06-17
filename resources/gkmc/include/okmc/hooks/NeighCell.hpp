// NeighCell.hpp

#ifndef GKMC_NeighCell_H
#define GKMC_NeighCell_H

#include <algorithm>
#include <array>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <string>
#include <vector>

#include <helper/dist.hpp>

namespace gkmc {

/**
 * Cell link list 3D neighbour search
 * */

template <class SysProps> class NeighCell {
  static_assert(SysProps::dims == 3, "NeighCell requires 3D box, for 2D check NeighCell2D");
  static constexpr int dims = 3;
  using Coords = typename SysProps::Coords;
  using value_type = typename Coords::value_type;
  using Species = typename SysProps::Species;
  using Pid = typename SysProps::Pid;

public:
  NeighCell(const Coords &boxDimensions, const value_type &cellLength,
            const bool &isMirror)
      : _boxDimensions(boxDimensions), _cellLength{cellLength},
        _isMirror{isMirror} {
    for (auto i = 0; i < dims; ++i) {
      _nCells[i] = std::ceil(_boxDimensions[i] / _cellLength);
    }
    _cells.resize(_nCells[0]);
    auto n = _nCells;
    std::for_each(begin(_cells), end(_cells), [&n](auto &y) {
      y.resize(n[1]);
      std::for_each(begin(y), end(y), [&n](auto &z) { z.resize(n[2]); });
    });
  }

  template <class Sys, class Particle>
  void added(Sys &, Particle &p, const Pid &pid) {
    std::array<int, dims> pos;
    const auto &c = p.coords();
    const Species &species = pid.first; // TODO
    for (auto i = 0; i < dims; ++i)
      pos[i] = c[i] / _cellLength;
    if (pos[0] < _nCells[0] && pos[1] < _nCells[1] && pos[2] < _nCells[2] &&
        pos[0] >= 0 && pos[1] >= 0 && pos[2] >= 0) {
      _cells[pos[0]][pos[1]][pos[2]][species].insert(pid.second);
    }
  }

  template <class Sys, class Particle>
  void removed(Sys &sys, Particle &p, const Pid &pid) {}

  template <class Sys, class Particle>
  void editedSize(Sys &, Particle &p, const Pid &pid, const int &pre) {}

  template <class Sys, class Particle>
  void editedCoords(Sys &, Particle &p, const Pid &pid, const Coords &pre) {
    std::array<int, dims> posCur;
    std::array<int, dims> posPre;
    bool flag = false;
    const Coords &cur = p.coords();
    const auto &index = pid.second;
    const auto &sp = pid.first;
    for (auto i = 0; i < dims; ++i) {
      posCur[i] = cur[i] / _cellLength;
      posPre[i] = pre[i] / _cellLength;
      if (posPre[i] != posCur[i])
        flag = true;
    }
    if (flag) {
      if (posPre[0] < _nCells[0] && posPre[1] < _nCells[1] &&
          posPre[2] < _nCells[2] && posPre[0] >= 0 && posPre[1] >= 0 &&
          posPre[2] >= 0)
        _cells[posPre[0]][posPre[1]][posPre[2]][sp].erase(index);
      if (posCur[0] < _nCells[0] && posCur[1] < _nCells[1] &&
          posCur[2] < _nCells[2] && posCur[0] >= 0 && posCur[1] >= 0 &&
          posCur[2] >= 0)
        _cells[posCur[0]][posCur[1]][posCur[2]][sp].insert(index);
    }
  }
  template <class Sys, class Particle>
  void edited(Sys &, Particle &p, const Pid &pid) {}

  template <class Sys, class Particle>
  void erased(Sys& sys, Particle &p, const Pid &pid) {
    const Coords &c = p.coords();
    const Species &species = pid.first; // TODO
    const size_t &index = pid.second;   // TODO
    std::array<int, dims> pos;
    for (auto i = 0; i < dims; ++i)
      pos[i] = c[i] / _cellLength;
    if (pos[0] < _nCells[0] && pos[1] < _nCells[1] && pos[2] < _nCells[2] &&
        pos[0] >= 0 && pos[1] >= 0 && pos[2] >= 0) {
      _cells[pos[0]][pos[1]][pos[2]][species].erase(index);
    }
    auto n = sys.template particles<Particle>(species).size();
    if (n && index != (n - 1)) {
      //const auto &x = sys.template particles<Particle>(species)[n - 1].coords();
      const auto &x = sys.template particle<Particle>(sys.cookPid(species, n - 1)).coords();
      for (auto i = 0; i < dims; ++i)
        pos[i] = x[i] / _cellLength;
      if (pos[0] < _nCells[0] && pos[1] < _nCells[1] && pos[2] < _nCells[2] &&
          pos[0] >= 0 && pos[1] >= 0 && pos[2] >= 0) {
        _cells[pos[0]][pos[1]][pos[2]][species].erase(n - 1);
        _cells[pos[0]][pos[1]][pos[2]][species].insert(index);
      }
    }
  }
  template <class Sys> void update(Sys &sys, double, double) {}

  template <class P1, class P2, class Sys, class CapRadiusFn>
  auto neighbour(Sys &sys, const typename Sys::Pid &pid,
                 CapRadiusFn &&capRadius,
                 const typename Sys::Species &sp) const {
    auto pred = (_isMirror) ? AreClose<typename Sys::Coords>{sys.props().boxDimensions()}
                            : AreClose<typename Sys::Coords>{};
    const auto &pivot = sys.template particle<P1>(pid);
    const auto &c = pivot.coords();
    std::array<int, dims> cellIndex;
    auto tol = 0.;
    auto& del = sys.template toDel<P2>(sp);
    for (auto i = 0; i < dims; ++i)
      cellIndex[i] = c[i] / _cellLength;
    for (auto i : indices) {
      auto ci = (cellIndex[0] + i) % _nCells[0];
      if (ci < 0)
        ci += _nCells[0];
      for (auto j : indices) {
        auto cj = (cellIndex[1] + j) % _nCells[1];
        if (cj < 0)
          cj += _nCells[1];
        for (auto k : indices) {
          auto ck = (cellIndex[2] + k) % _nCells[2];
          if (ck < 0)
            ck += _nCells[2];
          for (auto it : _cells[ci][cj][ck]) {
            if (sp != it.first)
              continue;
            for (auto jt : it.second) {
              //const auto &jPart = sys.template particles<P2>(it.first)[jt];
              if (del.find(jt) != std::end(del)) continue;
              auto jid = sys.cookPid(it.first, jt);
              const auto &jPart = sys.template particle<P2>(jid);
              if (pivot.id() == jPart.id()) continue;
              tol = capRadius(sys, pivot, jPart, pid, jid);
              if (pred(c, jPart.coords(), tol))
                return std::make_pair(true, jid);
            }
          }
        }
      }
    }
    return std::make_pair(false, typename Sys::Pid{});
  }

  // TODO: call the other overload from inside the loop
  template <class P1, class P2, class Sys, class CapRadiusFn>
  auto neighbour(Sys &sys, const typename Sys::Pid &pid,
                 CapRadiusFn &&capRadius) const {
    auto pred = (_isMirror) ? AreClose<typename Sys::Coords>{sys.props().boxDimensions()}
                            : AreClose<typename Sys::Coords>{};
    const auto &pivot = sys.template particle<P1>(pid);
    const auto &c = pivot.coords();
    std::array<int, dims> cellIndex;
    auto tol = 0.;
    for (auto i = 0; i < dims; ++i)
      cellIndex[i] = c[i] / _cellLength;
    for (auto i : indices) {
      auto ci = (cellIndex[0] + i) % _nCells[0];
      if (ci < 0)
        ci += _nCells[0];
      for (auto j : indices) {
        auto cj = (cellIndex[1] + j) % _nCells[1];
        if (cj < 0)
          cj += _nCells[1];
        for (auto k : indices) {
          auto ck = (cellIndex[2] + k) % _nCells[2];
          if (ck < 0)
            ck += _nCells[2];
          for (auto it : _cells[ci][cj][ck]) {
            if (!sys.template isSpecies<P2>(it.first)) continue;
            auto& del = sys.template toDel<P2>(it.first);
            for (auto jt : it.second) {
              if (del.find(jt) != std::end(del)) continue;
              auto jid = sys.cookPid(it.first, jt);
              const auto &jPart = sys.template particle<P2>(jid);
              if (pivot.id() == jPart.id()) continue;
              tol = capRadius(sys, pivot, jPart, pid, jid);
              if (pred(c, jPart.coords(), tol))
                return std::make_pair(true, jid);
            }
          }
        }
      }
    }
    return std::make_pair(false, typename Sys::Pid{});
  }

private:
  Coords _boxDimensions;
  value_type _cellLength;
  bool _isMirror;
  using Cells = std::vector<
      std::vector<std::vector<std::map<Species, std::unordered_set<size_t>>>>>;
  Cells _cells;
  std::array<int, dims> _nCells;
  std::array<int, 3> indices {{0, -1, 1}};
};

} // namespace !gkmc

#endif // !GKMC_NEIGHCELLSP_H

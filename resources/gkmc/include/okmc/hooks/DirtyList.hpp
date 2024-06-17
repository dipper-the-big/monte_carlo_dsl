// DirtyList.hpp

#ifndef GKMC_DirtyList_H
#define GKMC_DirtyList_H

#include <array>
#include <assert.h>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <helper/meta.hpp>

namespace gkmc {

namespace detail {
template <class Type1, class Cur, class Last, class Tup, class DirtyList,
          class Sys>
void updateNext(
    DirtyList &d, Sys &sys, bool reset,
    typename std::enable_if<!std::is_same<Cur, Last>::value>::type *dummy = 0) {
  using nextType = std::tuple_element_t<meta::typeIndex<Cur, Tup>() + 1, Tup>;
  (reset) ? d.template updateHelper<Type1>(sys)
          : d.template updateHelper<nextType>(sys);
}

template <class Type1, class Cur, class Last, class Tuple, class DirtyList,
          class Sys>
void updateNext(
    DirtyList &d, Sys &sys, bool reset,
    typename std::enable_if<std::is_same<Cur, Last>::value>::type *dummy = 0) {
  if (reset)
    d.template updateHelper<Type1>(sys);
}
} // !namespace detail

/**
 * Keeps a list of edited particles and calls reactions on them when updated
 */
template <class SysProps, class Reaction, class... Types> struct DirtyList {
  static constexpr int dims = SysProps::dims;
  using Coords = typename SysProps::Coords;
  using Pid = typename SysProps::Pid;
  using Species = typename SysProps::Species;
  using Tuple = std::tuple<Types...>;
  using Type1 = std::tuple_element_t<0, Tuple>;
  using TypeN = std::tuple_element_t<sizeof...(Types)-1, Tuple>;

  DirtyList(Reaction r, bool notAllTypes = false) : _reaction{r}, _notAll{notAllTypes} {}

  template <class T> static constexpr auto indexOfType() {
    return meta::typeIndex<T, Tuple>();
  }

  template <class T> constexpr auto &dirtMap() {
    assert (_notAll || (indexOfType<T>() < _dirty.size()));
    return std::get<indexOfType<T>()>(_dirty);
  }

  template <class T> auto dirtCount(const Species &sp) {
    size_t i = meta::typeIndex<T, Tuple>();
    if (i < _dirty.size()) {
      auto &dirt = _dirty[i]; // dirtMap<T>();
      auto it = dirt.find(sp);
      if (it != std::end(dirt))
        return it->second.size();
    }
    return std::size_t(0);
  }

  template <class Sys, class T> void added(const Sys& sys, const T &p, const Pid &pid) {
    _addDirty<indexOfType<T>()>(pid);
  }

  template <class Sys, class T>
  void editedCoords(const Sys& sys, const T &p, const Pid &pid, const Coords &pPrior) {
    _addDirty<indexOfType<T>()>(pid);
  }

  template <class Sys, class T>
  void editedSize(const Sys &, T &, const Pid &pid, const int &) {
    _addDirty<indexOfType<T>()>(pid);
  }

  template <class Sys, class T>
  void edited(const Sys& sys, const T &p, const Pid &pid, const T &pPrior) {
    _addDirty<indexOfType<T>()>(pid);
  }

  template <class Sys, class T>
  void edited(const Sys& sys, const T &p, const Pid &pid) {
    _addDirty<indexOfType<T>()>(pid);
  }

  template <class Sys, class T> void removed(const Sys &, const T &p, const Pid &pid) {
    _removeDirty<indexOfType<T>()>(pid);
  }

  template <class Sys, class T> void erased(const Sys &, const T &p, const Pid &pid) {
    _erasedDirty<indexOfType<T>()>(pid);
  }

  void blackList(Species sp) {
    _blackList.insert(sp);
  }

  template <class Sys> void update(Sys &sys, double t, double dt) {
    sys.commit();
    updateHelper<Type1>(sys);
  }

  template <class T, class Sys> auto updateHelper(Sys &sys) {
    auto &dirt = dirtMap<T>();
    auto noReaction = true; // assumption to begin with
    for (auto species = std::begin(dirt);
         species != std::end(dirt);) { // for all species
      auto temp = std::move(
          species->second); // dirty may be manipulated by reaction itself
      species->second.clear();
      for (auto pit = std::begin(temp); pit != std::end(temp);
           ++pit) { // for all SystemTuple
        auto pid = std::make_pair(species->first, *pit);
        if (_reaction(sys, pid, sys.template particle<T>(pid))) {
          // at this point species second might have something
          species->second = fillDirty(sys.template toDel<T>(species->first),
                                       species->second, ++pit, std::end(temp));
          sys.commit();
          //_checkSanity<Sys, T>(sys);
          noReaction = false;
          break;
        }
      }
      if (noReaction) {
        species = dirt.erase(species);
      } else {
        if (species->second.empty()) {
          dirt.erase(species); // TODO: check if needed
        }
        species = std::begin(dirt);
        break;
      }
    }
    detail::updateNext<Type1, T, TypeN, Tuple>(*this, sys, !noReaction);
  }

  // add only if not already exist
  template <class It>
  auto fillDirty(const std::set<size_t> &toDel, std::set<size_t> &v,
                  const It &first, const It &last) {
    std::set<size_t> res;
    v.insert(first, last);
    std::set_difference(begin(v), end(v), begin(toDel), end(toDel),
                        std::inserter(res, end(res)));
    return res;
  }



protected:
  template <class Sys, class T> void _checkSanity(Sys& sys) {
    auto& del = sys.template toDel<T>();
    for (auto& it : del) assert(it.second.size() == 0 && "del not clear yet");
    auto& dirt = dirtMap<T>();
    for (auto& sp : dirt) {
      auto sz = sys.template particles<T>(sp.first).size();
      for (auto& index : sp.second) {
        assert(index < sz && "index not less than sz");
      }
    }
  }
  template <size_t i> void _addDirty(const Pid &pid) {
    assert (_notAll || (i < _dirty.size()));
    if (i < sizeof...(Types)) {
      if (_blackList.find(pid.first) != std::end(_blackList)) return;
      auto &it = _dirty[i][pid.first]; // std::get<i>(_dirty)[pid.first];
      it.insert(pid.second);
    }
  }

  // called when pid is deleted and added to toDel
  // does not erase the species as it may invalidate the map iterator
  // while update()
  template <size_t i> void _removeDirty(const Pid &pid) {
    assert (_notAll || (i < _dirty.size()));
    if (i < sizeof...(Types)) {
      auto &dirt = _dirty[i]; // std::get<i>(_dirty);
      auto it = dirt.find(pid.first);
      if (it != std::end(dirt)) {
        auto jt = it->second.find(pid.second); // std::find(it->second.begin(),
                                               // it->second.end(), pid.second);
        if (jt != std::end(it->second)) {
          // std::swap(*jt, it->second[it->second.size() - 1]);
          // it->second.pop_back();
          it->second.erase(jt);
        }
      }
    }
  }

  template <size_t i> void _erasedDirty(const Pid &pid) {
    assert (_notAll || (i < _dirty.size()));
    if (i < sizeof...(Types)) {
      auto &dirt = _dirty[i]; // std::get<i>(_dirty);
      auto kt = dirt.find(pid.first);
      if (kt != std::end(dirt)) {
        std::set<size_t> temp;
        for (auto mt = kt->second.begin(); mt != kt->second.end();) {
          if (*mt > pid.second) {  // TODO: optimize since sorted
            temp.insert(*mt - 1);
            mt = kt->second.erase(mt);
          } else {
            ++mt;
          }
        }
        kt->second.insert(begin(temp), end(temp));
      }
    }
  }

  using DirtyT = std::array<std::unordered_map<Species, std::set<size_t>>,
                            sizeof...(Types)>;
  DirtyT _dirty; // // vector because modified in commit
  std::unordered_set<Species> _blackList;
  // using Reaction = std::function<bool(Sys &, const Pid&)>;
  // using Reactions = std::array<std::unordered_map<Species,
  // std::vector<Reaction>>, sizeof...(Types)>;
  Reaction _reaction;
  bool _notAll;
};
} // namespace !gkmc

#endif // !GKMC_DirtyList_H

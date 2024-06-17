// SystemTuple.hpp

#ifndef GKMC_SystemTuple_H
#define GKMC_SystemTuple_H

#include <algorithm>
#include <array>
#include <assert.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include <helper/meta.hpp>

namespace gkmc {

template <class SystemProps, class Hook, class... Types> class SystemTuple {
public:
  static constexpr int dims = SystemProps::dims;
  static_assert(dims > 0, "Negative value for dimension is not valid");
  using Props = SystemProps;
  using Coords = typename SystemProps::Coords;
  using value_type = typename SystemProps::real_type;
  using Pid = typename SystemProps::Pid;
  using Species = typename SystemProps::Species;
  using TypesTuple = std::tuple<Types...>;
  // using Type1 = std::tuple_element<0, Tuple>;
  template <class T>
  using particleM = std::unordered_map<Species, std::vector<T>>;
  using SystemTupleT = std::tuple<particleM<Types>...>;
  SystemTupleT _SystemTuple;
  SystemProps _sysProps;
  Hook _hook;
  size_t _pkidCount = 0;
  using DelT = std::array<std::unordered_map<Species, std::set<size_t>>,
                          sizeof...(Types)>;
  DelT _toDel;
  using Names =
      std::array<std::unordered_map<Species, std::string>, sizeof...(Types)>;
  Names _names;

public:
  template <class T>
  class iterator {
  public:
    using self_type = iterator;
    using ParticleType = T;
    using reference =  ParticleType&;
    using pointer =  ParticleType*;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = int;
    using it_type = typename std::vector<T>::iterator;
    iterator(it_type ptr) : ptr_(ptr) { }
    self_type operator++() { self_type i = *this; ptr_++; return i; }
    self_type operator++(int junk) { ptr_++; return *this; }
    reference operator*() { return *(ptr_); }
    it_type &operator->() { return ptr_; }
    bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
    bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
  private:
    it_type ptr_;
  };

  template <class T>
  auto begin(const Species &sp) {
    return iterator<T>(std::begin(particles<T>(sp)));
  }

  template <class T>
  auto end(const Species &sp) {
    return iterator<T>(std::end(particles<T>(sp)));
  }

  SystemTuple(SystemProps sysProps, Hook hook)
      : _sysProps{sysProps}, _hook{std::forward<Hook>(hook)} {}
  auto &props() { return _sysProps; };

  template <class T> const std::set<size_t> &toDel(const Species &sp) {
    assert(isSpecies<T>(sp));
    auto &del = std::get<indexOfType<T>()>(_toDel);
    return del[sp];
  }

  template <class T> const auto &toDel() {
    auto &del = std::get<indexOfType<T>()>(_toDel);
    return del;
  }

  template <size_t index> const auto &namesMap() const {
    return std::get<index>(_names);
  }

  template <class T> void speciesName(const Species &sp, const std::string &nm) {
    assert(isSpecies<T>(sp));
    auto &namesMap = std::get<indexOfType<T>()>(_names);
    namesMap[sp] = nm;
  }

  template <class T> const std::string& speciesName(const Species &sp) {
    auto &namesMap = std::get<indexOfType<T>()>(_names);
    return namesMap[sp];
  }

  template <class T> static constexpr size_t indexOfType() {
    return meta::typeIndex<particleM<T>, SystemTupleT>();
  }

  template <class T> constexpr auto mapOfType() -> particleM<T> & {
    return std::get<indexOfType<T>()>(_SystemTuple);
  }

  auto cookPid(const Species &sp, const size_t &i) const {
    return std::make_pair(sp, i);
  }
  template <class T> auto isSpecies(const Species &sp) {
    auto it = particles<T>().find(sp);
    return !(it == std::end(particles<T>()));
  }

  template <class T> const auto &getSpecies(const Pid &pid) {
    return pid.first;
  }

  template <class T> auto &particles(const Species &sp) {
    assert(isSpecies<T>(sp) && "No species for the particle.");
    return mapOfType<T>()[sp];
  }

  template <class T> const auto &particles() { return mapOfType<T>(); }

  template <class T> void addSpecies(const Species &sp) { mapOfType<T>()[sp]; }
  template <class T> void addSpecies(const Species &sp, const std::string &name) { 
    mapOfType<T>()[sp]; 
    auto &namesMap = std::get<indexOfType<T>()>(_names);
    namesMap[sp] = name;
  }
  template <class T> auto add(const Species &sp, const T &particle) {
    auto &pv = particles<T>(sp);
    pv.push_back(particle);
    pv[pv.size() - 1].id(_pkidCount++);
    auto pid = Pid{sp, pv.size() - 1};
    _hook.added(*this, pv[pv.size() - 1], pid);
    return pid;
  }

  /*!
  * Edit particle
  */
  template <class T> void edited(const Pid &pid) {
    auto &p = particle<T>(pid);
    _hook.edited(*this, p, pid);
  }

  /*!
  * Edit particle
  */
  template <class T> void edited(const Pid &pid, const T &priorParticle) {
    auto &p = particle<T>(pid);
    _hook.edited(*this, p, pid, priorParticle);
  }

  /*!
  * Change coordinates of a particle
  */
  template <class T> void changeCoords(const Pid &pid, const Coords &c) {
    auto &p = particle<T>(pid);
    auto prior = p.coords();
    p.coords(c);
    _hook.editedCoords(*this, p, pid, prior);
  }

  /*!
  * Change coordinates of a particle
  */
  template <class T> void changeSize(const Pid &pid, const int &c) {
    auto &p = particle<T>(pid);
    auto prior = p.size();
    p.size(c);
    _hook.editedSize(*this, p, pid, prior);
  }

  template <class T> auto &particle(const Pid &pid) {
    auto &v = particles<T>(pid.first);
    assert(pid.second < v.size() && "Wrong pid for particle.");
    /*
    if (pid.second < v.size()) {
      if (v.empty()) v.push_back(T{});
      return v[0];
    }
    */
    return v[pid.second];
  }
  /*!
  * Change coordinates of a particle and add wrap coordinate count
  * for any of the directions.
  */
  template <class T>
  void changeCoordsAdd(const Pid &pid, const value_type &c, const int &add,
                       const size_t &index) {
    auto &p = particle<T>(pid);
    auto prior = p.coords();
    p.coords(c, index);
    p.addCount(add, index);
    _hook.editedCoords(*this, p, pid, prior);
  }

  /*!
  *
  */
  template <class T> auto remove(const Pid &pid) {
    auto &del = std::get<indexOfType<T>()>(_toDel);
    del[pid.first].insert(pid.second);
    auto &p = particle<T>(pid);
    _hook.removed(*this, p, pid);
  }

  void commit() { (bool[]){commitHelper<Types>()...}; }

  /*!
   * Deletes all the SystemTuple in toDel.
   * Also, changes dirty Pids accordingly.
   * TODO: check if remove idiom can be used
   */
  template <class T> bool commitHelper() {
    auto &del = std::get<indexOfType<T>()>(_toDel);
    auto &pm = mapOfType<T>();
    for (auto &it : del) {
      for (auto jt = std::rbegin(it.second); jt != std::rend(it.second); ++jt) {
        _hook.erased(*this, pm[it.first][*jt], Pid{it.first, *jt});
        if (*jt != pm[it.first].size() - 1) {
          pm[it.first][*jt] = pm[it.first][pm[it.first].size() - 1];
        }
        pm[it.first].pop_back();
      }
    }
    del.clear();
    return true;
  }

  /*!
   * commits and executes the reactions recursively until dirty SystemTuple exist.
   */
  void update(double t, double dt) { _hook.update(*this, t, dt); }
};

template <size_t index, class P, class Nmap>
bool printSystemTupleHelper(P &pMap, const Nmap& names) {
  for (auto &it : pMap) {
    if (it.second.empty()) continue;
    auto nm = names.find(it.first);
    std::cout << "Type: " << index << "\tSpecies: ";
    if (nm != std::end(names)) std::cout << nm->second;
    std::cout <<" -- \n";
    for (auto &jt : it.second) {
      printParticle(jt);
      std::cout << '\n';
    }
    std::cout << '\n';
  }
  return false;
}

template <size_t i, class Sys> void printSystemTupleExpander(const Sys &) {}
template <size_t i, class Sys, class Type, class... Types>
void printSystemTupleExpander(const Sys &ps, Type &t, Types &... ts) {
  printSystemTupleHelper<i>(t, ps.template namesMap<i>());
  printSystemTupleExpander<i + 1>(ps, ts...);
}

template <class Prop, class Hook, class... Types>
auto printSystemTuple(SystemTuple<Prop, Hook, Types...> &ps) {
  printSystemTupleExpander<0>(ps, ps.template particles<Types>()...);
}

} // namespace !gkmc

#endif // !GKMC_SystemTuple_HremoveParticle

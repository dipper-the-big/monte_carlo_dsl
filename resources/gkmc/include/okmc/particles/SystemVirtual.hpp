// SystemVirtual.hpp

#ifndef GKMC_SystemVirtual_H
#define GKMC_SystemVirtual_H

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
#include <okmc/reactions/ReactionCollection.hpp>
#include <okmc/particles/SimpleParticle.hpp>

namespace gkmc {

template <class SystemProps> class SystemVirtual {
public:
  static constexpr int dims = SystemProps::dims;
  static_assert(dims > 0, "Negative value for dimension is not valid");
  using value_type = typename Sysprops::real_type;
  using Props = SystemProps;
  using Coords = typename SystemProps::Coords;
  using Pid = typename SystemProps::Pid;
  using Species = typename SystemProps::Species;
  using Type = SimpleParticleVirtual<dims, 1>;
  using TypesTuple = std::tuple<Type>;
  using ParticlePtr = std::unique_ptr<Type>;
  template <class T>
  using particleM = std::unordered_map<Species, std::vector<T>>;
  using SystemTupleT = particleM<ParticlePtr>;
  using Hook = ReactionCollection<SystemVirtual<Props>, Type>&;
  SystemTupleT _SystemTuple;
  SystemProps _sysProps;
  Hook _hook;
  size_t _pkidCount = 0;
  using DelT = std::unordered_map<Species, std::set<size_t>>;
  DelT _toDel;
  using Names = std::unordered_map<Species, std::string>;
  Names _names;

public:

  class iterator {
  public:
    using self_type = iterator;
    using reference = Type&;
    using pointer = Type*;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = int;
    using it_type = typename std::vector<ParticlePtr>::iterator;
    iterator(it_type ptr) : ptr_(ptr) { }
    self_type operator++() { self_type i = *this; ptr_++; return i; }
    self_type operator++(int junk) { ptr_++; return *this; }
    reference operator*() { return *(ptr_); }
    auto& operator->() { return *(ptr_); }
    // pointer operator->() { return ptr_->get(); }
    bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
    bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
  private:
    it_type ptr_;
  };

  SystemVirtual(SystemProps sysProps, Hook hook)
      : _sysProps{sysProps}, _hook{hook} {}
  auto &props() { return _sysProps; };

  template <class T> const std::set<size_t> &toDel(const Species &sp) {
    assert(isSpecies<T>(sp));
    return _toDel[sp];
  }

  template <class T> const auto &toDel() {
    return _toDel;
  }

  template <size_t index> const auto &namesMap() const {
    return _names;
  }

  template <class T> void name(const Species &sp, const std::string &nm) {
    assert(isSpecies<T>(sp));
    _names[sp] = nm;
  }

  template <class T> static constexpr size_t indexOfType() {
    return 0;
  }

  template <class T> constexpr auto& mapOfType() /*-> particleM<T> & */{
    return _SystemTuple;
  }

  template <class T>
  auto begin(const Species &sp) {
    return iterator(std::begin(SystemTuple<T>(sp)));
  }

  template <class T>
  auto end(const Species &sp) {
    return iterator(std::end(SystemTuple<T>(sp)));
  }

  auto cookPid(const Species &sp, const size_t &i) const {
    return std::make_pair(sp, i);
  }
  template <class T> auto isSpecies(const Species &sp) {
    auto it = SystemTuple<T>().find(sp);
    return !(it == std::end(SystemTuple<T>()));
  }

  template <class T> const auto &getSpecies(const Pid &pid) {
    return pid.first;
  }

  template <class T> auto &SystemTuple(const Species &sp) {
    assert(isSpecies<T>(sp) && "No species for the particle.");
    return mapOfType<T>()[sp];
  }

  template <class T> const auto &SystemTuple() { return mapOfType<T>(); }

  template <class T> void addSpecies(const Species &sp) { mapOfType<T>()[sp]; }
  template <class T> void addSpecies(const Species &sp, const std::string &name) { 
    //SystemTuple<T>()[sp]; 
    mapOfType<T>()[sp];
    _names[sp] = name;
  }
  template <class T> auto add(const Species &sp, const T &particle) {
    auto &pv = SystemTuple<T>(sp);
    Type* x = new T(particle);
    pv.emplace_back(std::unique_ptr<Type>{x});
    pv[pv.size() - 1]->id(_pkidCount++);
    auto pid = Pid{sp, pv.size() - 1};
    _hook.added(*this, *(pv[pv.size() - 1]), pid);
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
  * Change size of a particle
  */
  template <class T> void changeSize(const Pid &pid, const int &c) {
    auto &p = particle<T>(pid);
    auto prior = p.size();
    p.size(c);
    _hook.editedSize(*this, p, pid, prior);
  }

  template <class T> auto &particle(const Pid &pid) {
    auto &v = SystemTuple<T>(pid.first);
    assert(pid.second < v.size() && "Wrong pid for particle.");
    /*
    if (pid.second < v.size()) {
      if (v.empty()) v.push_back(T{});
      return v[0];
    }
    */
    return *(v[pid.second]);
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
    auto &del = _toDel;
    del[pid.first].insert(pid.second);
    auto &p = particle<T>(pid);
    _hook.removed(*this, p, pid);
  }

  void commit() { (bool[]){commitHelper<Type>()}; }

  /*!
   * Deletes all the SystemTuple in toDel.
   * Also, changes dirty Pids accordingly.
   * TODO: check if remove idiom can be used
   */
  template <class T> bool commitHelper() {
    auto &del = _toDel;
    auto &pm = mapOfType<T>();
    for (auto &it : del) {
      for (auto jt = std::rbegin(it.second); jt != std::rend(it.second); ++jt) {
        _hook.erased(*this, *(pm[it.first][*jt]), Pid{it.first, *jt});
        if (*jt != pm[it.first].size() - 1) {
          pm[it.first][*jt] = std::move(pm[it.first][pm[it.first].size() - 1]);
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

} // namespace !gkmc

#endif // !GKMC_SystemTuple_HremoveParticle
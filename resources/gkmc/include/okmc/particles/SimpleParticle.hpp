// SimpleParticle.hpp

#ifndef GKMC_SimpleParticle_H
#define GKMC_SimpleParticle_H

#include <array>

namespace gkmc {

template <class C, int sz = 1> class SimpleParticle {
public:
  using Coords = C;
  using value_type = typename Coords::value_type;
  SimpleParticle(const Coords &c) : _coords{c} {}
  SimpleParticle() = default;
  const auto &id() const { return _id; }
  const auto &coords() const { return _coords; }
  void id(const size_t &id) { _id = id; }
  void coords(const Coords &c) { _coords = c; }
  void coords(value_type val, size_t i) { _coords[i] = val; }
  bool operator==(const SimpleParticle &p) const {
    return _coords == p._coords && _id == p._id;
  }
  static constexpr int size() { return sz; }
  void size(int val) {
    assert (val == sz && "setting size different from preset size in simple particle");
  }
private:
  size_t _id = 0;
  Coords _coords;
};

template <class Base> class PbcParticle : public Base {
public:
  using Coords = typename Base::Coords;
  static constexpr auto dims = std::tuple_size<Coords>::value; //Coords::size();
  PbcParticle(const Coords &c, int sz) : Base{c}, _count{} { Base::size(sz); }
  PbcParticle(const Coords &c) : Base{c}, _count{} {}
  void addCount(int add, size_t index) { _count[index] += (size_t)add; }
  const auto &count() const { return _count; }

private:
  std::array<int, dims> _count;
};

template <class Base> class SizeParticle : public Base {
public:
  using Coords = typename Base::Coords;
  SizeParticle(const Coords &c, int size = 1)
      : Base{c}, _size{size} {}
  void size(int val) { _size = val; }
  int size() const { return _size; }

private:
  int _size;
};

template <class C, int sz = 1> class SimpleParticleVirtual {
public:
  using Coords = C;
  using value_type = typename Coords::value_type;
  SimpleParticleVirtual(const Coords &c) : _coords{c} {}
  SimpleParticleVirtual() = default;
  virtual ~SimpleParticleVirtual() = default;
  const auto &id() const { return _id; }
  const auto &coords() const { return _coords; }
  void id(const size_t &id) { _id = id; }
  void coords(const Coords &c) { _coords = c; }
  void coords(value_type val, size_t i) { _coords[i] = val; }
  virtual bool operator==(const SimpleParticleVirtual &p) const {
    return _coords == p._coords && _id == p._id;
  }
  virtual int size() const { return sz; }
  virtual void addCount(int add, size_t index) {}
  virtual void size(int val) {
    assert (val == sz && "setting size different from preset size in simple particle");
  }
private:
  size_t _id = 0;
  Coords _coords;
};

template <class Coords, int sz>
void printParticle(const SimpleParticle<Coords, sz>& particle) {
  std::cout<<particle.id()<<'\t'<<sz<<'\t';
  for (auto& it : particle.coords()) {
    std::cout<<it<<'\t'; 
  }
}

template <class Base>
void printParticle(const PbcParticle<Base>& particle) {
  const Base& b = particle;
  printParticle(b);
  for (auto& it : particle.count()) {
    std::cout<<it<<'\t'; 
  }
}

template <class Base>
void printParticle(const SizeParticle<Base>& particle) {
  const Base& b = particle;
  printParticle(b);
  std::cout<<particle.size()<<'\t'; 
}

} // !namespace gkmc

#endif // !GKMC_SimpleParticle_H

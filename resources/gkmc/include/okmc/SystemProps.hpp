// SystemProps.hpp

#ifndef GKMC_SystemProps_H
#define GKMC_SystemProps_H

#include <array>
#include <map>

namespace gkmc {

template<int D, class Sp, class realT = double>
class SystemProps {
public:
  static constexpr int dims = D;
  using real_type = realT;
  using Coords = std::array<real_type, dims>;
  using Species = Sp;
  using Pid = std::pair<Species, std::size_t>;

  void box(int unitCells, real_type latticeConstant) {
    _latticeConstant = latticeConstant;
    for (auto &it : _unitCells)
      it = unitCells;
    for (auto &it : _boxDimensions)
      it = latticeConstant * unitCells;
  }
  void box(Coords unitCells, real_type latticeConstant) {
    _latticeConstant = latticeConstant;
    _unitCells = unitCells;
    for (auto i = 0; i < dims; ++i) {
      _boxDimensions[i] = latticeConstant * unitCells[i];
    }
  }
  void boxDimensions(const Coords& boxDimensions) { _boxDimensions = boxDimensions; }
  void temperature(real_type val) { _temperature = val; }
  void latticeConstant(real_type lc) { _latticeConstant = lc; }

  // getters
  const auto &unitCells() const { return _unitCells; }
  const auto &latticeConstant() const { return _latticeConstant; }
  const auto &boxDimensions() const { return _boxDimensions; }
  const auto &temperature() const { return _temperature; }

private:
  Coords _boxDimensions;
  Coords _unitCells;
  real_type _latticeConstant;
  real_type _temperature;
};
} // namespace !gkmc

#endif // !GKMC_SystemProps_H

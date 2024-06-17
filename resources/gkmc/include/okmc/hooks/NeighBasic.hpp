// NeighBasic.hpp

#ifndef GKMC_NEIGHBASIC_H
#define GKMC_NEIGHBASIC_H

#include <vector>

#include <helper/dist.hpp>
#include <okmc/hooks/Dummy.hpp>

namespace gkmc {

/**
 * Brute force neighbour search.
 * */
class NeighBasic : public Dummy {
private:
  bool _isMirror;

public:

  NeighBasic(const bool &isMirror) : _isMirror{isMirror} {};
// TODO: call the other overload from inside the loop
  template <class P1, class P2, class Sys, class CapRadiusFn>
  auto neighbour(Sys &SystemTuple, const typename Sys::Pid &pid,
                 CapRadiusFn &&capRadius) const {
    // const auto &all = SystemTuple.template particles<P2>();
    auto &all = SystemTuple.template particles<P2>();
    const auto &p = SystemTuple.template particle<P1>(pid);
    auto pred = (_isMirror)
                    ? AreClose<typename Sys::Coords>{SystemTuple.props().boxDimensions()}
                    : AreClose<typename Sys::Coords>{};
    //for (const auto &lot : all) {
    for (auto &lot : all) {
      size_t i = 0;
      auto& del = SystemTuple.template toDel<P2>(lot.first);
      auto it = SystemTuple.template begin<P2>(lot.first);
      auto last = SystemTuple.template end<P2>(lot.first);
      for (; it != last; ++it) {
        if (it->id() == p.id() || del.find(i) != std::end(del)) {
          ++i;
          continue;
        }
        auto tol = capRadius(SystemTuple, p, *it, pid, SystemTuple.cookPid(lot.first, i));
        if (pred(p.coords(), it->coords(), tol)) {
          return make_pair(true, SystemTuple.cookPid(lot.first, i));
        }
        ++i;
      }
    }
    return make_pair(false, typename Sys::Pid{});
  }

  template <class P1, class P2, class Sys, class CapRadiusFn>
  auto neighbour(Sys &SystemTuple, const typename Sys::Pid &pid,
                 CapRadiusFn &&capRadius,
                 const typename Sys::Species &sp) {
    auto pred = (_isMirror)
                    ? AreClose<typename Sys::Coords>{SystemTuple.props().boxDimensions()}
                    : AreClose<typename Sys::Coords>{};
    const auto &p = SystemTuple.template particle<P1>(pid);
    // const auto &lot = SystemTuple.template particles<P2>(sp);
    auto& del = SystemTuple.template toDel<P2>(sp);
    size_t i = 0;
    auto it = SystemTuple.template begin<P2>(sp);
    auto last = SystemTuple.template end<P2>(sp);
    for (; it != last; ++it) {
      if (it->id() == p.id() || del.find(i) != std::end(del)) {
        ++i;
        continue;
      }
      auto tol = capRadius(SystemTuple, p, *it, pid, SystemTuple.cookPid(sp, i));
      if (pred(p.coords(), it->coords(), tol)) {
        return make_pair(true, SystemTuple.cookPid(sp, i));
      }
      ++i;
    }
    return make_pair(false, typename Sys::Pid{});
  }
};

} // namespace !gkmc

#endif // !GKMC_NeighBasic_H

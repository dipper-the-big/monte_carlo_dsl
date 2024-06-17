// ReactionInactiveCounter.hpp

#ifndef GKMC_ReactionInactiveCounter_H
#define GKMC_ReactionInactiveCounter_H

namespace gkmc {

/**
 * A wrapper that counts inactive streak of a reaction. Used in the old Kai's notes
 * problem in the predicate condition for bkl.
 * */
template <class Reaction>
class ReactionInactiveCounter {
public:
  ReactionInactiveCounter(Reaction r, size_t& counter)
      : _reaction{r}, _curInactiveSteps{counter} {}
  template <class SystemTuple, class Particle>
  auto operator()(SystemTuple& ps, const typename SystemTuple::Pid &pid, Particle& p) {
    bool res = _reaction(ps, pid, p);
    if (res) _curInactiveSteps = 0;
    return res;
  }
private:
  Reaction _reaction;
  size_t& _curInactiveSteps;
};
}
#endif
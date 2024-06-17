#include <okmc/hooks/Dummy.hpp>
#include <okmc/hooks/Both.hpp>

#ifndef GKMC_ReactionCollection_H
#define GKMC_ReactionCollection_H

namespace gkmc {

/**
 * Reaction collection that can be used for runtime polymorphism in reactions.
 * */
template <template <class> class System, class Props, class Particle, class Other = Dummy>
class ReactionCollection : public Both<DirtyList<Props, int, Particle>, Other> {
public:
  using Sys = System<ReactionCollection&>;
  //using Props = Sys::Props;
  using Dirt = DirtyList <Props, int, Particle>;
  using Base = Both<Dirt, Other>;

  using Pid = typename Sys::Pid;
  using Species = typename Sys::Species;
  using Reaction = std::function<bool(Sys &, const Pid&, Particle &)>;

  ReactionCollection(Other other = Other{}) : Base{Dirt{0}, other} {}//_dirty{0} { }


  void addSpecies(Species sp) {
    _reactions[sp];
  }

  void addReaction(Reaction reaction, Species sp) {
    _reactions[sp].push_back(reaction);
  }

  void addReaction(Reaction reaction) {
    for (auto& it : _reactions) {
      it.second.push_back(reaction);
    }
  }
  void update(Sys &sys, double t, double dt) {
    sys.commit();
    auto &dirt = Base::first().template dirtMap<Particle>();
    auto noReaction = true; // assumption to begin with
    for (auto species = std::begin(dirt); species != std::end(dirt);) { // for all species
      auto temp = std::move(species->second); // dirty may be manipulated by reaction itself
      species->second.clear();
      for (auto pit = std::begin(temp); pit != std::end(temp); ++pit) { // for all SystemTuple
        auto pid = std::make_pair(species->first, *pit);
        if (_reactions.find(species->first) == std::end(_reactions)) continue;
        for (auto& rit : _reactions[species->first]) {
          if (rit(sys, pid, sys.template particle<Particle>(pid))) {
            // at this point species second might have something
            species->second = Base::first().fillDirty(sys.template toDel<Particle>(species->first),
                                       species->second, ++pit, std::end(temp));
            sys.commit();
            //_checkSanity<Sys, Particle>(sys);
            noReaction = false;
            break;
          }
        }
        if (noReaction == false) break;
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
    Base::second().update(sys, t, dt);
  }

private: 
  std::unordered_map<Species, std::vector<Reaction>> _reactions;
};

}

#endif
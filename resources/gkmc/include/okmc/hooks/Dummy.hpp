// Dummy.hpp

#ifndef GKMC_Dummy_H
#define GKMC_Dummy_H

namespace gkmc {

/**
 *  A no operation hook that can also be used for inheriting from
*/
class Dummy {
public:
  Dummy() = default;
  template <class Sys, class Particle>
  void added(const Sys &, Particle &, const typename Sys::Pid &) {}
  template <class Sys, class Particle>
  void edited(const Sys &, Particle &, const typename Sys::Pid &) {}
  template <class Sys, class Particle>
  void editedCoords(const Sys &, Particle &, const typename Sys::Pid &,
                    const typename Sys::Coords &) {}
  template <class Sys, class Particle>
  void removed(const Sys &, Particle &, const typename Sys::Pid &) {}
  template <class Sys, class Particle>
  void erased(const Sys &, Particle &, const typename Sys::Pid &) {}
  template <class Sys> void update(Sys &, double, double) {}
  template <class Sys, class Particle>
  void editedSize(const Sys &, Particle &, const typename Sys::Pid &,
                    const int &) {}
};

} // namespace !gkmc

#endif // !GKMC_Dummy_H

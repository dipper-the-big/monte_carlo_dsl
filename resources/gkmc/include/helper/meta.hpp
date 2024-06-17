#ifndef GKMC_Meta_H
#define GKMC_Meta_H

namespace gkmc {
namespace meta {

// gettting index of the type from a tuple.
template <size_t i, class T, class Tuple> struct Index;

template <size_t i, class T, class... Types>
struct Index<i, T, std::tuple<T, Types...>> {
  static constexpr std::size_t value = i;
};

template <size_t i, class T>
struct Index<i, T, std::tuple<>> {
  static constexpr std::size_t value = i + 1;  // Err. condition
};

template <size_t i, class T, class U, class... Types>
struct Index<i, T, std::tuple<U, Types...>>
    : Index<i + 1, T, std::tuple<Types...>> {};

template <class T, class Tuple> constexpr auto typeIndex() {
  return Index<0, T, Tuple>::value;
}
}
}

#endif //!GKMC_Meta_H
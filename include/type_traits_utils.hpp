#ifndef ZONOTOPE_TYPE_TRAITS_UTILS_HPP_
#define ZONOTOPE_TYPE_TRAITS_UTILS_HPP_

#include <limits>
#include <type_traits>

namespace zonotope {
  
template<typename A, typename B = A>
using enable_if_integer_t = typename std::enable_if<std::numeric_limits<A>::is_integer, B>::type;

template<typename A, typename B = A>
using enable_if_not_integer_t = typename std::enable_if<!std::numeric_limits<A>::is_integer, B>::type;

} // namespace zonotope

#endif // ZONOTOPE_TYPE_TRAITS_UTILS_HPP_
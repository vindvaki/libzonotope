#ifndef ZONOTOPE_HALFSPACES_HPP_
#define ZONOTOPE_HALFSPACES_HPP_

#include "zonotope.hpp"

#include "combination_traversal.hpp"
#include "combination_kernel_container.hpp"
#include "event_point_2.hpp"

namespace zonotope {

template<typename Zonotope_data_t,
         typename Output_functor_t,
         typename Kernel_number_t>
struct Halfspaces_traversal_output_functor {
  
  const Zonotope_data_t& z;
  Output_functor_t& f;
  
  Halfspaces_traversal_output_functor(const Zonotope_data_t& z, Output_functor_t& f)
    : z(z), f(f)
    { }

  bool operator() (const Combination_kernel_container<Zonotope_data_t, Kernel_number_t>& c) {
    if ( c.size() == ( z.dimension - 2 ) ) {
      handle_event_points(z, c, f);
      return true;
    }
    return false;
  }
};

template <typename Zonotope_data_t,
          typename Output_functor_t,
          typename Kernel_number_t>
void zonotope_halfspaces(const Zonotope_data_t& z, Output_functor_t& f) {
  Halfspaces_traversal_output_functor<Zonotope_data_t, Output_functor_t, Kernel_number_t> g(z, f);
  Combination_kernel_container<Zonotope_data_t, Kernel_number_t> empty_combination(z, z.dimension - 2);
  traverse_combinations(empty_combination, g);
}

} // namespace zonotope

#endif // ZONOTOPE_HALFSPACES_HPP_

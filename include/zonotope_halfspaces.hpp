#ifndef ZONOTOPE_HALFSPACES_HPP_
#define ZONOTOPE_HALFSPACES_HPP_

#include "combination_traversal.hpp"
#include "hyperplane.hpp"
#include "zonotope_halfspaces_output_functor.hpp"
#include "combination_kernel_container.hpp"
#include "preprocess_generators.hpp"

// External dependencies
#include <vector>
#include <set>
#include <gmpxx.h>

namespace zonotope {

/**
 * @brief Generic construction of the set of halfspaces
 *
 */
template <typename User_number_t,
          typename Internal_number_t = mpz_class,
          typename Halfspaces_container_t = std::set<Hyperplane<User_number_t> > >
void zonotope_halfspaces (
  const std::vector<std::vector<User_number_t> >& generators_in,
  Halfspaces_container_t& halfspaces )
{
  typedef Combination_kernel_container<Internal_number_t> Combination_container_t;

  typedef Container_output_functor<Halfspaces_container_t,
                                   Hyperplane<Internal_number_t>,
                                   Hyperplane<User_number_t> >
      Container_output_functor_t;

  Container_output_functor_t Halfspaces_container_output_fn ( halfspaces );

  typedef Zonotope_halfspaces_output_functor<Internal_number_t,
                                             Combination_container_t,
                                             Container_output_functor_t>
      Traversal_output_functor_t;

  const int d = generators_in[0].size();

  std::vector<std::vector<Internal_number_t> > internal_generators;

  preprocess_generators<User_number_t, Internal_number_t> (generators_in,
                                                           internal_generators);

  Combination_container_t empty_combination (internal_generators, d-1);

  Traversal_output_functor_t Traversal_output_fn (internal_generators,
                                                  Halfspaces_container_output_fn);

  traverse_combinations<Combination_container_t, Traversal_output_functor_t>
    (empty_combination, Traversal_output_fn);
  // we only traverse up to (d-2)-combinations because after that,
  // Traversal_output_fn takes over and traverses the (d-1)-child-combinations
  // in a manner specific to the halfspace traversal.
}

} // namespace zonotope

#endif // ZONOTOPE_HALFSPACES_HPP_

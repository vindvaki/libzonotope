#ifndef ZONOTOPE_HALFSPACES_HPP_
#define ZONOTOPE_HALFSPACES_HPP_

#include "hyperplane.hpp"
#include "combination_traversal.hpp"
#include "zonotope_halfspaces_output_functor.hpp"
#include "combination_kernel_container.hpp"

#include <vector>
#include <set>

/**
 * @brief Generic construction of the set of halfspaces
 *
 * @tparam Combination_container A combination container that also
 *                               maintains the kernel of the
 *                               combination
 */
template <typename NT,
          typename Combination_container = Combination_kernel_container<NT> >
std::set<Hyperplane<NT> >
zonotope_halfspaces (const std::vector<std::vector<NT> >& generators) {

  typedef Zonotope_halfspaces_output_functor<NT, Combination_container> Output_functor;

  const int d = generators[0].size();

  Combination_container empty_combination (generators, d-1);  
  Output_functor zonotope_halfspaces_output (generators);

  traverse_combinations<Combination_container, Output_functor>
    (empty_combination, d-2, zonotope_halfspaces_output);
  // we only traverse up to (d-2)-combinations because after that, the
  // output functor takes over and traverses the
  // (d-1)-child-combinations in a manner specific to the halfspace
  // traversal.
  
  return zonotope_halfspaces_output.halfspaces;
}
 
#endif

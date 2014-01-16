#ifndef ZONOTOPE_VOLUME_HPP_
#define ZONOTOPE_VOLUME_HPP_

#include "hyperplane.hpp"
#include "combination_traversal.hpp"
#include "zonotope_volume_output_functor.hpp"
#include "combination_volume_container.hpp"

#include <vector>

/**
 * @brief Generic construction of the set of halfspaces
 *
 * @tparam Combination_container A combination container that also
 *                               maintains the kernel of the
 *                               combination
 */
template <typename NT,
          typename Combination_container = Combination_volume_container<NT> >
NT
zonotope_volume (const std::vector<std::vector<NT> >& generators) {

  typedef Zonotope_volume_output_functor<NT, Combination_container> Output_functor;

  const int d = generators[0].size();

  Combination_container empty_combination (generators, d);
  Output_functor zonotope_volume_output (generators);

  traverse_combinations<Combination_container, Output_functor>
    (empty_combination, d, zonotope_volume_output);
  
  return zonotope_volume_output.volume;
}
 
#endif

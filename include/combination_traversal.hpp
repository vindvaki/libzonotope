#ifndef COMBINATION_TRAVERSAL_HPP_
#define COMBINATION_TRAVERSAL_HPP_

#include <algorithm>
#include <vector>
#include <iostream>

namespace zonotope {

/**
 * @brief A generic depth-first combination traversal algorithm.
 *
 * @tparam Combination_container A container type for combinations. It
 *         must support the operations .size(), .neighbor_upper_bound(),
 *         .extend(int), .get_kernel(), and .is_valid()
 * 
 * @tparam Output_functor A functor type that decides what to do with
 *                        the combinations traversed (it could, for example,
 *                        store them in a set).
 *
 * @param current_combination The root of the current traversal subtree
 *
 * @param max_size A tight upper bound on the allowed size of a combination
 *
 * @param output The output functor.
 */
template <typename Combination_container,
          typename Output_functor>
void traverse_combinations (const Combination_container& current_combination,
                            const int max_size,
                            Output_functor& output)
{
  output( current_combination ); // the current combination is valid
                                 // and the output functor decides
                                 // what to do with it

  if ( current_combination.size() == max_size ) {
    return;
  }

  const int ELEMENT_MAX = current_combination.neighbor_upper_bound();

  for ( int i = 1 + current_combination.back(); i < ELEMENT_MAX; ++i ) {
    Combination_container child_combination( current_combination );
    child_combination.extend(i);
    if ( child_combination.is_valid() ) {
      traverse_combinations(child_combination, max_size, output);
    }
  }
}

} // namespace zonotope

#endif

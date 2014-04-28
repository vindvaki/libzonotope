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
 *                               should implement at least all the
 *                               same interface as Combination_base.
 * 
 * @tparam Output_functor A functor type that decides what to do with
 *                        the combinations traversed (it could, for
 *                        example, store them in a set). Returns true
 *                        iff the current combination a leaf.
 *
 * @param current_combination The root of the current traversal subtree
 *
 * @param max_size A tight upper bound on the allowed size of a combination
 *
 * @param output The output functor.
 */
template <typename Combination_container,
          typename Output_functor>
void traverse_combinations (
  const Combination_container& current_combination,
  Output_functor& output)
{
  if ( output(current_combination) ) {
    // the current combination is a leaf and has been handled
    return;
  }

  for ( int i = current_combination.next_elements_begin();
        i < current_combination.next_elements_end();
        ++i )
  {
    Combination_container child_combination( current_combination );
    child_combination.extend(i);
    if ( child_combination.is_valid() ) {
      traverse_combinations(child_combination, output);
    }
  }
}

} // namespace zonotope

#endif

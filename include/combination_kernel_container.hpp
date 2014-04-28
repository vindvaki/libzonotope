#ifndef COMBINATION_KERNEL_CONTAINER_HPP_
#define COMBINATION_KERNEL_CONTAINER_HPP_

#include "linalg.hpp"
#include "combination_base.hpp"

#include <algorithm>
#include <vector>

namespace zonotope {

/**
 * @brief A combination container for incremental kernel updates
 * 
 */
template <typename NT>
struct Combination_kernel_container : Combination_base
{
  /**
   * The generators of the zonotope
   */
  const std::vector<std::vector<NT> >& generators;

  /**
   * The kernel of generators[combination]
   */
  std::vector<std::vector<NT> > kernel;

  Combination_kernel_container( const std::vector<std::vector<NT> >& generators,
                                const int MAX_SIZE )
    : Combination_base (MAX_SIZE, generators.size())
    , generators (generators)
    , kernel (identity_matrix<NT>(generators[0].size()))
    { }

  void extend(const int i) {
    Combination_base::extend(i);
    update_kernel<NT>(kernel, generators[i]);
  }
  
  /**
   * @brief True iff the combination is independent
   */
  bool is_valid() const {
    return ( (size() + kernel.size()) == generators[0].size() );
  }

};

} // namespace zonotope

#endif

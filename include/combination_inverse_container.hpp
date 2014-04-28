#ifndef COMBINATION_INVERSE_CONTAINER_HPP_
#define COMBINATION_INVERSE_CONTAINER_HPP_

#include "linalg.hpp"
#include "combination_base.hpp"

#include <algorithm>
#include <vector>

namespace zonotope {

/**
 * @brief A combination container for incremental inverse and volume updates
 * 
 */
template <typename NT>
struct Combination_inverse_container : Combination_base
{
  /**
   * The generators of the zonotope
   */
  const std::vector<std::vector<NT> >& generators;

  /**
   * The product of elementary row operation marices to diagonalize
   * transpose(generators[combination])
   */
  std::vector<std::vector<NT> > inverse;

  /**
   * The value of 1 / det(inverse). When k = d,
   * generators[combination] is a nonsingular d-by-d matrix and
   * determinant = det(generators[combination])
   */
  NT determinant;


  Combination_inverse_container( const std::vector<std::vector<NT> >& generators,
                                 const int MAX_SIZE )
    : Combination_base (MAX_SIZE, generators.size())
    , generators (generators)
    , inverse (identity_matrix<NT>(generators[0].size()))
    , determinant (1)
    { }


  void extend(const int i) {
    Combination_base::extend(i);
    update_inverse<NT>(generators, elements, i, inverse, determinant);
  }
  /**
   * @brief True iff the combination is independent
   */
  bool is_valid() const {
    return ( determinant != 0 );
  }

};

} // namespace zonotope

#endif

#ifndef COMBINATION_INVERSE_CONTAINER_HPP_
#define COMBINATION_INVERSE_CONTAINER_HPP_

#include "update_inverse.hpp"
#include "combination_base.hpp"

namespace zonotope {

/**
 * @brief A combination container for incremental inverse and determinant updates
 */

template <typename Zonotope_data_t,
          typename Inverse_number_t_>
struct Combination_inverse_container : Combination_base
{
  typedef Inverse_number_t_ 
    Inverse_number_t;

  typedef Row_matrix<Inverse_number_t, 
                     Zonotope_data_t::Dimension_at_compile_time, 
                     Zonotope_data_t::Dimension_at_compile_time> 
    Matrix_t;

  const Zonotope_data_t& z;

  /**
   * The product of elementary row operation marices to diagonalize
   * transpose(generators[combination])
   */
  Matrix_t inverse;

  /**
   * The value of 1 / det(inverse). When k = d,
   * generators[combination] is a nonsingular d-by-d matrix and
   * determinant = det(generators[combination])
   */
  Inverse_number_t determinant;

  Inverse_number_t scaling;

  Combination_inverse_container( const Zonotope_data_t& z, const int MAX_SIZE )
    : Combination_base (MAX_SIZE, z.num_generators-1)
    , z (z)
    , inverse (Matrix_t::Identity(z.dimension, z.dimension))
    , determinant (1)
    , scaling (1)
    { }

  void extend(const int i) {
    Combination_base::extend(i);
    std::tie(determinant, inverse) = 
      update_inverse(
        this->size()-1,
        z.generators(Inverse_number_t(0)).col(i),
        inverse,
        determinant);

    scaling *= static_cast<Inverse_number_t>(z.scaling(Inverse_number_t(0), i));
  }

  bool is_valid() const {
    return ( determinant != 0 );
  }

};

} // namespace zonotope

#endif

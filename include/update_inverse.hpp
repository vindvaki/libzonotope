#ifndef UPDATE_INVERSE_HPP_
#define UPDATE_INVERSE_HPP_

#include "eigen_utils.hpp"
#include <tuple>

namespace zonotope {

/**
 *
 * @param generators The vector pool indexed by the combination
 *
 * @param combination The current combination of indices
 *
 * @param next_element The element to be appended to the combination
 *                     (i.e. the index of the next vector in the
 *                     correpsonding vector-combination).
 *
 * @param inverse The accumulated row operations to diagonalize
 *                generators[combination].
 *
 * @param determinant `1 / det(inverse)`
 *
 * @post both inverse and determinant have been updated to represent
 *       [combination, next_element].
 *
 */
template <typename Number_t, 
          typename Derived_vector,
          typename Derived_inverse>
std::tuple< Number_t, Derived_matrix<Derived_inverse> >
update_inverse(const int k,
               const typename Eigen::MatrixBase<Derived_vector>& next_column,
               const typename Eigen::MatrixBase<Derived_inverse>& inverse_in,
               const Number_t& determinant_in) 
{

  const int d = inverse_in.rows();

  // init lambda
  Derived_col_vector<Derived_inverse> lambda = inverse_in * next_column;

  // pivot
  int pivot_row = -1;
  for ( int i = k; i < d; ++i ) {
    if ( lambda[i] != 0 ) {
      pivot_row = i;
      break;
    }
  }

  if ( pivot_row == -1 ) {
    // adding x makes the vector combination singular
    return std::make_tuple(Number_t(0), inverse_in);
  }

  // update the inverse
  Derived_matrix<Derived_inverse> inverse_out = inverse_in;

  if ( k != pivot_row ) {
    std::swap( lambda(pivot_row), lambda(k) );
    inverse_out.row(pivot_row).swap(inverse_out.row(k));
  }

  for ( int i = 0; i < d; ++i ) {
    if ( i != k ) {
      inverse_out.row(i) = ((lambda(k) * inverse_out.row(i)) - (lambda(i) * inverse_out.row(k))) / determinant_in;
    }
  }

  // update the determinant
  auto determinant_out = lambda(k);
  if ( k != pivot_row ) {
    determinant_out *= -1;
  }

  return std::make_tuple(determinant_out, inverse_out);
}

} // namezpace zonotope

#endif // UPDATE_INVERSE_HPP_

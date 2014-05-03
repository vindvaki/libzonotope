#ifndef UPDATE_INVERSE_HPP_
#define UPDATE_INVERSE_HPP_

#include "linalg.hpp"

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
template <typename NT, typename Vector_t>
void update_inverse( const std::vector<int>& combination,
                     const Vector_t& x,
                     std::vector<Vector_t>& inverse,
                     NT& determinant) {

  const int k = combination.size();
  const int d = inverse.size();

  // init lambda
  Vector_t lambda (d);
  for ( int i = 0; i < d; ++i ) {
    lambda[i] = 0;
    for ( int j = 0; j < d; ++j ) {
      lambda[i] += inverse[i][j] * x[j];
    }
  }

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
    determinant = 0;
    return;
  }

  if ( k != pivot_row ) {
    std::swap( lambda[pivot_row], lambda[k] );
    std::swap( inverse[pivot_row], inverse[k] );
  }

  // update the inverse
  for ( int i = 0; i < d; ++i ) {
    for ( int j = 0; j < d; ++j ) {
      if ( i != k ) {
        inverse[i][j] *= lambda[k];
        inverse[i][j] -= lambda[i] * inverse[k][j];
        inverse[i][j] /= determinant;
      }
    }
  }

  // update the determinant
  determinant = lambda[k];

  if ( k != pivot_row ) {
    determinant *= -1;
  }
}

} // namezpace zonotope

#endif // UPDATE_INVERSE_HPP_

#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "standardize_vector.hpp"

#include <vector>
#include <utility>
#include <algorithm>

#include <iostream>
#include <utility>

namespace zonotope {

/**
 * @brief The dot product of two vectors.
 */
template <typename NT, typename Vector_t = std::vector<NT> >
inline NT dot( const Vector_t& a, const Vector_t& b )
{
  const int d = a.size();
  NT result = 0;
  for ( int i = 0; i < d; ++i ) {
    result += a[i] * b[i];
  }
  return result;
}

/**
 * @brief Update a kernel basis when a vector is added to a combination
 *
 * @param kernel the basis of the kernel to be restricted

 * @param v the vector to be added to the combination (used to
 *        restrict the kernel).
 *
 * @post kernel has been updated to account for the new vector. If v
 *       is already in the span of the combination, then kernel
 *       remains untouched.
 */

template <typename NT, typename Vector_t = std::vector<NT> >
void update_kernel( std::vector<Vector_t>& kernel, const Vector_t& v )
{
  const int d = kernel[0].size();
  const int k = kernel.size();

  Vector_t x = v;
  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    x[i] = dot<NT, Vector_t> ( kernel[i], v );
    if ( x[i] != 0 ) {
      j = i;
    }
  }

  if ( j == -1 ) {
    // v is already in the combination
    return;
  }

  std::swap(kernel[k-1], kernel[j]);
  std::swap(x[k-1], x[j]);

  for ( int i = 0; i < k - 1; ++i ) {
    for ( int r = 0; r < d; ++r ) {
      kernel[i][r] = x[k-1]*kernel[i][r] - x[i]*kernel[k-1][r];
    }
    standardize_vector<NT, Vector_t> ( kernel[i] );
  }
  kernel.pop_back();
}

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
template <typename NT, typename Vector_t = std::vector<NT> >
void update_inverse( const std::vector<Vector_t>& generators,
                     const std::vector<int>& combination,
                     const int next_element,
                     std::vector<Vector_t>& inverse,
                     NT& determinant) {

  const int k = combination.size();
  const int d = generators[0].size();

  const Vector_t& x = generators[next_element];

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

} // namespace zonotope

#endif // LINALG_HPP_

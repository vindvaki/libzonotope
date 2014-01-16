#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "standardize_vector.hpp"

#include <vector>
#include <utility>
#include <algorithm>

/**
 * @brief The dot product of two vectors.
 */
template <typename NT, typename NT_a = NT, typename NT_b = NT >
inline NT dot( const std::vector<NT_a>& a, const std::vector<NT_b>& b )
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
 * 
 * @return A multiplicative delta for the absolute value of the
 *         determinant.
 */

template <typename NT>
const NT update_kernel( std::vector<std::vector<NT> >& kernel,
                        const std::vector<NT>& v )
{
  using std::vector;
  using std::swap;

  const int d = kernel[0].size();
  const int k = kernel.size();
  
  vector<NT> x (k);
  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    x[i] = dot<NT> ( kernel[i], v );
    if ( x[i] != 0 ) {
      j = i;
    }
  }

  if ( j == -1 ) {
    // v is already in the combination
    return 1;
  }

  swap(kernel[k-1], kernel[j]);
  swap(x[k-1], x[j]);

  NT absolute_determinant_delta = x[k-1];
  
  for ( int i = 0; i < k - 1; ++i ) {
    for ( int r = 0; r < d; ++r ) {
      kernel[i][r] = x[k-1]*kernel[i][r] - x[i]*kernel[k-1][r];
    }
    absolute_determinant_delta *= standardize_vector ( kernel[i] );
  }
  kernel.pop_back();
  
  return absolute_determinant_delta;
}

#endif

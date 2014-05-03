#ifndef UPDATE_KERNEL_HPP_
#define UPDATE_KERNEL_HPP_

#include <vector>

#include "standardize_vector.hpp"
#include "linalg.hpp"

namespace zonotope {

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

template <typename Vector_ZZ>
void update_kernel_zz( std::vector<Vector_ZZ>& kernel, const Vector_ZZ& v )
{
  const int d = kernel[0].size();
  const int k = kernel.size();

  auto x = v;
  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    x[i] = dot( kernel[i], v );
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
    standardize_vector( kernel[i] );
  }
  kernel.pop_back();
}


template <typename Vector_FF>
void update_kernel_ff( std::vector<Vector_FF>& kernel, const Vector_FF& v )
{
  const int d = kernel[0].size();
  const int k = kernel.size();

  auto x = v;
  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    x[i] = dot( kernel[i], v );
    if ( x[i] != 0 ) {
      if ( ( j == -1 ) || ( abs(x[j]) > abs(x[i]) ) ) {
        j = i;
      }
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
      kernel[i][r] -= (x[i] * kernel[k-1][r]) /x[k-1];
    }
  }
  kernel.pop_back();
}



} // namespace zonotope

#endif // UPDATE_KERNEL_HPP_

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
 * @param input_kernel the basis of the kernel to be restricted
 
 * @param v the vector to be added to the combination (used to
 *        restrict the kernel).
 *
 * @return A kernel for the extended combination. If v is already in
 *         the span of the combination, a copy of the input kernel is
 *         returned.
 */

template <typename NT>
std::vector<std::vector<NT> > update_kernel( const std::vector<std::vector<NT> >& input_kernel,
                                             const std::vector<NT>& v )
{
  using std::vector;
  using std::swap;

  const int d = input_kernel[0].size();
  const int k = input_kernel.size();
  vector<vector<NT> > output_kernel = input_kernel;
  
  vector<NT> x (k);
  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    x[i] = dot<NT> ( input_kernel[i], v );
    if ( x[i] != 0 ) {
      j = i;
    }
  }

  if ( j == -1 ) {
    // v is already in the combination
    return output_kernel;
  }

  swap(output_kernel[k-1], output_kernel[j]);
  swap(x[k-1], x[j]);
  
  for ( int i = 0; i < k - 1; ++i ) {
    for ( int r = 0; r < d; ++r ) {
      output_kernel[i][r] = x[k-1]*output_kernel[i][r] - x[i]*output_kernel[k-1][r];
    }
    standardize_vector ( output_kernel[i] );
  }
  output_kernel.pop_back();
  
  return output_kernel;
}

#endif

#ifndef STANDARDIZE_VECTOR_HPP_
#define STANDARDIZE_VECTOR_HPP_

#include <vector>
#include <array>

namespace zonotope {

/**
 * Standardize a vector of integers
 *
 * @tparam Vector_t An iterble, mutable container of integers of type
 *                  Vector_t::value_type. There must exist a function
 *                  
 *                     gcd :: T -> T -> T
 *
 *                  for T = Vector_t::value_type, that returns the
 *                  greatest common divisor.
 * 
 * @param numbers The input vector of integers.
 * 
 * @post The numbers vector has been reduced by the gcd of its entries
 */

template <typename Vector_t>
inline typename Vector_t::value_type
standardize_vector (Vector_t& numbers ) {
  auto a = numbers[0];
  
  for ( const auto& b : numbers ) {
    if ( b == 0 ) {
      continue;
    }
    a = gcd(a, b);
 
    if ( a == 1 ) {
      return a;
    }
  }
  // a = gcd(numbers)
 
  if ( a == 0 ) {
    // NOTE: This is a degenerate case, where all the numbers are 0
    return -1;
  }
  for ( auto& b : numbers ) { 
    b /= a;
  }
  
  // numbers = numbers / a
  return a;
}

} // namespace zonotope

#endif // STANDARDIZE_VECTOR_HPP_

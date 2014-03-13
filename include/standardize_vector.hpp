#ifndef STANDARDIZE_VECTOR_HPP_
#define STANDARDIZE_VECTOR_HPP_

#include <vector>
#include <array>
#include <gmpxx.h>

namespace zonotope {

template <typename Number_t>
inline void _gcd(Number_t& a, const Number_t& b);

template <>
inline void _gcd<mpz_class>(mpz_class& a, const mpz_class& b) {
  mpz_gcd( a.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t() );
}

/**
 * Standardize a vector of integers
 *
 * @param numbers The input vector of integers
 * @post The numbers vector has been reduced by the gcd of its entries
 */

template <typename Number_t, typename Vector_t>
inline Number_t standardize_vector ( Vector_t& numbers ) {
  Number_t a = numbers[0];
  for ( const Number_t& b : numbers ) {
    if ( a == 0 ) {
      if ( b == 0 ) {
        continue;
      }
      a = b;
    }
    _gcd<Number_t> ( a, b );
    // a = gcd(a,b) > 0
    if ( a == 1 ) {
      return a;
    }
  }
  // a = gcd(numbers)
 
  if ( a == 0 ) {
    // NOTE: This is a degenerate case, where all the numbers are 0
    return -1;
  }
  for ( Number_t& b : numbers ) { 
    b /= a;
  }
  // numbers = numbers / a
  return a;
}

/**
 * Standardize a vector of mpq_class rationals.
 *
 * @param numbers The vector of rationals
 * @post The vector has been reduced to the smallest co-directional integral vector.
 */
template<>
inline mpq_class standardize_vector<mpq_class, std::vector<mpq_class> >( std::vector<mpq_class>& numbers ) 
{
  mpz_class a = numbers[0].get_num();
  mpz_class b = numbers[0].get_den();

  for ( const mpq_class& c : numbers ) { 

    mpz_gcd ( a.get_mpz_t(), a.get_mpz_t(), c.get_num_mpz_t() );
    mpz_lcm ( b.get_mpz_t(), b.get_mpz_t(), c.get_den_mpz_t() );

    if ( a == 1 && b == 1 ) {
      return a;
    }
  }
  
  for ( mpq_class& c : numbers ) {
    c.get_num() /= a;
    c.get_num() *= (b / c.get_den());
    c.get_den() = 1;
  }
  // numbers = numbers / (a / b)
  
  return (a / b);
}

} // namespace zonotope

#endif // STANDARDIZE_VECTOR_HPP_

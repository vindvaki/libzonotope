#ifndef STANDARDIZE_VECTOR_HPP_
#define STANDARDIZE_VECTOR_HPP_

#include <vector>
#include <gmpxx.h>

template<typename Number_t>
inline void standardize_vector ( std::vector<Number_t>& );

/**
 * Standardize a vector of mpz_class integers.
 *
 * @param numbers The input vector of integers
 * @post The numbers vector has been reduced by the gcd of its entries
 */
template<>
inline void standardize_vector<mpz_class> ( std::vector<mpz_class>& numbers ) 
{
  mpz_class a = numbers[0];
  for ( const mpz_class& b : numbers ) {
    if ( a == 0 ) {
      if ( b == 0 ) {
        continue;
      }
      a = b;
    }
    mpz_gcd ( a.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t() );
    // a = gcd(a,b) > 0
    if ( a == 1 ) {
      return;
    }
  }
  // a = gcd(numbers)
 
  if ( a == 0 ) {
    // NOTE: This is a degenerate case, where all the numbers are 0
    return;
  }
  for ( mpz_class& b : numbers ) { 
    b /= a;
  }
  // numbers = numbers / a
}

/**
 * Standardize a vector of mpz_class rartionals.
 *
 * @param numbers The vector of rationals
 * @post The vector has been reduced to the smallest co-directional integral vector.
 */
template<>
inline void standardize_vector<mpq_class>( std::vector<mpq_class>& numbers ) 
{
  mpz_class a = numbers[0].get_num();
  mpz_class b = numbers[0].get_den();

  for ( const mpq_class& c : numbers ) { 

    mpz_gcd ( a.get_mpz_t(), a.get_mpz_t(), c.get_num_mpz_t() );
    mpz_lcm ( b.get_mpz_t(), b.get_mpz_t(), c.get_den_mpz_t() );

    if ( a == 1 && b == 1 ) {
      return;
    }
  }

  for ( mpq_class& c : numbers ) {
    c.get_num() /= a;
    c.get_num() *= (b / c.get_den());
    c.get_den() = 1;
  }
}

#endif

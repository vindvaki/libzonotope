#ifndef STANDARDIZE_VECTOR_CGAL_HPP_
#define STANDARDIZE_VECTOR_CGAL_HPP_

#include "standardize_vector.hpp"

#include <CGAL/Gmpz.h>
#include <CGAL/number_utils.h>

/**
 * Same as for mpz_class, but using CGAL::Gmpz
 */

template<>
inline void standardize_vector<CGAL::Gmpz> ( std::vector<CGAL::Gmpz>& numbers ) 
{
  CGAL::Gmpz a = numbers[0];
  for ( const CGAL::Gmpz& b : numbers ) {
    if ( a == 0 ) {
      if ( b == 0 ) {
        continue;
      }
      a = b;
    }
    a = CGAL::gcd( a, b );
    if ( a == 1 ) {
      return;
    }
  }
  // a = gcd(numbers)

  if ( a == 0 ) {
    return;
  }
  for ( CGAL::Gmpz& b : numbers ) {
    b /= a;
  }
  // numbers = numbers / a
}

#endif

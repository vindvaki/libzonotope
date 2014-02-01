#ifndef TYPE_CASTING_FUNCTOR_HPP_
#define TYPE_CASTING_FUNCTOR_HPP_

#include "hyperplane.hpp"

#include <cassert>
#include <gmpxx.h>

template <typename Input_t, typename Output_t>
struct Type_casting_functor {
  Output_t operator() (const Input_t& val) const {
    return Output_t(val);
  }
};

template <>
struct Type_casting_functor<mpz_class, long> {
  long operator() (const mpz_class& val) const {
    assert( val.fits_slong_p() );
    return val.get_si();
  }
};

template <>
struct Type_casting_functor<Hyperplane<mpz_class>, Hyperplane<long> > {
  Type_casting_functor<mpz_class, long> Cast_mpz_to_long;

  Hyperplane<long> operator() (const Hyperplane<mpz_class>& h_mpz ) {
    const int d = h_mpz.normal.size();
      
    Hyperplane<long> h_long (d);
        
    h_long.offset = Cast_mpz_to_long(h_mpz.offset);
    for ( int i = 0; i < d; ++i ) {
      h_long.normal[i] = Cast_mpz_to_long(h_mpz.normal[i]);
    }
    return h_long;
  }
};

#endif // TYPE_CASTING_FUNCTOR_HPP_

#ifndef TYPE_CASTING_FUNCTOR_HPP_
#define TYPE_CASTING_FUNCTOR_HPP_

#include <cassert>
#include <gmpxx.h>

namespace zonotope {

/**
 * Here we define a a utility functor for type-safe conversion between
 * "compatible" types (like integer -> fraction, or integer_type_a ->
 * integer_type_b). Internally, this is (for example) used to convert
 * from arbitrary precision (provided by gmp) to built in fixed
 * precision types. Adding support for other types is as easy as
 * providing a new specialization of the templated interface.
 */


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
struct Type_casting_functor<mpz_class, double> {
  Type_casting_functor<mpz_class, long> cast_to_long;
  double operator() (const mpz_class& val) const {
    return double(cast_to_long(val));
  }
};

template<>
struct Type_casting_functor<std::vector<long>, std::vector<mpz_class> > {
  std::vector<mpz_class> operator() (const std::vector<long>& v_long ) const {
    const int d = v_long.size();
    std::vector<mpz_class> v_mpz (d);
    for ( int i = 0; i < d; ++i ) {
      v_mpz[i] = v_long[i];
    }
    return v_mpz;
  }
};

template<>
struct Type_casting_functor<std::vector<std::vector<long> >, std::vector<std::vector<mpz_class> > > {

  Type_casting_functor<std::vector<long>, std::vector<mpz_class> > cast;

  std::vector<std::vector<mpz_class> > operator() (const std::vector<std::vector<long> >& v_long ) {
    const int d = v_long.size();
    std::vector<std::vector<mpz_class> > v_mpz (d);
    for ( int i = 0; i < d; ++i ) {
      v_mpz[i] = cast(v_long[i]);
    }
    return v_mpz;
  }
};

} // namespace zonotope

#endif // TYPE_CASTING_FUNCTOR_HPP_

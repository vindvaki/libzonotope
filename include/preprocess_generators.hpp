#ifndef PREPROCESS_GENERATORS_HPP_
#define PREPROCESS_GENERATORS_HPP_

#include <gmpxx.h>
#include <vector>

namespace zonotope {

/**
 * Converts the generator matrix from a user facing type to an internal number
 * type. The internal number type should be integral, and the conversion should
 * be of the form
 *
 *     generators_out = cast(scaling_factor * generators_in)
 *
 * with scaling_factor as small as possible.
 *
 */
template <typename User_number_t,
          typename Internal_number_t>
void preprocess_generators (
  const std::vector<std::vector<User_number_t> >& generators_in,
  std::vector<std::vector<Internal_number_t> >& generators_out,
  Internal_number_t& scaling_factor )
{

  const int d = generators_in[0].size();
  const int n = generators_in.size();

  generators_out = std::vector<std::vector<Internal_number_t> > (n, std::vector<Internal_number_t> (d) );
  scaling_factor = 1;

  for ( int k = 0; k < n; ++k ) {
    for ( int i = 0; i < d; ++i ) {
      generators_out[k][i] = generators_in[k][i];
    }
  }
}

template <typename User_number_t,
          typename Internal_number_t>
void preprocess_generators (
  const std::vector<std::vector<User_number_t> >& generators_in,
  std::vector<std::vector<Internal_number_t> >& generators_out)
{
  Internal_number_t scaling_factor;
  preprocess_generators(generators_in, generators_out, scaling_factor);
}

template <>
void
preprocess_generators<mpq_class, mpz_class> (
  const std::vector<std::vector<mpq_class> >& generators_in,
  std::vector<std::vector<mpz_class> >& generators_out,
  mpz_class& scaling_factor)
{
  const int d = generators_in[0].size();
  const int n = generators_in.size();
  generators_out = std::vector<std::vector<mpz_class> > (n, std::vector<mpz_class> (d) );

  mpz_class generators_den_lcm = generators_in[0][0].get_den();
  for ( const auto& v : generators_in ) {
    for ( const mpq_class& x : v ) {
      mpz_lcm ( generators_den_lcm.get_mpz_t(),
                generators_den_lcm.get_mpz_t(),
                x.get_den_mpz_t() );
    }
  }

  for ( int k = 0; k < n; ++k ) {
    for ( int i = 0; i < d; ++i ) {
      const mpq_class& x = generators_in[k][i];
      generators_out[k][i] = (generators_den_lcm * x.get_num()) / x.get_den();
    }
  }

  scaling_factor = generators_den_lcm;
}

template <>
void
preprocess_generators<double, mpz_class> (
    const std::vector<std::vector<double> >& generators_in,
    std::vector<std::vector<mpz_class> >& generators_out,
    mpz_class& scaling_factor ) {

  const int d = generators_in[0].size();
  const int n = generators_in.size();
  std::vector<std::vector<mpq_class> > generators_mpq = std::vector<std::vector<mpq_class> > (n, std::vector<mpq_class> (d));

  for ( int k = 0; k < n; ++k ) {
    for ( int i = 0; i < d; ++i ) {
      generators_mpq[k][i] = generators_in[k][i];
    }
  }
  preprocess_generators<mpq_class, mpz_class> (generators_mpq, generators_out, scaling_factor);
}

} // namespace zonotope

#endif // PREPROCESS_GENERATORS_HPP_

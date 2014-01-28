#ifndef ZONOTOPE_VOLUME_HPP_
#define ZONOTOPE_VOLUME_HPP_

#include "hyperplane.hpp"
#include "combination_traversal.hpp"
#include "zonotope_volume_output_functor.hpp"
#include "combination_inverse_container.hpp"

// External dependencies
#include <vector>
#include <gmpxx.h>



/**
 * @brief Generic construction of the set of halfspaces
 *
 * @tparam Combination_container A combination container that also
 *                               maintains the kernel of the
 *                               combination
 */
template <typename NT = mpz_class,
          typename Combination_container = Combination_inverse_container<NT> >
NT
zonotope_volume (const std::vector<std::vector<NT> >& generators) {

  typedef Zonotope_volume_output_functor<NT, Combination_container> Output_functor;

  const int d = generators[0].size();

  Combination_container empty_combination (generators, d);
  Output_functor zonotope_volume_output (generators);

  traverse_combinations<Combination_container, Output_functor>
    (empty_combination, d, zonotope_volume_output);
  
  return zonotope_volume_output.volume;
}

template <>
long zonotope_volume (const std::vector<std::vector<long> >& generators) {
  const int n = generators.size();
  const int d = generators[0].size();

  std::vector<std::vector<mpz_class> >
    generators_mpz (n, std::vector<mpz_class> (d));

  for ( int k = 0; k < n; ++k ) {
    for ( int i = 0; i < d; ++i ) {
      generators_mpz[k][i] = generators[k][i];
    }
  }
  mpz_class volume_mpz = zonotope_volume<mpz_class> (generators_mpz);
  
  return volume_mpz.get_si();
}
 
#endif








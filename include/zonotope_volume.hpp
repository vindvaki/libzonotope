#ifndef ZONOTOPE_VOLUME_HPP_
#define ZONOTOPE_VOLUME_HPP_

#include "combination_traversal.hpp"
#include "zonotope_volume_output_functor.hpp"
#include "combination_inverse_container.hpp"
#include "preprocess_generators.hpp"
#include "type_casting_functor.hpp"

#include <vector>
#include <gmpxx.h>
#include <cmath>

namespace zonotope {

template <typename T>
T _pow(const T& x, unsigned int n) {}

template <>
mpz_class _pow<mpz_class>(const mpz_class& x, unsigned int n) {
  mpz_class y;
  mpz_pow_ui(y.get_mpz_t(), x.get_mpz_t(), n);
  return y;
}

template <typename User_number_t = mpz_class,
          typename Internal_number_t = mpz_class>
User_number_t zonotope_volume (const std::vector<std::vector<User_number_t> >& generators) {

  typedef Combination_inverse_container<Internal_number_t> Combination_container_t;
  typedef Zonotope_volume_output_functor<Internal_number_t, Combination_container_t> Output_functor_t;
  Type_casting_functor<Internal_number_t, User_number_t> Cast_to_user_type;

  const int d = generators[0].size();

  std::vector<std::vector<Internal_number_t> > internal_generators;
  Internal_number_t scaling_factor;
  preprocess_generators(generators, internal_generators, scaling_factor);

  Combination_container_t empty_combination (internal_generators, d);
  Output_functor_t zonotope_volume_output (internal_generators);

  traverse_combinations<Combination_container_t, Output_functor_t>
    (empty_combination, zonotope_volume_output);

  User_number_t volume = Cast_to_user_type(zonotope_volume_output.volume);
  scaling_factor = _pow<Internal_number_t> (scaling_factor, d);
  volume /= Cast_to_user_type(scaling_factor);
  
  return volume;
}

} // namespace zonotope
 
#endif // ZONOTOPE_VOLUME_HPP_







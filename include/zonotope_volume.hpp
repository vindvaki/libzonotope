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

template<typename Output_number_t,
         typename Internal_number_t,
         typename Zonotope_data_t>
Output_number_t 
zonotope_volume (const Zonotope_data_t& z) 
{
  typedef Combination_inverse_container<Zonotope_data_t, Internal_number_t> 
    Combination_container_t;

  typedef Zonotope_volume_output_functor<Output_number_t, Combination_container_t>
    Output_functor_t;

  Combination_container_t 
    empty_combination (z, z.dimension);

  Output_functor_t 
    zonotope_volume_output (z.dimension);

  traverse_combinations<Combination_container_t, Output_functor_t>
    (empty_combination, zonotope_volume_output);

  Output_number_t volume = static_cast<Output_number_t> (zonotope_volume_output.volume);

  return volume;
}

} // namespace zonotope
 
#endif // ZONOTOPE_VOLUME_HPP_







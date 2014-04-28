#ifndef ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_

#include <vector>
#include <iostream>

#include "output_functor_base.hpp"

namespace zonotope {

template <typename NT, typename Combination_container>
struct Zonotope_volume_output_functor : Output_functor_base<NT>
{
  using typename Output_functor_base<NT>::Generator_container_t;
  
  NT volume;

  Zonotope_volume_output_functor (const Generator_container_t& generators)
    : Output_functor_base<NT> (generators)
    , volume (NT(0))
  { }

  bool operator() (const Combination_container& combination) {
    if ( combination.size() == (this->d) ) {
      volume += abs(combination.determinant);
      return true;
    }
    return false;
  }
};

} // namespace zonotope

#endif // ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_

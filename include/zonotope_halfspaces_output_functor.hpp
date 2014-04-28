#ifndef ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

#include "hyperplane.hpp"
#include "event_point_2.hpp"
#include "type_casting_functor.hpp"
#include "container_output_functor.hpp"
#include "output_functor_base.hpp"

#include <cassert>
#include <vector>

namespace zonotope {

template <typename NT,
          typename Combination_container,
          typename Halfspaces_container_output_functor>
struct Zonotope_halfspaces_output_functor : Output_functor_base<NT>
{
  using typename Output_functor_base<NT>::Generator_container_t;
  
  Halfspaces_container_output_functor& Output_fn;

  Zonotope_halfspaces_output_functor (
    const Generator_container_t& generators,
    Halfspaces_container_output_functor& Output_fn )
    : Output_functor_base<NT>(generators),
      Output_fn( Output_fn )
  {}

  bool operator() (const Combination_container& combination) {

    using std::vector;

    if ( combination.size() == (this->d)-2 ) {

      handle_event_points<NT,
                          vector<NT>,
                          vector<vector<NT> >,
                          Halfspaces_container_output_functor>
        ( combination.back(),
          combination.elements,
          combination.kernel[0],
          combination.kernel[1],
          this->generators,
          Output_fn );
      
      return true;
    }
    return false;
  }
};

} // namespace zonotope

#endif // ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

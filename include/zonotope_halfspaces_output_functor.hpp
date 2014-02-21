#ifndef ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

#include "hyperplane.hpp"
#include "event_point_2.hpp"
#include "type_casting_functor.hpp"
#include "container_output_functor.hpp"

#include <cassert>
#include <vector>

namespace zonotope {

template <typename NT,
          typename Combination_container,
          typename Halfspaces_container_output_functor>
struct Zonotope_halfspaces_output_functor
{

  const std::vector<std::vector<NT> >& generators;

  const int n;
  const int d;

  Halfspaces_container_output_functor& Output_fn;

  Zonotope_halfspaces_output_functor (
    const std::vector<std::vector<NT> >& generators,
    Halfspaces_container_output_functor& Output_fn )
    : generators( generators ),
      n ( generators.size() ),
      d ( generators[0].size() ),
      Output_fn( Output_fn )
  {}

  void operator() (const Combination_container& combination) {

    using std::vector;

    if ( combination.size() == d-2 ) {

      handle_event_points<NT,
                          vector<NT>,
                          vector<vector<NT> >,
                          Halfspaces_container_output_functor>
        ( combination.back(),
          combination.elements,
          combination.kernel[0],
          combination.kernel[1],
          generators,
          Output_fn );
    }
  }
};

} // namespace zonotope

#endif // ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

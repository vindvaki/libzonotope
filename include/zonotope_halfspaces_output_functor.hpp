#ifndef ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

#include "hyperplane.hpp"
#include "event_point_2.hpp"
#include "type_casting_functor.hpp"
#include "container_output_functor.hpp"

#include <cassert>
#include <vector>

template <typename NT,
          typename Combination_container,
          typename Halfspaces_container_output_functor>
struct Zonotope_halfspaces_output_functor
{

  const std::vector<std::vector<NT> >& generators;

  const int n;
  const int d;
  
  Halfspaces_container_output_functor& Output_fn;

  std::vector<NT> generator_sum; 

  Zonotope_halfspaces_output_functor (
    const std::vector<std::vector<NT> >& generators,
    Halfspaces_container_output_functor& Output_fn )
    : generators( generators ),
      n ( generators.size() ),
      d ( generators[0].size() ),
      Output_fn( Output_fn )
  {
    generator_sum = std::vector<NT> ( d, NT ( 0 ) );
    for ( const std::vector<NT>& v : generators ) {
      for ( int i = 0; i < d; ++i ) {
        generator_sum[i] += v[i];
      }
    }
  }

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
          generator_sum,
          Output_fn );
    }
  }
};

#endif

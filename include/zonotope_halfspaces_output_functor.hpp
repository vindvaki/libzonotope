#ifndef ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

#include "hyperplane.hpp"
#include "event_point_2.hpp"
#include "type_casting_functor.hpp"

#include <cassert>
#include <vector>


template <typename Container_t, typename Input_t, typename Output_t = Input_t>
struct Container_output_functor {
  Container_t& data;
  Type_casting_functor<Input_t, Output_t> Cast_type;

  Container_output_functor(Container_t& data) : data(data) {}

  void operator() (const Input_t& val) {
    data.insert( Cast_type(val) );
  }
};

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
      
      // construct a boolean map of the combination for faster
      // membership lookup
      vector<bool> combination_map(n, false);
      for ( int i : combination.elements ) {
        combination_map[i] = true;
      }

      handle_event_points<NT,
                          vector<NT>,
                          vector<vector<NT> >,
                          Halfspaces_container_output_functor,
                          vector<bool> >
        ( combination.back(),
          combination_map,
          combination.kernel[0],
          combination.kernel[1],
          generators,
          generator_sum,
          Output_fn );
    }
  }
};

#endif

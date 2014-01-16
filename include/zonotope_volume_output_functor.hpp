#ifndef ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_

#include "hyperplane.hpp"
#include "event_point_2.hpp"

#include <vector>
#include <set>


template <typename NT, typename Combination_container>
struct Zonotope_volume_output_functor {
  const std::vector<std::vector<NT> >& generators;

  const int n;
  const int d;
  
  NT volume;

  Zonotope_volume_output_functor (const std::vector<std::vector<NT> >& generators)
    : generators( generators ),
      n ( generators.size() ),
      d ( generators[0].size() ),
      volume ( 0 )
  { }

  void operator() (const Combination_container& combination) {
    using std::vector;

    if ( combination.size() == d ) {
      volume += combination.absolute_determinant;
    }
  }
};

#endif

#ifndef ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_

#include <vector>
#include <iostream>

namespace zonotope {

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
      if ( combination.determinant > 0 ) {
        volume += combination.determinant;
      } else {
        volume -= combination.determinant;
      }
    }
  }
};

} // namespace zonotope

#endif

#ifndef ZONOTOPE_HALFSPACES_HPP_
#define ZONOTOPE_HALFSPACES_HPP_

#include "hyperplane.hpp"
#include "combination_traversal.hpp"
#include "combination_container_iml.hpp"
#include "event_point_2.hpp"

#include <vector>
#include <set>

template <typename NT>
struct Zonotope_halfspaces_output_functor {
  const std::vector<std::vector<NT> >& generators;

  const int n;
  const int d;

  std::set<Hyperplane<NT> > halfspaces;

  std::vector<NT> generator_sum;

  Zonotope_halfspaces_output_functor (const std::vector<std::vector<NT> >& generators)
    : generators( generators ),
      n ( generators.size() ),
      d ( generators[0].size() )
  {

    generator_sum = std::vector<NT> ( d, NT ( 0 ) );
    for ( int i = 0; i < n; ++i ) {
      for ( int r = 0; r < d; ++r ) {
        generator_sum[r] += generators[i][r];
      }
    }
  }

  void operator() (const Combination_kernel_container<NT>& combination) {
    using std::vector;
    using std::set;

    if ( combination.size() == d-2 ) {

      vector<bool> combination_map(n, false);
      for ( int i : combination.elements ) {
        combination_map[i] = true;
      }

      handle_event_points<NT,
                          vector<NT>,
                          vector<vector<NT> >,
                          set<Hyperplane<NT> >,
                          vector<bool> >
        ( combination.back(),
          combination_map,
          combination.kernel[0],
          combination.kernel[1],
          generators,
          generator_sum,
          halfspaces );
    }
  }
};

template <typename NT>
std::set<Hyperplane<NT> >
zonotope_halfspaces_iml (const std::vector<std::vector<NT> >& generators) {
  const int d = generators[0].size();
  
  Zonotope_halfspaces_output_functor<NT> zonotope_halfspaces_output (generators);
  Combination_container_iml<NT> empty_combination (generators, d-1);

  traverse_combinations<Combination_kernel_container<NT>,
                        Zonotope_halfspaces_output_functor<NT> >
    (empty_combination, d-2, zonotope_halfspaces_output);
  
  return zonotope_halfspaces_output.halfspaces;
}


#endif

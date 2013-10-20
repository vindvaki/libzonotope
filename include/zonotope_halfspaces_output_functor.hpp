#ifndef ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_HALFSPACES_OUTPUT_FUNCTOR_HPP_

#include "hyperplane.hpp"
#include "event_point_2.hpp"

#include <vector>
#include <set>

//
// TODO: Pass an output functor to handle_event_points instead of a
//       set.  That way, we can (for example) solve some optimization
//       problems without explicitly storing the halfspaces, and
//       output faster when duplicate halfspaces are allowed.
// 


template <typename NT, typename Combination_container>
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

  void operator() (const Combination_container& combination) {
    using std::vector;
    using std::set;

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

#endif

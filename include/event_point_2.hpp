#ifndef EVENT_POINT_2_HPP_
#define EVENT_POINT_2_HPP_

#include "standardize_vector.hpp"
#include "linalg.hpp"
#include "hyperplane.hpp"
#include "compare_by_angle.hpp"

#include <algorithm>
#include <vector>

namespace zonotope {

/**
 * @brief A structure to order vectors in the plane by angle around the origin.
 */
template <typename Number_t>
struct Event_point_2 {

  /**
   * Corresponds to sign*generators[index]
   */
  
  int generator_index;
  int sign;
  
  Number_t x; ///< The planar x-coordinate
  Number_t y; ///< The planar y-coordinate

  Event_point_2 ( int index, int sign, Number_t x, Number_t y ) :
    generator_index ( index ),
    sign ( sign ),
    x ( x ),
    y ( y ) {}

  /**
   *  Compare two event points by angle in `[0, 2*pi)`
   */
  inline bool operator< ( const Event_point_2<Number_t>& other ) const {
    return compare_by_angle<Event_point_2<Number_t> > (*this, other);
  }

};

/**
 * @brief Handle the last step of the zonotope H-rep. construction (the planar view)
 */
template <typename Number_t,
          typename Vector_t,
          typename Generator_container,
          typename Halfspaces_output_functor,
          typename Hyperplane_t = Hyperplane<Number_t> >
inline void handle_event_points (
  const int largest_index,
  const std::vector<int>& current_combination,
  const Vector_t& c0,
  const Vector_t& c1,
  const Generator_container& generators,
  Halfspaces_output_functor& output_fn )
{
  using std::vector;

  const int n = generators.size();
  const int d = generators[0].size();

  // A vector that projects to the inequality offset
  Vector_t offset_vector( d );

  // The event points in the plane spanned by c0, c1
  vector<Event_point_2<Number_t> > event_points;


  // construct a boolean map of the combination for faster
  // membership lookup
  vector<bool> is_elem(n, false);
  for ( int i : current_combination ) {
    is_elem[i] = true;
  }

  // Generate the event points
  for ( int i = 0; i < n; ++i ) {
    if ( is_elem[i] ) {
      // i is in the current combination, so it generators[i] projects to the
      // origin in the plane spanned by c0, c1.
      continue;
    }
    // i is in the complementary combination
    const Vector_t& v = generators[i];
    const Number_t x = dot<Number_t>( c0, v );
    const Number_t y = dot<Number_t>( c1, v );

    if ( x != 0 || y != 0 ) {
      // i corresponds to a nontrivial event
      event_points.push_back( Event_point_2<Number_t> ( i,  1,  x,  y ) );
      event_points.push_back( Event_point_2<Number_t> ( i, -1, -x, -y ) );

      if ( y < 0 || ( y == 0 && x < 0 ) ) {
        // v = generators[i] is below the x-axis
        for ( int r = 0; r < d; ++r ) {
          offset_vector[r] += v[r];
        }
      }
    }
  }

  std::sort ( event_points.begin(), event_points.end() );
  // The event points are sorted in counterclockwise order around the origin,
  // by angle in [0, 2*pi)

  // Rotate a halfplane in counterclockwise order round the origin
  // The initial halfplane is everything below the x-axis, and offset_vector
  // is the sum of those generators.
  for ( const auto& event : event_points ) {
    int i = event.generator_index;
    
    for ( int r = 0; r < d; ++r ) {
      offset_vector[r] += event.sign * generators[i][r];
    }
    if ( i > largest_index ) {
      Hyperplane_t h (d);
      for ( int r = 0; r < d; ++r ) {
        h.normal[r] = -event.y * c0[r] + event.x * c1[r];
      }
      standardize_vector<Number_t, Vector_t> ( h.normal );

      h.offset = -dot<Number_t, Vector_t>(h.normal, offset_vector);
      output_fn(h);
    }
  }
}

} // namespace zonotope

#endif // EVENT_POINT_2_HPP_

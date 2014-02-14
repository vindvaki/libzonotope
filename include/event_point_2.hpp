#ifndef EVENT_POINT_2_HPP_
#define EVENT_POINT_2_HPP_

#include "standardize_vector.hpp"
#include "linalg.hpp"
#include "hyperplane.hpp"

#include <algorithm>
#include <vector>


/**
 * @brief A structure to order vectors in the plane by angle around the origin.
 */
template <typename Number_t>
struct Event_point_2 {

  /**
   * @brief An index to identify the point with a generator.
   *
   * - `i >= 0` corresponds to `generator[i]`
   * - `i < 0` corresponds to `-generator[-i-1]`
   */
  int i;

  Number_t x;                                               ///< The planar x-coordinate
  Number_t y;                                               ///< The planar y-coordinate

  Event_point_2 ( int i, Number_t x, Number_t y ) : i ( i ), x ( x ), y ( y ) {}

  /**
   *  Compare two event points by angle in `[0, 2*pi)`
   */
  inline bool operator< ( const Event_point_2<Number_t>& other ) const {
    if ( y >= 0 && other.y <= 0 ) {
      return true;
    }
    if ( y <= 0 && other.y >= 0 ) {
      return false;
    }

    // both points lie on the same side of the x-axis
    return y * other.x < x * other.y ;
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
  const Vector_t& generator_sum,
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
      event_points.push_back( Event_point_2<Number_t> ( i, x, y ) );
      event_points.push_back( Event_point_2<Number_t> ( -i - 1, -x, -y ) );

      if ( y <= 0 ) {
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
  for ( typename vector<Event_point_2<Number_t> >::const_iterator event_it = event_points.begin();
        event_it != event_points.end();
        ++event_it ) {
    const Event_point_2<Number_t>& event = *event_it;

    if ( event.i < 0 ) {
      // a generator has just left the halfplane
      for ( int r = 0; r < d; ++r ) {
        offset_vector[r] -= generators[-1 - event.i][r];
      }
    } else {
      // a generator has just entered the halfplane
      for ( int r = 0; r < d; ++r ) {
        offset_vector[r] += generators[event.i][r];
      }

      if ( event.i > largest_index ) {
        // (current_combination, event.i) is a new combination, so we must
        // consider its pair of inequalities

        Hyperplane_t h(d);
        Hyperplane_t h_antipode(d);

        for ( int r = 0; r < d; ++r ) {
          h.normal[r] = -event.y * c0[r] + event.x * c1[r] ;
        }
        standardize_vector( h.normal );

        for ( int r = 0; r < d; ++r ) {
          h_antipode.normal[r] = -h.normal[r];
          h.offset -= h.normal[r] * offset_vector[r];
          h_antipode.offset += h.normal[r] * generator_sum[r];
        }
        h_antipode.offset += h.offset;

        h.combination = current_combination;
        h.combination.push_back( event.i );
        h_antipode.combination = h.combination;

        output_fn( h );
        output_fn( h_antipode );
      }
    }
  }
}

#endif

#ifndef EVENT_POINT_2_HPP_
#define EVENT_POINT_2_HPP_

#include "linalg.hpp"
#include "zonotope.hpp"

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
  typedef Number_t value_type;
  
  int generator_index;
  int sign;
  
  value_type x; ///< The planar x-coordinate
  value_type y; ///< The planar y-coordinate

  Event_point_2 ( int index, int sign, const value_type& x, const value_type& y ) :
    generator_index ( index ),
    sign ( sign ),
    x ( x ),
    y ( y ) {}

  /**
   *  Compare two event points by angle in `[0, 2*pi)`
   */
  inline bool operator< ( const Event_point_2<Number_t>& other ) const {
    return compare_by_angle<Event_point_2<value_type> > (*this, other);
  }

};

/**
 * @brief Handle the last step of the zonotope H-rep. construction (the planar view)
 */
template <typename Zonotope_data_t,
          typename Output_functor_t,
          typename Kernel_number_t>
inline void handle_event_points (
  const Zonotope_data_t& z,
  const Combination_kernel_container<Zonotope_data_t, Kernel_number_t>& c,
   Output_functor_t& outfn )
{
  typedef typename Output_functor_t::value_type Hyperplane_t;
  
  typedef Event_point_2<Kernel_number_t> EP;
  
  const int n = z.num_generators;
  const int d = z.dimension;

  // A vector that projects to the inequality offset

  typename Hyperplane_t::Vector_t offset_vector (d);
  
  // The event points in the plane spanned by c.kernel[0], c.kernel[1]
  std::vector<EP> event_points;

  // construct a boolean map of the c for faster
  // membership lookup
  std::vector<bool> is_elem(n, false);
  for ( int i : c.elements ) {
    is_elem[i] = true;
  }

  // Generate the event points
  for ( int i = 0; i < n; ++i ) {
    if ( is_elem[i] ) {
      // i is in the current c, so it generators[i] projects to the
      // origin in the plane spanned by c.kernel[0], c.kernel[1].
      continue;
    }

    // i is in the complementary c
    const auto& v = z.generators(c.kernel[0][0])[i];
    const auto x = dot( c.kernel[0], v );
    const auto y = dot( c.kernel[1], v );

    if ( x != 0 || y != 0 ) {
      // i corresponds to a nontrivial event
      event_points.push_back( EP ( i,  1,  x,  y ) );
      event_points.push_back( EP ( i, -1, -x, -y ) );

      if ( y < 0 || ( y == 0 && x < 0 ) ) {
        // v = generators[i] is below the x-axis
        for ( int r = 0; r < d; ++r ) {
          // TODO: Use same type as Hyperplane_t
          offset_vector[r] += z.generators(offset_vector[r])[i][r];
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
      offset_vector[r] += event.sign * z.generators(offset_vector[r])[i][r];
    }
    if ( i > c.elements.back() ) {
      // TODO: Identify where to 
      
      Hyperplane_t h (d);

      for ( int r = 0; r < d; ++r ) {
        Kernel_number_t tmp = -event.y * c.kernel[0][r] + event.x * c.kernel[1][r];
        h.normal[r] = static_cast<typename Hyperplane_t::Number_t>(tmp);
      }

      h.offset = -dot(h.normal, offset_vector);

      // standardize the representation
      if ( h.offset > 0 ) {
        for ( int r = 0; r < d; ++r ) {
          h.normal[r] /= h.offset;
        }
      } else if ( h.offset < 0 ) {
        for ( int r = 0; r < d; ++r ) {
          h.normal[r] /= -h.offset;
        }
        h.offset = -1;
      }
      
      outfn(h);
    }
  }
}

} // namespace zonotope

#endif // EVENT_POINT_2_HPP_

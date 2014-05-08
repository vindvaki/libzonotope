#ifndef ZONOTOPE_EVENT_POINT_2_HPP_
#define ZONOTOPE_EVENT_POINT_2_HPP_

#include "zonotope.hpp"
#include "standardize_vector.hpp"
#include "compare_by_angle.hpp"

#include <type_traits>
#include <algorithm>
#include <vector>

namespace zonotope {

/**
 * @brief A structure to order vectors in the plane by angle around the origin.
 */
template <typename Number_t_>
struct Event_point_2 {

  /**
   * Corresponds to sign*generators[index]
   */
  typedef Number_t_ Number_t;

  int index;
  int sign;

  Number_t x;
  Number_t y;

  Event_point_2 (int index,
                 int sign,
                 const Number_t& x,
                 const Number_t& y)
    : index(index)
    , sign(sign)
    , x(x)
    , y(y)
    { }

  /**
   *  Compare two event points by angle in `[0, 2*pi)`
   */
  inline bool operator< ( const Event_point_2<Number_t>& other ) const {
    return compare_by_angle(*this, other);
  }

};

template <typename H_in, typename H_out>
inline
typename std::enable_if<std::is_same<H_in, H_out>::value, H_out>::type
convert_hyperplane_types_t (H_in& h_in) {
  auto standardizer = standardize_vector(h_in.data());
  assert( standardizer > 0 );
  return h_in;
}

template <typename H_in, typename H_out>
inline
typename std::enable_if<!(std::is_same<H_in, H_out>::value), H_out>::type
convert_hyperplane_types_t ( H_in& h_in) {
  H_out h_out(h_in.combination());

  const int d = h_in.data().size();

  auto standardizer_in = standardize_vector(h_in.data());
  assert( standardizer_in > 0 );

  for ( int i = 0; i < d; ++i ) {
    h_out.data(i) = static_cast<typename H_out::Number_t>(h_in.data(i));
  }

  auto standardizer_out = standardize_vector(h_out.data());
  assert( standardizer_out > 0 );

  return h_out;
}

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
  typedef typename Output_functor_t::value_type Output_hyperplane_t;
  typedef typename Output_hyperplane_t::Number_t Output_number_t;
  typedef Hyperplane<Zonotope_data_t::Dimension_at_compile_time, Kernel_number_t> Hyperplane_t;
  typedef Event_point_2<Kernel_number_t> EP;

  const int n = z.num_generators;
  const int d = z.dimension;

  // A vector that projects to the inequality offset
  typedef Col_vector<Kernel_number_t, Zonotope_data_t::Dimension_at_compile_time>  Offset_vector_t;
  Offset_vector_t offset_vector (Offset_vector_t::Zero(d,1));

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
    const auto& v = z.generators_uniform(c.kernel(0,0)).col(i);
    const auto x = c.kernel.col(0).dot(v);
    const auto y = c.kernel.col(1).dot(v);

    if ( x != 0 || y != 0 ) {
      // i corresponds to a nontrivial event
      event_points.push_back( EP ( i,  1,  x,  y ) );
      event_points.push_back( EP ( i, -1, -x, -y ) );

      if ( y < 0 || ( y == 0 && x < 0 ) ) {
        // v = generators[i] is below the x-axis
        offset_vector += z.generators_uniform(offset_vector(0)).col(i);
      }
    }
  }

  std::sort ( event_points.begin(), event_points.end() );
  // The event points are sorted in counterclockwise order around the origin,
  // by angle in [0, 2*pi)

  // Rotate a halfplane in counterclockwise order round the origin
  // The initial halfplane is everything below the x-axis, and offset_vector
  // is the sum of those generators.

  auto current_elements (c.elements);

  for ( const auto& event : event_points ) {
    int i = event.index;

    offset_vector += event.sign * z.generators_uniform(offset_vector(0)).col(i);

    if ( i > c.elements.back() ) {
      current_elements.push_back(i);

      Hyperplane_t h (current_elements);
      h.normal = -event.y * c.kernel.col(0) + event.x * c.kernel.col(1);
      h.offset() = -h.normal.dot(offset_vector);

      outfn(convert_hyperplane_types_t<Hyperplane_t, Output_hyperplane_t>(h));

      current_elements.pop_back();
    }
  }
}

} // namespace zonotope

#endif // ZONOTOPE_EVENT_POINT_2_HPP_

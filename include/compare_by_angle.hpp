#ifndef COMPARE_BY_ANGLE_HPP_
#define COMPARE_BY_ANGLE_HPP_

/**
 *  Compare two points by angle in `[0, 2*pi)`
 */
template <typename Point_2>
inline bool compare_by_angle ( const Point_2& a, const Point_2& b ) {

  if ( a.y == 0 ) {
    if ( b.y == 0 ) {
      return a.x >= 0;
    }

    if ( a.x >= 0 ) {
      return true;
    }
    // a.x < 0

    return b.y < 0;
  }
  // a.y != 0

  if ( b.y == 0 ) {
    if ( b.x >= 0 ) {
      return false;
    }
    // b.x < 0
    return a.y > 0;
  }
  // b.y != 0

  if ( a.y > 0 && b.y < 0 ) {
    return true;
  }
  if ( a.y < 0 && b.y > 0 ) {
    return false;
  }
  // both points lie on the same side of the x-axis
  return a.y * b.x < a.x * b.y ;
}

#endif // COMPARE_BY_ANGLE_HPP_

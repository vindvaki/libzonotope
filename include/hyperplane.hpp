#ifndef HYPERPLANE_HPP_
#define HYPERPLANE_HPP_

#include <vector>


//
// TODO: Store a supporting combination along with the hyperplane
//       (canonical if possible).
// 


/**
 * Represents a "hyperplane pair" `(normal, offset)`, where `normal` is
 * a vector and `offset` is a number. It is a utility class with
 * various uses. It could for example, be used to (conceptually)
 * represent the set
 *
 *   { x : dot(normal, x) + offset >= 0 },
 *
 * and if the user enforces some fixed dimension and the condition
 * `(offset==1)`, then we have a one-to-one representation of halfspaces
 * in the given dimension.
 */
template <typename Number_t>
struct Hyperplane {
  
  typedef std::vector<Number_t> Vector_t;
  Number_t offset;
  Vector_t normal;
  
  /**
   * @brief Construct the trivial hyperplane in d dimensions
   */
  Hyperplane( const typename Vector_t::size_type d ) : 
    offset( Number_t(0) ), 
    normal( d, Number_t(0) )
  {}
  
  /**
   * @brief Construct a hyperplane from a given offset and normal
   *
   * __Note:__ It is the responsibility of the user to ensure the
   * hyperplanes are in some standard format, if such is desired.
   * Without a standard format, the built in comparison operators are
   * of limited use.
   *
   * For example:
   *
   * - A standard format for equivalence up to a positive constant
   *   multiple can represent halfspaces with orientation.
   *   
   * - A standard format for equivalence up to any constant multiple
   *   can represent hyperplane equations.
   */
  Hyperplane( const Number_t& offset, const Vector_t& normal ) : 
    offset ( offset ), 
    normal ( normal ) {}
  
  /** 
   * @brief Compare two hyperplanes lexicographically
   */
  bool operator< ( const Hyperplane& other ) const {
    if ( offset < other.offset ) {
      return true;
    }
    
    if ( ( offset == other.offset ) && ( normal < other.normal ) ) {
      return true;
    }
    
    return false;
  }

  /**
   * @brief Compare two hyperplanes for (cartesian) coordinate-wise equality
   *
   * Two hyperplanes are considered equal if and only if they have the
   * same dimensions and the same cartesian coordinates.
   */
  bool operator== ( const Hyperplane& other ) const {
    if ( ( offset == other.offset ) && ( dimension() == other.dimension() ) ) {
      const int d = dimension();
      for ( int i = 0; i < d; ++i ) {
        if ( normal[i] != other.normal[i] ) {
          return false;
        }
      }
      return true;
    }
    return false;
  }
  
  typename Vector_t::size_type dimension() const {
    return normal.size();
  }

};


#endif // HYPERPLANE_HPP_

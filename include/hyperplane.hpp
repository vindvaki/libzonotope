#ifndef HYPERPLANE_HPP_
#define HYPERPLANE_HPP_

#include "type_casting_functor.hpp"

#include <vector>

namespace zonotope {

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
template <typename Vector_t_>
struct Hyperplane {

  typedef Vector_t_ Vector_t;
  typedef typename Vector_t::value_type Number_t;

  typedef typename Vector_t::size_type size_type;
  
  Number_t offset;
  Vector_t normal;

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

  Hyperplane( const size_type d )
    : offset ( 0 )
    , normal ( Vector_t(d) )
    {}

  /**
   * This hyperplane type doesn't store the combination. 
   */
  Hyperplane( const std::vector<int>& combination )
    : Hyperplane(combination.size() + 1)
    { }

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
    return (offset == other.offset) && (normal == other.normal);
  }

  size_type dimension() const {
    return normal.size();
  }

};


template <typename Vector_t_>
struct Hyperplane_and_combination {

  typedef Vector_t_ Vector_t;
  typedef typename Vector_t::value_type Number_t;

  typedef typename Vector_t::size_type size_type;
  
  Number_t offset;
  Vector_t normal;

  std::vector<int> combination;

  Hyperplane_and_combination( const std::vector<int>& combination )
    : offset ( 0 )
    , normal ( Vector_t( combination.size()+1 ) )
    , combination( combination )
    { }

  bool operator< ( const Hyperplane_and_combination& other ) const {
    return combination < other.combination;
  }

  bool operator== ( const Hyperplane_and_combination& other ) const {
    return combination == other.combination;
  }

  size_type dimension() const {
    return normal.size();
  }

};

} // namespace zonotope

#endif // HYPERPLANE_HPP_

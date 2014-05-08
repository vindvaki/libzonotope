#ifndef ZONOTOPE_HYPERPLANE_HPP_
#define ZONOTOPE_HYPERPLANE_HPP_

#include "eigen_utils.hpp"
#include "type_casting_functor.hpp"

#include <vector>
#include <limits>

namespace zonotope {

/**
 * Wraps a vector of d+1 elements to represent a hyperplane in d-space
 */
template <int D, typename Number_t_>
class Hyperplane {

public:
  enum { Data_dimension = (-1 * ( D == -1 )) + ((D+1) * ( D >= 0 )) };
  typedef Number_t_ Number_t;
  typedef Col_vector<Number_t, Data_dimension> Data_vector_t;


private:
  enum { Eigen_needs_to_align = (sizeof(Data_vector_t)%16 == 0) };
  Data_vector_t data_;
  std::vector<int> combination_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(Eigen_needs_to_align)

  Number_t& offset() const {
    return const_cast<Number_t&>(data_(0));
  }

  Eigen::Block<Data_vector_t, D, 1> normal;

  Data_vector_t& data() const {
    return const_cast<Data_vector_t&>(data_);
  }

  Number_t& data(const int i) const {
    return const_cast<Number_t&>(data_(i));
  }

  std::vector<int>& combination() const {
    return const_cast<std::vector<int>&>(combination_);
  }

  int dimension() const {
    return size()-1;
  }

  int size() const {
    return data_.size();
  }

  Hyperplane( const std::vector<int>& combination )
    : data_( Data_vector_t::Zero( combination.size() + 2, 1 ) )
    , combination_( combination )
    , normal( Eigen::Block<Data_vector_t, D, 1> (data_, 1, 0, dimension(), 1) )
  { }

  /**
   * @brief Compare two hyperplanes lexicographically, first by their data and,
   *        in case of equal data, compare them by their underlying combinations.
   */

  bool compare_by_data_first( const Hyperplane& other ) const {
    for ( int i = 0; i < size(); ++i ) {
      // Loop invariant: data(0..i-1) == other.data(0..i-1)

      if ( data(i) < other.data(i) ) {
        return true;
      }
      if ( data(i) > other.data(i) ) {
        return false;
      }
    }
    // data() == other.data()

    return combination() < other.combination();
  }

  bool compare_by_combination_first( const Hyperplane& other ) const {
    if ( combination() < other.combination() ) {
      return true;
    }
    if ( combination() > other.combination() ) {
      return false;
    }
    // combination() == other.combination()

    for ( int i = 0; i < size(); ++i ) {
      // Loop invariant: data(0..i-1) == other.data(0..i-1)

      if ( data(i) < other.data(i) ) {
        return true;
      }
      if ( data(i) > other.data(i) ) {
        return false;
      }
    }
    // data() == other.data()

    return false;
  }

  bool operator< ( const Hyperplane& other ) const {
    return compare_by_data_first(other);
  }

  /**
   * @brief Compare two hyperplanes for (cartesian) coordinate-wise equality
   *
   * Two hyperplanes are considered equal if and only if they have the
   * same dimensions and the same cartesian coordinates.
   */
  bool operator== ( const Hyperplane& other ) const {
    return data == other.data();
  }
};

} // namespace zonotope

#endif // ZONOTOPE_HYPERPLANE_HPP_

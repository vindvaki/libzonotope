#ifndef STANDARDIZE_VECTOR_HPP_
#define STANDARDIZE_VECTOR_HPP_

#include "type_traits_utils.hpp"

#include <eigen3/Eigen/Eigen>

namespace zonotope {

/**
 * Standardize a vector of integers
 *
 * @tparam Vector_t An iterble, mutable container of integers of type
 *                  Vector_t::value_type. There must exist a function
 *
 *                     gcd :: T -> T -> T
 *
 *                  for T = Vector_t::value_type, that returns the
 *                  greatest common divisor.
 *
 * @param numbers The input vector of integers.
 *
 * @post The numbers vector has been reduced by the gcd of its entries
 */


template <class Derived>
enable_if_integer_t<typename Derived::Scalar>
standardize_vector(Eigen::MatrixBase<Derived> const & numbers_const) {
  auto& numbers = const_cast< Eigen::MatrixBase<Derived>& >(numbers_const);
  typename Derived::Scalar a = 0;
  const int d = numbers.size();
  for ( int i = 0; i < d; ++i ) {
    if ( numbers(i) == 0 ) {
      continue;
    }
    a = abs(gcd(a, numbers(i)));
    if ( a == 1 ) {
      return a;
    }
  }

  if ( a != 0 ) {
    numbers /= a;
  }

  // numbers = numbers / a
  return a;
}

template <class Derived>
enable_if_not_integer_t<typename Derived::Scalar>
standardize_vector (Eigen::MatrixBase<Derived> const & numbers_const) {
  auto& numbers = const_cast< Eigen::MatrixBase<Derived>& >(numbers_const);
  auto norm = numbers.template lpNorm<1>();
  // NOTE: We use the L1 norm (sum of absolute values) to preserve rationality
  if ( norm > 0 ) {
    numbers /= norm;
  }
  return norm;
}


} // namespace zonotope

#endif // STANDARDIZE_VECTOR_HPP_

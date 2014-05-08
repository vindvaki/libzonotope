#ifndef ZONOTOPE_UPDATE_KERNEL_HPP_
#define ZONOTOPE_UPDATE_KERNEL_HPP_

#include "type_traits_utils.hpp"
#include "eigen_utils.hpp"
#include "standardize_vector.hpp"

namespace zonotope {

/**
 * @brief Update a kernel basis when a vector is added to a combination
 *
 * @param kernel the basis of the kernel to be restricted

 * @param v the vector to be added to the combination (used to
 *        restrict the kernel).
 *
 * @post kernel has been updated to account for the new vector. If v
 *       is already in the span of the combination, then kernel
 *       remains untouched.
 */


template <typename Derived_kernel_ZZ,
          typename Derived_vector_ZZ>
enable_if_integer_t<typename Derived_kernel_ZZ::Scalar, 
                    Derived_col_matrix<Derived_kernel_ZZ> >
update_kernel(const Eigen::MatrixBase<Derived_kernel_ZZ>& kernel_in,
              const Eigen::MatrixBase<Derived_vector_ZZ>& v )
{
  const int d = kernel_in.rows();
  const int k = kernel_in.cols();

  Derived_col_matrix<Derived_kernel_ZZ>
    kernel_out = kernel_in.leftCols(k-1);
  
  Derived_row_vector<Derived_kernel_ZZ>
    x = kernel_in.transpose() * v;
  
  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    if ( x(i) != 0 ) {
      j = i;
      // don't break because we'd like to avoid a swap if possible
    }
  }

  if ( j == -1 ) {
    // v is already in the combination
    return kernel_in;
  }
  
  if ( j != k-1 ) {
    kernel_out.col(j) = kernel_in.col(k-1);
    std::swap(x(k-1), x(j));
  }
  
  for ( int i = 0; i < k-1; ++i ) {
    kernel_out.col(i) = x(k-1)*kernel_out.col(i) - x(i)*kernel_in.col(j);
    standardize_vector( kernel_out.col(i) );
  }
  
  return kernel_out;
}

template <typename Derived_kernel_FF,
          typename Derived_vector_FF>
enable_if_not_integer_t<typename Derived_kernel_FF::Scalar, 
                        Derived_col_matrix<Derived_kernel_FF> >
update_kernel(const Eigen::MatrixBase<Derived_kernel_FF>& kernel_in,
                 const Eigen::MatrixBase<Derived_vector_FF>& v)
{
  const int d = kernel_in.rows();
  const int k = kernel_in.cols();
  
  Derived_col_matrix<Derived_kernel_FF>
    kernel_out = kernel_in.leftCols(k-1);
  
  Derived_row_vector<Derived_kernel_FF>
    x = kernel_in.transpose() * v;

  int j = -1;
  for ( int i = 0; i < k; ++i ) {
    if ( x(i) != 0 ) {
      if ( ( j == -1 ) || ( abs(x(j)) > abs(x(i)) ) ) {
        j = i;
      }
    }
  }

  if ( j == -1 ) {
    // v is already in the combination
    return kernel_in;
  }

  if ( j != k-1 ) {
    kernel_out.col(j) = kernel_in.col(k-1);
    std::swap(x(k-1), x(j));
  }

  for ( int i = 0; i < k - 1; ++i ) {
    kernel_out.col(i) -= x(i) * kernel_in.col(j) / x(k-1);
  }
  
  return kernel_out;
}



} // namespace zonotope

#endif // ZONOTOPE_UPDATE_KERNEL_HPP_

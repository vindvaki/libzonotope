#ifndef ZONOTOPE_EIGEN_UTILS_
#define ZONOTOPE_EIGEN_UTILS_

#include <eigen3/Eigen/Eigen>

namespace zonotope {

// 
// Primitive shorthands
//
template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
using Col_matrix = Eigen::Matrix<Scalar,
                                 RowsAtCompileTime,
                                 ColsAtCompileTime,
                                 Eigen::ColMajor,
                                 RowsAtCompileTime,
                                 ColsAtCompileTime>;

template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
using Row_matrix = Eigen::Matrix<Scalar,
                                 RowsAtCompileTime,
                                 ColsAtCompileTime,
                                 Eigen::RowMajor,
                                 RowsAtCompileTime,
                                 ColsAtCompileTime>;


template <typename Scalar, int RowsAtCompileTime>
using Col_vector = Col_matrix<Scalar, RowsAtCompileTime, 1>;

template <typename Scalar, int ColsAtCompileTime>
using Row_vector = Row_matrix<Scalar, 1, ColsAtCompileTime>;

//
// Derived shorthands
//

// inherit everything
template <class Derived>
using Derived_matrix = Col_matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>;

// Deriving a column major matrix, we only inherit the number of rows
template <class Derived>
using Derived_col_matrix = Col_matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Eigen::Dynamic>;

template <class Derived>
using Derived_col_vector = Col_vector<typename Derived::Scalar, Derived::RowsAtCompileTime>;

template <class Derived>
using Derived_row_vector = Row_vector<typename Derived::Scalar, Derived::ColsAtCompileTime>;

} // namespace zonotope

#endif

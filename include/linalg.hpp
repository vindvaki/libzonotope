#ifndef LINALG_HPP_
#define LINALG_HPP_

#include <vector>

namespace zonotope {


/**
 * @brief The dot product of two vectors.
 */
template <typename Vector_t>
inline typename Vector_t::value_type
dot(const Vector_t& a, const Vector_t& b)
{
  const int d = a.size();
  typename Vector_t::value_type result = 0;
  for ( int i = 0; i < d; ++i ) {
    result += a[i] * b[i];
  }
  return result;
}

template <typename Vector_t>
std::vector<Vector_t>
init_column_matrix(int cols, const Vector_t& zero_vector)
{
  return std::vector<Vector_t> (cols, zero_vector);
}

template <typename Matrix_t>
Matrix_t
make_identity_matrix(const Matrix_t& other)
{
  Matrix_t result = other;
  for ( typename Matrix_t::size_type col = 0;  col < other.size(); ++col ) {
    for ( typename Matrix_t::value_type::size_type row = 0; row < other[0].size(); ++row ) {
      result[col][row] = 0;
    }
    result[col][col] = 1;
  }
  return result;
}


template <typename T, typename Matrix_t>
Matrix_t
matrix_from_buffer(int d, int n, const T* buffer)
{
  Matrix_t result (n, typename Matrix_t::value_type(d));
  for ( int k = 0; k < n; ++k ) {
    for ( int i = 0; i < d; ++i ) {
      result[k][i] = buffer[k*d + i];
    }
  }
  return result;
}

} // namespace zonotope

#endif // LINALG_HPP_

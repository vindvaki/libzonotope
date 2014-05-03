#ifndef ZONOTOPE_HPP_
#define ZONOTOPE_HPP_

#include "linalg.hpp"
#include "hyperplane.hpp"

#include <tuple>
#include <vector>
#include <cassert>

#include <boost/multiprecision/gmp.hpp>

namespace zonotope {

using boost::multiprecision::gcd;
using boost::multiprecision::lcm;

/**
 * General case: Wraps a fixed size array type
 *
 * The values are always initalized to T(0);
 */
template <int Compile_time_size>
struct Vector_dimension_wrapper {
  template<typename T>
  struct type : std::array<T, Compile_time_size> {

    type() : std::array<T, Compile_time_size>() {
      for ( auto& val : *this ) {
        val = T(0);
      }
    }

    /**
     * Reduces boilerplate and helps us catch size mismatch (albeit
     * only at runtime).
     */
    type(const int declared_runtime_size) : type() {
      assert( declared_runtime_size == Compile_time_size ); 
    }
  };
};

/**
 * Special case: Wraps a dynamically allocated vector
 */
template <>
struct Vector_dimension_wrapper<-1> {
  template<typename T>
  struct type : std::vector<T> {
    type(const int d) : std::vector<T> (d, T(0)) {}
  };
};

template <int D = -1,
          typename Integer_t = boost::multiprecision::mpz_int,
          typename Rational_t = boost::multiprecision::mpq_rational,
          typename Fast_float_t = boost::multiprecision::mpf_float>
class Zonotope_data
{
public:
  template <typename T>
  using Vector_t = typename Vector_dimension_wrapper<D>::type<T>;

  template <typename T>
  using Matrix_t = typename std::vector< Vector_t<T> >;

  typedef Integer_t    ZZ;
  typedef Rational_t   QQ;
  typedef Fast_float_t FF;

  typedef Vector_t<ZZ> Vector_ZZ;
  typedef Vector_t<QQ> Vector_QQ;
  typedef Vector_t<FF> Vector_FF;
  
  typedef Matrix_t<ZZ> Matrix_ZZ;
  typedef Matrix_t<QQ> Matrix_QQ;
  typedef Matrix_t<FF> Matrix_FF;

  const int dimension;
  const int num_generators;

  //
  // constants used for linear algebra
  // 
  const Vector_ZZ zero_vector_zz;     // `dimension`-element zero vector
  const Matrix_ZZ zero_matrix_zz;     // square `dimension`-by-`dimension zero matrix
  const Matrix_ZZ identity_matrix_zz; // square ... identity matrix
  
  const Vector_FF zero_vector_ff;
  const Matrix_FF zero_matrix_ff;
  const Matrix_FF identity_matrix_ff;
  
  //
  // combinatorially equivalent representations
  // 
  const Matrix_QQ generators_qq; // faithful representation of input
  const Matrix_ZZ generators_zz; // generator-wise integral scaling of qq
  const Matrix_FF generators_ff; // fast floating point precision approximation of qq

  //
  // type based lookup
  //
  inline const Matrix_QQ& generators(const QQ& val) const { return generators_qq; }
  inline const Matrix_ZZ& generators(const ZZ& val) const { return generators_zz; }
  inline const Matrix_FF& generators(const FF& val) const { return generators_ff; }


  //
  // Constructors
  // 


  template<typename T>
  Zonotope_data(int d, int n, const T* generators_buffer)
    : Zonotope_data(matrix_from_buffer<T, Matrix_QQ>(d, n, generators_buffer))
  {
    assert( ( d > 0 ) && ( ( D == d ) || ( D == -1 ) ) );
  }

  template<typename T>
  Zonotope_data(int n, const T* generators_buffer)
    : Zonotope_data(matrix_from_buffer<T, Matrix_QQ>(D, n, generators_buffer))
  {
    static_assert(D > 0, "Input dimension must be constant.");
  }
  
private:
  
  Zonotope_data(const Matrix_QQ& input_generators_qq)
    : dimension(input_generators_qq[0].size())
    , num_generators(input_generators_qq.size())
      
      // zz
    , zero_vector_zz(dimension)
    , zero_matrix_zz(init_column_matrix(dimension, zero_vector_zz))
    , identity_matrix_zz(make_identity_matrix(zero_matrix_zz))
      
      // ff
    , zero_vector_ff(dimension)
    , zero_matrix_ff(init_column_matrix(dimension, zero_vector_ff))
    , identity_matrix_ff(make_identity_matrix(zero_matrix_ff))
      
      // generators
    , generators_qq(input_generators_qq)
    , generators_zz(init_generators_zz(generators_qq))
    , generators_ff(init_generators_ff(generators_qq))
  {
    static_assert( (D == -1) || (D > 0), "Invalid input dimension" );
  }

  Matrix_FF init_generators_ff(const Matrix_QQ& generators_qq) {
    const int d = generators_qq[0].size();
    const int n = generators_qq.size();
    
    Matrix_FF generators_ff = init_column_matrix(n, Vector_FF(d));
    for ( int k = 0; k < n; ++k ) {
      for ( int i = 0; i < d; ++i ) {
        generators_ff[k][i] = static_cast<FF>(generators_qq[k][i]);
      }
    }
    return generators_ff;
  }

  /**
   * Scales each generator to Vector_ZZ.
   */
  Matrix_ZZ init_generators_zz(const Matrix_QQ& generators_qq) {
    const int d = generators_qq[0].size();
    const int n = generators_qq.size();
    
    Matrix_ZZ generators_zz = init_column_matrix(generators_qq.size(), Vector_ZZ(d));
    for ( int k = 0; k < n; ++k ) {
      const auto& col_qq = generators_qq[k];
      auto&       col_zz = generators_zz[k];

      ZZ num_gcd = 0;
      ZZ den_lcm = 1;
      for ( const auto& val : generators_qq[k] ) {
        if ( val == 0 ) {
          continue;
        }
        num_gcd = gcd(num_gcd, static_cast<ZZ>(numerator(val)));
        den_lcm = lcm(den_lcm, static_cast<ZZ>(denominator(val)));
      }
      assert( num_gcd != 0 ); // fails iff generators_qq[k] is all zeros
    
      for ( int i = 0; i < d; ++i ) {
        QQ tmp = (col_qq[i] * den_lcm) / num_gcd;
        assert( tmp == numerator(tmp) );
        col_zz[i] = static_cast<ZZ>(numerator(tmp));
      }
    }
    return generators_zz;
  }
};

} // namespace zonotope

#endif // ZONOTOPE_HPP_

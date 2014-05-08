#ifndef ZONOTOPE_HPP_
#define ZONOTOPE_HPP_

#include "eigen_utils.hpp"
#include "hyperplane.hpp"

#include <cassert>

#include <cmath>
#include <boost/multiprecision/gmp.hpp>
#include <eigen3/Eigen/Eigen>

namespace zonotope {

using boost::multiprecision::gcd;
using boost::multiprecision::lcm;

using boost::multiprecision::mpz_int;
using boost::multiprecision::mpq_rational;
using boost::multiprecision::mpf_float;

template <int D = Eigen::Dynamic,
          typename ZZ_ = mpz_int,
          typename QQ_ = mpq_rational,
          typename FF_ = double>
class Zonotope_data
{
public:

  enum { Dimension_at_compile_time = D };

  template <typename T>
  using Vector_t = Col_vector<T, D>;

  template <typename T>
  using Matrix_t = Col_matrix<T, D, Eigen::Dynamic>;

  typedef ZZ_ ZZ;
  typedef QQ_ QQ;
  typedef FF_ FF;

  typedef Vector_t<ZZ> Vector_ZZ;
  typedef Vector_t<QQ> Vector_QQ;
  typedef Vector_t<FF> Vector_FF;

  typedef Matrix_t<ZZ> Matrix_ZZ;
  typedef Matrix_t<QQ> Matrix_QQ;
  typedef Matrix_t<FF> Matrix_FF;

  const int dimension;
  const int num_generators;

private:
  std::vector<ZZ> scaling_;

public:
  //
  // combinatorially equivalent representations
  //

  const Matrix_QQ generators_qq;         // faithful representation of input
  const Matrix_ZZ generators_zz;         // generator-wise integral scaling of qq
  const Matrix_ZZ generators_zz_uniform; // uniform integral scaling of qq
  const Matrix_FF generators_ff;         // fast floating point precision approximation of qq

  const ZZ scaling(const QQ& type_dummy, int i) const {
    return ZZ(1);
  }

  const ZZ scaling(const FF& type_dummy, int i) const {
    return ZZ(1);
  }

  const ZZ& scaling(const ZZ& type_dummy, int i) const {
    return scaling_[i];
  }

  //
  // type based lookup
  //
  inline const Matrix_QQ& generators(const QQ& val) const { return generators_qq; }
  inline const Matrix_ZZ& generators(const ZZ& val) const { return generators_zz; }
  inline const Matrix_FF& generators(const FF& val) const { return generators_ff; }

  inline const Matrix_QQ& generators_uniform(const QQ& val) const { return generators_qq; }
  inline const Matrix_ZZ& generators_uniform(const ZZ& val) const { return generators_zz_uniform; }
  inline const Matrix_FF& generators_uniform(const FF& val) const { return generators_ff; }



  //
  // Constructors
  //

  template<typename T>
  Zonotope_data(const int d, const int n, const T* generators_buffer)
    : Zonotope_data(convert_to_QQ(d, n, generators_buffer))
  {
    assert( ( d > 0 ) && ( ( D == d ) || ( D == -1 ) ) );
  }

  template<typename T>
  Zonotope_data(int n, const T* generators_buffer)
    : Zonotope_data(D, n, generators_buffer)
  {
    static_assert(D > 0, "Input dimension must be constant.");
  }

private:

  template <typename T>
  Matrix_QQ convert_to_QQ(const int d, const int n, const T* buffer) {
    Matrix_QQ m_qq(d, n);
    for ( int k = 0; k < n; ++k ) {
      for ( int i = 0; i < d; ++i ) {
        m_qq(i, k) = static_cast<QQ>(buffer[k*d + i]);
      }
    }
    return m_qq;
  }

  Zonotope_data(const Matrix_QQ& input_generators_qq)
    : dimension(input_generators_qq.rows())
    , num_generators(input_generators_qq.cols())

      // generators
    , scaling_(num_generators, ZZ(1))
    , generators_qq(input_generators_qq)
    , generators_zz(init_generators_zz())
    , generators_zz_uniform(init_generators_zz_uniform())
    , generators_ff(init_generators_ff())
  {
    static_assert( (D == -1) || (D > 0), "Invalid input dimension" );
  }

  Matrix_FF init_generators_ff() {

    const int d = generators_qq.rows();
    const int n = generators_qq.cols();

    Matrix_FF generators_ff (d, n);

    for ( int k = 0; k < n; ++k ) {
      for ( int i = 0; i < d; ++i ) {
        generators_ff(i, k) = static_cast<FF>(generators_qq(i, k));
      }
    }
    return generators_ff;
  }

  /**
   * Scales each generator to Vector_ZZ.
   */
  Matrix_ZZ init_generators_zz() {
    const int d = generators_qq.rows();
    const int n = generators_qq.cols();

    Matrix_ZZ generators_zz (d, n);

    for ( int k = 0; k < n; ++k ) {
      ZZ& den_lcm = scaling_[k];

      for ( int i = 0; i < d; ++i ) {
        const QQ& val = generators_qq(i, k);
        if ( val == 0 ) {
          continue;
        }
        den_lcm = lcm(den_lcm, static_cast<ZZ>(denominator(val)));
      }

      assert (den_lcm != 0);

      for ( int i = 0; i < d; ++i ) {
        QQ tmp = (generators_qq(i, k) * den_lcm);
        assert( tmp == numerator(tmp) );
        generators_zz(i, k) = static_cast<ZZ>(numerator(tmp));
      }
    }
    return generators_zz;
  }

  Matrix_ZZ init_generators_zz_uniform() {
    const int d = generators_qq.rows();
    const int n = generators_qq.cols();

    Matrix_ZZ generators_zz_uniform (d, n);

    ZZ matrix_den_lcm = 1;

    for ( int k = 0; k < n; ++k ) {
      for ( int i = 1; i < d; ++i ) {
        matrix_den_lcm = lcm(denominator(generators_qq(i,k)), matrix_den_lcm);
      }
    }

    for ( int k = 0; k < n; ++k ) {
      for ( int i = 0; i < d; ++i ) {
        auto num = numerator(generators_qq(i,k));
        auto den = denominator(generators_qq(i,k));
        generators_zz_uniform(i, k) = num * (matrix_den_lcm / den);
      }
    }
    return generators_zz_uniform;
  }

};

} // namespace zonotope

#endif // ZONOTOPE_HPP_

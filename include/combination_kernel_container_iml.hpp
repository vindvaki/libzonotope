#ifndef COMBINATION_KERNEL_CONTAINER_IML_HPP_
#define COMBINATION_KERNEL_CONTAINER_IML_HPP_

#include "linalg.hpp"

#include <algorithm>
#include <vector>

#include <gmpxx.h>
#include <iml.h>
#include <cassert>

namespace zonotope {

/**
 * @brief A combination container using IML for kernel computations
 * 
 */
struct Combination_kernel_container_IML {

  /**
   * The generators of the zonotope
   */
  const std::vector<std::vector<mpz_class> >& generators;

  /**
   * An integer in 1..n
   */
  const int MAX_SIZE;
  
  /**
   * A sorted stack of up to MAX_SIZE integers in 0..n-1, that
   * can still be extended to include MAX_SIZE integers.
   */
  std::vector<int> elements;

  /**
   * The kernel of generators[combination]
   */
  std::vector<std::vector<mpz_class> > kernel;
  

  /**
   * @brief Construct an empty combination with room for up to
   *        MAX_SIZE elements.
   */
  Combination_kernel_container_IML( const std::vector<std::vector<mpz_class> >& generators,
                                    const int MAX_SIZE ) :
    generators( generators ),
    MAX_SIZE( MAX_SIZE )
  {
    using std::vector;
    
    const int d = generators[0].size();
    
    kernel = vector<vector<mpz_class> > ( d, vector<mpz_class>( d, 0 ) );
    for ( int i = 0; i < d; ++i ) {
      kernel[i][i] = 1;
    }
  }

  /**
   * @brief Extend the combination to include i
   *
   * @param i an integer in (top() + 1) .. neighbor_upper_bound()
   */
  void extend(const int new_element) {
    elements.push_back(new_element);
    
    const int d = generators[0].size();
    const int k = size();

    mpz_t* combination_vectors = new mpz_t[ k * d ];
    { // initialize combination_vectors
      int j = 0;
      for ( int i : elements ) {
        for ( int r = 0; r < d; ++r ) {
          mpz_init(combination_vectors[j]);
          mpz_set(combination_vectors[j], generators[i][r].get_mpz_t());
          ++j;
        }
      }
    }

    mpz_t *kernel_iml;
    long new_kernel_dimension = nullspaceMP( k, d, combination_vectors, &kernel_iml );

    if ( new_kernel_dimension == kernel.size() ) {
      // the current kernel is still valid
      goto cleanup_iml_extend_combination;
    }

    { // update the stored kernel
      kernel.pop_back();
      int j = 0;
      for ( int i = 0; i < new_kernel_dimension; ++i ) {
        for ( int r = 0; r < d; ++r ) {
          kernel[i][r] = mpz_class( kernel_iml[j] );
          ++j;
        }
      }
    }

  cleanup_iml_extend_combination:
    if ( kernel_iml != NULL ) {
      for ( int j = 0; j < d * new_kernel_dimension; ++j ) {
        mpz_clear( kernel_iml[j] );
      }
      free( kernel_iml );
    }

    for (int j = 0; j < k*d; ++j) {
      mpz_clear(combination_vectors[j]);
    }

    delete [] combination_vectors;
  }

  /**
   * @brief The largest value in the combination, or -1 if empty
   */
  int back() const {
    if (!elements.size()) {
      return -1;
    }
    return elements.back();
  }

  std::vector<int>::size_type size() const {
    return elements.size();
  }

  int neighbor_upper_bound() const {
    const int n = generators.size();
    const int k = elements.size();
    
    return (n - (MAX_SIZE - k) + 1);
  }

  /**
   * @brief true iff i is in the combination
   *
   * Performs an O(log(k)) search for i in a combination of size k.
   */
  bool find(int i) const {
    return std::binary_search(elements.begin(), elements.end(), i);
  }

  /**
   * @brief true iff the combination is independent
   */
  bool is_valid() const {
    const int d = generators[0].size();
    return ( size() + kernel.size() == d );
  }
};

} // namespace zonotope

#endif

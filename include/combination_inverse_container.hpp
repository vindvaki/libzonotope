#ifndef COMBINATION_INVERSE_CONTAINER_HPP_
#define COMBINATION_INVERSE_CONTAINER_HPP_


#include "linalg.hpp"

#include <algorithm>
#include <vector>

/**
 * @brief A combination container for incremental inverse and volume updates
 * 
 */
template <typename NT>
struct Combination_inverse_container {
  /**
   * The generators of the zonotope
   */
  const std::vector<std::vector<NT> >& generators;

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
   * The product of elementary row operation marices to diagonalize
   * transpose(generators[combination])
   */
  std::vector<std::vector<NT> > inverse;

  /**
   * The value of 1 / det(inverse). When k = d,
   * generators[combination] is a nonsingular d-by-d matrix and
   * determinant = det(generators[combination])
   */
  NT determinant;

  /**
   * @brief Construct an empty combination with room for up to
   *        MAX_SIZE elements.
   */
  Combination_inverse_container( const std::vector<std::vector<NT> >& generators,
                                const int MAX_SIZE ) :
    generators( generators ),
    MAX_SIZE( MAX_SIZE ),
    determinant( 1 )
  {
    using std::vector;
    const int d = generators[0].size();
    inverse = vector<vector<NT> > ( d, vector<NT>( d, 0 ) );
    for ( int i = 0; i < d; ++i ) {
      inverse[i][i] = 1;
    }
  }

  /**
   * @brief Extend the combination to include i
   *
   * @param i an integer in (top() + 1) .. neighbor_upper_bound()
   */
  void extend(const int i) {
    update_inverse<NT> ( generators, elements, i, inverse, determinant );
    elements.push_back(i);
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
    return generators.size();
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
    return ( determinant != 0 );
  }

};



#endif

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include <vector>

template <typename NT>
std::vector<std::vector<NT> >
random_generators (int d, int n,
                   const NT& lo, const NT& hi,
                   const NT& seed);

/*
 * 
 * IMPLEMENTATIONS
 *
 */

#include <random>

template<>
std::vector<std::vector<long> >
random_generators (int d, int n,
                   const long& lo, const long& hi,
                   const long& seed) {
  std::vector<std::vector<long> > generators (n, std::vector<long> (d));
  std::uniform_int_distribution<long> dis (lo, hi);
  
  std::mt19937 gen(seed);
  
  for ( int k = 0; k < n; ++k ) {
    for ( int i = 0; i < d; ++i ) {
      generators[k][i] = dis(gen);
    }
  }

  return generators;
}


#endif /* TEST_UTILS_HPP_ */


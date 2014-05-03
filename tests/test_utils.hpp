#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include <vector>

template <typename NT>
std::vector<NT>
random_buffer(typename std::vector<NT>::size_type count,
              const NT& lo,
              const NT& hi,
              const NT& seed);

/*
 * 
 * IMPLEMENTATIONS
 *
 */

#include <random>

template<>
std::vector<long>
random_buffer(typename std::vector<long>::size_type count,
              const long& lo,
              const long& hi,
              const long& seed)
{
  std::uniform_int_distribution<long> dis (lo, hi);
  std::mt19937 gen(seed);
  
  std::vector<long> result (count);  
  for ( auto& val : result ) {
    val = dis(gen);
  }
  return result;
}

#endif /* TEST_UTILS_HPP_ */


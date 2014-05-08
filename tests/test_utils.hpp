#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include <vector>

namespace zonotope {

template <typename NT>
std::vector<NT>
random_buffer(typename std::vector<NT>::size_type count,
              const NT& lo,
              const NT& hi,
              const long& seed);

/*
 *
 * IMPLEMENTATIONS
 *
 */

//
// random_buffer
//

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

template<>
std::vector<double>
random_buffer(typename std::vector<double>::size_type count,
              const double& lo,
              const double& hi,
              const long& seed)
{
  std::uniform_real_distribution<> dis (lo, hi);
  std::mt19937 gen(seed);

  std::vector<double> result (count);
  for ( auto& val : result ) {
    val = dis(gen);
  }
  return result;
}

} // namespace zonotope

#endif /* TEST_UTILS_HPP_ */

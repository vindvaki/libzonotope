#ifndef COMBINATION_KERNEL_CONTAINER_HPP_
#define COMBINATION_KERNEL_CONTAINER_HPP_

#include "update_kernel.hpp"
#include "combination_base.hpp"

#include <algorithm>
#include <vector>

#include <gmpxx.h>

namespace zonotope {

/**
 * @brief A combination container for incremental kernel updates
 * 
 */
template <typename Zonotope_data_t, typename Kernel_number_t>
struct Combination_kernel_container : Combination_base
{ };

template <typename Zonotope_data_t>
struct Combination_kernel_container<Zonotope_data_t, typename Zonotope_data_t::ZZ> : Combination_base
{

  typedef typename Zonotope_data_t::ZZ        Kernel_number_t;
  typedef typename Zonotope_data_t::Matrix_ZZ Matrix_t;
  
  const Zonotope_data_t& z;
  
  Matrix_t kernel;

  Combination_kernel_container( const Zonotope_data_t& z,
                                const int MAX_SIZE )
    : Combination_base (MAX_SIZE, z.num_generators)
    , z(z)
    , kernel (z.identity_matrix_zz)
    { }

  void extend(const int i) {
    Combination_base::extend(i);
    update_kernel_zz(kernel, z.generators_zz[i]);
  }
  
  bool is_valid() const {
    return ( (size() + kernel.size()) == z.dimension );
  }
};

template <typename Zonotope_data_t>
struct Combination_kernel_container<Zonotope_data_t, typename Zonotope_data_t::FF> : Combination_base
{
  typedef typename Zonotope_data_t::Matrix_FF Matrix_t;
  typedef typename Zonotope_data_t::FF        Kernel_number_t;
  
  const Zonotope_data_t& z;

  Matrix_t kernel;

  Combination_kernel_container( const Zonotope_data_t& z,
                                const int MAX_SIZE )
    : Combination_base (MAX_SIZE, z.num_generators)
    , z(z)
    , kernel (z.identity_matrix_ff)
    { }

  void extend(const int i) {
    Combination_base::extend(i);
    update_kernel_ff(kernel, z.generators_ff[i]);
  }
  
  bool is_valid() const {
    return ( (size() + kernel.size()) == z.dimension );
  }
};

} // namespace zonotope

#endif

#ifndef COMBINATION_KERNEL_CONTAINER_HPP_
#define COMBINATION_KERNEL_CONTAINER_HPP_

#include "eigen_utils.hpp"
#include "update_kernel.hpp"
#include "combination_base.hpp"

#include <gmpxx.h>

namespace zonotope {

/**
 * @brief A combination container for incremental kernel updates
 */

template <typename Zonotope_data_t, typename Kernel_number_t_>
struct Combination_kernel_container : Combination_base
{
  typedef Kernel_number_t_ Kernel_number_t;
  typedef Col_matrix<Kernel_number_t, Zonotope_data_t::Dimension_at_compile_time, Eigen::Dynamic> Matrix_t;

  const Zonotope_data_t& z;
  Matrix_t kernel;

  Combination_kernel_container( const Zonotope_data_t& z, const int MAX_SIZE )
    : Combination_base(MAX_SIZE, z.num_generators - 1)
    , z(z)
    , kernel(Matrix_t::Identity(z.dimension, z.dimension))
    { }

  bool is_valid() const {
    return ( (size() + kernel.cols()) == z.dimension );
  }

  void extend(const int i) {
    Combination_base::extend(i);
    kernel = update_kernel(kernel, z.generators(kernel(0,0)).col(i));
  }
};

} // namespace zonotope

#endif

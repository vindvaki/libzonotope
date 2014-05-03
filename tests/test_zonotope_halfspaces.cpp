#include "zonotope_halfspaces.hpp"
#include "container_output_functor.hpp"
#include "test_utils.hpp"

#include <gmpxx.h>
#include <iostream>
#include <set>
#include <list>

int main(int argc, char** argv) {
  using namespace std;
  using namespace zonotope;
  
  const int d = 6; // atol(argv[1]);
  const int n = atol(argv[1]);
  
  vector<long> values = random_buffer(d*n, -100L, 100L, 0L);

  typedef Zonotope_data<d> Zonotope_data_defaults;

  typedef Zonotope_data_defaults::ZZ ZZ;
  typedef double FF;
  
  typedef Zonotope_data<d,
                        Zonotope_data_defaults::ZZ,
                        Zonotope_data_defaults::QQ,
                        FF>
    Zonotope_data_FF;

  Zonotope_data_FF z (n, &values[0]);
  
  typedef Hyperplane<Zonotope_data_FF::Vector_FF> Hyperplane_FF;

  typedef std::set<Hyperplane_FF> Container_t;

  //
  // Compute using FF kernel
  //
  Container_t
    halfspaces_kernel_ff;

  Container_output_functor<Container_t>
    halfspaces_kernel_ff_outfn (halfspaces_kernel_ff); 
  
  zonotope_halfspaces<Zonotope_data_FF,
                      Container_output_functor<Container_t>,
                      FF>  (z, halfspaces_kernel_ff_outfn);
  
  // NOTE: the output functor does not depend on the kernel type!

  //
  // Compute using ZZ kernel
  //
  
  Container_t
    halfspaces_kernel_zz;
  
  Container_output_functor<Container_t>
    halfspaces_kernel_zz_outfn (halfspaces_kernel_zz);
  
  zonotope_halfspaces<Zonotope_data_FF,
                      Container_output_functor<Container_t>,
                      ZZ> (z, halfspaces_kernel_zz_outfn);

  //
  // Compute the average error
  //
  FF error = 0;
 
  return 0;
}

#include "zonotope.hpp"
#include "zonotope_volume.hpp"
#include "test_utils.hpp"

#include <gmpxx.h>
#include <iostream>

int main(int argc, char** argv) {
  using namespace std;
  using namespace zonotope;
  
  assert(argc == 2);

  const int d = 6; 
  const int n = atol(argv[1]);

  typedef mpz_int      ZZ;
  typedef mpq_rational QQ;
  typedef double       FF;

  typedef Zonotope_data<d, ZZ, QQ, FF> Zonotope_data_t;

  vector<double> values = random_buffer(d*n, -1.0, 1.0, 0L);

  Zonotope_data_t z (n, &(values[0]));

  auto volume_dd = zonotope_volume<double, FF> (z);
  auto volume_zz = zonotope_volume<double, ZZ> (z);

  cout << "n = " << n << "\n"
       << "d = " << d << "\n"
       << "volume_dd = " << volume_dd << "\n"
       << "volume_zz = " << volume_zz << "\n"
       << "abs(volume_zz - volume_dd) = " << std::abs(volume_zz - volume_dd) << "\n";
 
  return 0;
}

#include "zonotope_volume.hpp"
#include "type_casting_functor.hpp"
#include "test_utils.hpp"


#include <gmpxx.h>
#include <iostream>


#include <cstdlib>

int main(int argc, char** argv) {
  using namespace std;
  
  const int d = atol(argv[1]);
  const int n = atol(argv[2]);

  zonotope::Type_casting_functor<vector<vector<long> >, vector<vector<mpz_class> > > cast;

  vector<vector<long> > generators_long = random_generators(d, n, -100L, 100L, 0L);
  vector<vector<mpz_class> > generators_mpz = cast(generators_long);

  mpz_class volume = zonotope::zonotope_volume<mpz_class> (generators_mpz);
  long volume_long = zonotope::zonotope_volume (generators_long);
  cout << "n = " << n << "\n"
       << "d = " << d << "\n"
       << "volume_gmpz = " << volume << "\n"
       << "volume_long = " << volume_long << "\n\n";
 
  return 0;
}

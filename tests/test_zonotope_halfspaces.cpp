#include "zonotope_halfspaces.hpp"
#include "type_casting_functor.hpp"
#include "test_utils.hpp"

#include <gmpxx.h>
#include <iostream>


int main(int argc, char** argv) {
  using namespace std;
  
  const int d = atol(argv[1]);
  const int n = atol(argv[2]);

  Type_casting_functor<vector<vector<long> >, vector<vector<mpz_class> > > cast;

  vector<vector<long> > generators_long = random_generators(d, n, -100L, 100L, 0L);
  vector<vector<mpz_class> > generators_mpz = cast(generators_long);

  set<Hyperplane<mpz_class> > halfspaces_mpz;
  zonotope_halfspaces<mpz_class> (generators_mpz, halfspaces_mpz);

  set<Hyperplane<long> > halfspaces_long;
  zonotope_halfspaces<long> (generators_long, halfspaces_long);

  cout << "n=" << n << " "
       << "d=" << d << " "
       << "ieqs=" << halfspaces_mpz.size() << " "
       << "ieqs_long=" << halfspaces_long.size() << "\n";
 
  return 0;
}

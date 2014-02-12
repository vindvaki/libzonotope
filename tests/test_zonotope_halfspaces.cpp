#include "zonotope_halfspaces.hpp"

#include <gmpxx.h>
#include <iostream>

int main() {
  using namespace std;
  
  int n_tests;
  cin >> n_tests;
  for ( int t = 0; t < n_tests; ++t ) {
    int n, d;
    cin >> n >> d;

    vector<vector<mpz_class> >
      generators_mpz ( n, vector<mpz_class> ( d ) );

    vector<vector<long> >
      generators_long ( n, vector<long> ( d ) );

    for ( int k = 0; k < n; ++k ) {
      for ( int i = 0; i < d; ++i ) {
        cin >> generators_long[k][i];
        generators_mpz[k][i] = generators_long[k][i];
      }
    }

    set<Hyperplane<mpz_class> > halfspaces_mpz;
    zonotope_halfspaces<mpz_class> (generators_mpz, halfspaces_mpz);

    set<Hyperplane<long> > halfspaces_long;
    zonotope_halfspaces<long> (generators_long, halfspaces_long);

    cout << "n=" << n << " "
         << "d=" << d << " "
         << "ieqs=" << halfspaces_mpz.size() << " "
         << "ieqs_long=" << halfspaces_long.size() << "\n";
  }
 
  return 0;
}

#include "../include/zonotope_volume.hpp"

#include <gmpxx.h>
#include <iostream>

int main() {
  using namespace std;
  
  typedef mpz_class NT;
  
  int n_tests;

  cin >> n_tests;
  for ( int t = 0; t < n_tests; ++t ) {
    int n, d;
    cin >> n >> d;
    vector<vector<NT> > generators ( n, vector<NT> ( d ) );
    vector<vector<long> > generators_long ( n, vector<long> ( d ) );
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < d; ++j ) {
        cin >> generators_long[i][j];
        generators[i][j] = generators_long[i][j];
      }
    }
    NT volume = zonotope_volume<NT> (generators);
    long volume_long = zonotope_volume (generators_long);
    cout << "n = " << n << "\n"
         << "d = " << d << "\n"
         << "volume_gmpz = " << volume << "\n"
         << "volume_long = " << volume_long << "\n\n";
  }
 
  return 0;
}

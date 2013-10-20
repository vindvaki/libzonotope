#include "zonotope_halfspaces.hpp"

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
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < d; ++j ) {
        cin >> generators[i][j];
      }
    }

    set<Hyperplane<NT> > halfspaces = zonotope_halfspaces<NT> (generators);
    cout << "n=" << n << " "
         << "d=" << d << " "
         << "ieqs=" << halfspaces.size() << "\n";
  }

 
  return 0;
}

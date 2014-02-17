#include "vertex_enum.hpp"

#include <vector>
#include <set>
#include <list>

#include <iostream>

#include <CGAL/Gmpzf.h>

typedef long Input_number_t;
typedef CGAL::Gmpzf Exact_number_t;

int main() {
  using namespace std;
  using namespace zonotope;

  int n_tests;
  cin >> n_tests;
  for ( int t = 0; t < n_tests; ++t ) {
    int n, d;
    cin >> n >> d;
    vector<vector<Input_number_t> > generators (n, vector<Input_number_t> (d, Input_number_t(0)));
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < d; ++j ) {
        cin >> generators[i][j];
      }
    }

    Zonotope_vertex_enumerator<Input_number_t, Exact_number_t> z (generators);

    cout << n << " " << d << " " << z.vertices.size() << "\n";
  }

  return 0;
}

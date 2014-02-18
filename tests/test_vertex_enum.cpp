#include "vertex_enum.hpp"
#include "zonotope_vertex_adjacency_oracle_CGAL.hpp"

#include <vector>
#include <set>
#include <list>
#include <iostream>

#include <CGAL/Gmpzf.h>

int main() {
  using namespace std;
  using namespace zonotope;

  typedef Zonotope_vertex_adjacency_oracle_CGAL<long, CGAL::Gmpzf>
      Adjacency_oracle_t;

  int n_tests;
  cin >> n_tests;
  for ( int t = 0; t < n_tests; ++t ) {
    int n, d;
    cin >> n >> d;
    vector<vector<long> > generators (n, vector<long> (d));
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < d; ++j ) {
        cin >> generators[i][j];
      }
    }

    vector<vector<long> > vertices =
        zonotope_vertices<long, Adjacency_oracle_t > (generators);

    cout << n << " " << d << " " << vertices.size() << "\n";
  }

  return 0;
}

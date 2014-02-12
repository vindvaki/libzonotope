#include "cell_enum.hpp"
#include <vector>
#include <iostream>
#include <CGAL/Gmpz.h>

typedef CGAL::Gmpz Number_t;

using namespace std;

struct Output_counter_functor {
  int count;

  Output_counter_functor() : count(0) {}

  void operator() (const Cell_t& c) {
    ++count;
  }
};


struct Standard_output_functor {

  void operator() (const Cell_t& c) {
    for (int i = 0; i < c.size(); ++i) {
      cout << c[i] << " ";
    }
    cout << "\n";
  }

};


int main() {

  int n_tests;
  cin >> n_tests;
  for ( int t = 0; t < n_tests; ++t ) {
    int n, d;
    cin >> n >> d;
    vector<Number_t> b (n, Number_t(0));
    vector<vector<Number_t> > A (n, vector<Number_t> (d, Number_t(0)));
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < d; ++j ) {
        cin >> A[i][j];
      }
    }

    Standard_output_functor output;
    Output_counter_functor output_counter;
    cout << n << " " << d << " ";

    enumerate_cells_reverse_search<Number_t, Output_counter_functor> (A, b, output_counter);
    cout << output_counter.count;

    cout << "\n\n";
  }

  return 0;
}

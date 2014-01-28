#include "../include/zonotope_volume.hpp"

extern "C" {
  long _zonotope_volume(int d, int n, long* _generators) {

    std::vector<std::vector<long> >
      generators (n, std::vector<long> (d));
    
    for ( int k = 0; k < n; ++k ) {
      for ( int i = 0; i < d; ++i ) {
        generators[k][i] = _generators[i + k*n];
      }
    }
    return zonotope_volume(generators);
  }
}






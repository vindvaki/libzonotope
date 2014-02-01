#include "../include/zonotope_volume.hpp"
#include "../include/zonotope_halfspaces.hpp"


/**
 * Convert the column-major m-by-n matrix `arr` from a flat array to
 * an m-by-n matrix. If `transpose == true`, then instead output the
 * transpose of `arr`. The contents of `mat` are entirely overwritten
 * (along with its dimensions).
 * 
 */
template <typename T>
static void array_to_matrix_2(const int m,
                       const int n,
                       const bool transpose,
                       const T* arr,
                       std::vector<std::vector<T> >& mat) {

  mat = std::vector<std::vector<T> > (m, std::vector<T> (n));

  if ( transpose ) {
    for ( int i = 0; i < m; ++i ) {
      for ( int j = 0; j < n; ++j ) {
        mat[i][j] = arr[j + i*n];
      }
    }
  } else {
    for ( int i = 0; i < m; ++i ) {
      for ( int j = 0; j < n; ++j ) {
        mat[i][j] = arr[i + j*m];
      }
    }
  }
}


/**
 *
 * Exported interface
 * 
 */

extern "C" {
  long zonotope_volume_long(int d, int n, long* generators) {
    std::vector<std::vector<long> > _generators;
    array_to_matrix_2(n, d, true, generators, _generators);
    return zonotope_volume(_generators);
  }

  long* zonotope_halfspaces_long(int d, int n, long* generators) {
    std::vector<std::vector<long> > _generators;
    array_to_matrix_2(n, d, true, generators, _generators);
    std::set<Hyperplane<long> > _halfspaces;
    zonotope_halfspaces(_generators, _halfspaces);
    
    long* halfspaces = new long[(d+1) * _halfspaces.size()];

    std::set<Hyperplane<long> >::size_type count = 0;
    for ( const Hyperplane<long>& h : _halfspaces ) {
      halfspaces[count++] = h.offset;
      for ( int i = 0; i < d; ++i ) {
        halfspaces[count++] = h.normal[i];
      }
    }
    
    return halfspaces;
  }
}

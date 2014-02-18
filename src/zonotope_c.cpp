#include "zonotope_volume.hpp"
#include "zonotope_halfspaces.hpp"

#include "vertex_enum.hpp"
#include "zonotope_vertex_adjacency_oracle_CGAL.hpp"
#include <CGAL/Gmpzf.h>


#include <cstdlib>


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
#include "zonotope_c.h"

long zonotope_volume_long(int d, int n, long* generators) {
  using namespace zonotope;

  std::vector<std::vector<long> > _generators;
  array_to_matrix_2(n, d, true, generators, _generators);
  return zonotope_volume(_generators);
}

long zonotope_halfspaces_long(const int d, const int n, const long* generators, long** halfspaces) {
  using namespace zonotope;
  using std::vector;
  using std::set;

  vector<vector<long> > _generators;
  array_to_matrix_2(n, d, true, generators, _generators);
  set<Hyperplane<long> > _halfspaces;
  zonotope_halfspaces(_generators, _halfspaces);

  (*halfspaces) = (long*)malloc(sizeof(long) * (d+1) * _halfspaces.size());

  set<Hyperplane<long> >::size_type count = 0;
  for ( const Hyperplane<long>& h : _halfspaces ) {
    (*halfspaces)[count++] = h.offset;
    for ( int i = 0; i < d; ++i ) {
      (*halfspaces)[count++] = h.normal[i];
    }
  }

  return _halfspaces.size();
}

long zonotope_vertices_long(const int d, const int n, const long* generators, long** vertices) {
  using namespace zonotope;
  using std::vector;

  vector<vector<long> > _generators;
  array_to_matrix_2(n, d, true, generators, _generators);

  typedef Zonotope_vertex_adjacency_oracle_CGAL<long, CGAL::Gmpzf>
      Adjacency_oracle_t;

  vector<vector<long> > _vertices =
      zonotope_vertices<long, Adjacency_oracle_t> (_generators);

  (*vertices) = (long*)malloc(sizeof(long) * d * _vertices.size());

  vector<vector<long> >::size_type count = 0;
  for ( const auto& v : _vertices ) {
    for ( int i = 0; i < d; ++i ) {
      (*vertices)[count++] = v[i];
    }
  }

  return _vertices.size();
}

} // extern "C"

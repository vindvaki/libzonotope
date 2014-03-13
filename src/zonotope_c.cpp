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
 * Generic pointer functions to reduce boilerplate
 */

template <typename Number_t>
static Number_t
zonotope_volume_ptr(int d, int n, const Number_t* generators)
{
  using namespace zonotope;
  std::vector<std::vector<Number_t> > _generators;
  array_to_matrix_2(n, d, true, generators, _generators);
  return zonotope_volume(_generators);
}

template <typename Number_t>
static long
zonotope_halfspaces_ptr(const int d, const int n, const Number_t* generators, Number_t** halfspaces)
{
  using namespace zonotope;
  using std::vector;
  using std::set;

  vector<vector<Number_t> > _generators;
  array_to_matrix_2(n, d, true, generators, _generators);
  set<Hyperplane<Number_t> > _halfspaces;
  zonotope_halfspaces<Number_t>(_generators, _halfspaces);

  (*halfspaces) = (Number_t*)malloc(sizeof(Number_t) * (d+1) * _halfspaces.size());

  typename set<Hyperplane<Number_t> >::size_type count = 0;
  for ( const Hyperplane<Number_t>& h : _halfspaces ) {
    (*halfspaces)[count++] = h.offset;
    for ( int i = 0; i < d; ++i ) {
      (*halfspaces)[count++] = h.normal[i];
    }
  }

  return _halfspaces.size();
}

template <typename Number_t>
static long
zonotope_vertices_ptr(const int d, const int n, const Number_t* generators, Number_t** vertices)
{
  using namespace zonotope;
  using std::vector;

  vector<vector<Number_t> > _generators;
  array_to_matrix_2(n, d, true, generators, _generators);

  typedef Zonotope_vertex_adjacency_oracle_CGAL<Number_t, CGAL::Gmpzf>
      Adjacency_oracle_t;

  vector<vector<Number_t> > _vertices =
      zonotope_vertices<Number_t, Adjacency_oracle_t> (_generators);

  (*vertices) = (Number_t*)malloc(sizeof(Number_t) * d * _vertices.size());

  typename vector<vector<Number_t> >::size_type count = 0;
  for ( const auto& v : _vertices ) {
    for ( int i = 0; i < d; ++i ) {
      (*vertices)[count++] = v[i];
    }
  }

  return _vertices.size();
}


//
//
// Exported interface
//
//

extern "C" {
#include "zonotope_c.h"

//
// Volume
//
long zonotope_volume_long(int d, int n, const long* generators) {
  return zonotope_volume_ptr<long>(d, n, generators);
}

double zonotope_volume_double(int d, int n, const double* generators) {
  return zonotope_volume_ptr<double>(d, n, generators);
}

//
// Halfspaces
//
long zonotope_halfspaces_long(const int d, const int n, const long* generators, long** halfspaces) {
  return zonotope_halfspaces_ptr<long>(d, n, generators, halfspaces);
}

long zonotope_halfspaces_double(const int d, const int n, const double* generators, double** halfspaces) {
  return zonotope_halfspaces_ptr<double>(d, n, generators, halfspaces);
}

//
// Vertices
//

long zonotope_vertices_long(const int d, const int n, const long* generators, long** vertices) {
  return zonotope_vertices_ptr<long>(d, n, generators, vertices);
}

long zonotope_vertices_double(const int d, const int n, const double* generators, double** vertices) {
  return zonotope_vertices_ptr<double>(d, n, generators, vertices);
}

} // extern "C"

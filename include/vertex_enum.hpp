#ifndef VERTEX_ENUM_HPP_
#define VERTEX_ENUM_HPP_

// inspired by cell_enum.hpp

#include "linalg.hpp"

// STL dependencies
#include <vector>
#include <stack>
#include <utility>
#include <unordered_set>

// For debugging
#include <iostream>

// TODO: Make sure that we have the correct vertices for sum([0,v] : v in V)

namespace zonotope {

template <typename Adjacency_oracle_t>
struct Sign_flip_functor {

  const Adjacency_oracle_t& adj_;

  Sign_flip_functor(const Adjacency_oracle_t& adj) : adj_(adj) {}

  void operator() (typename Adjacency_oracle_t::Vector_t& v,
                   typename Adjacency_oracle_t::Cell_t& sign_vector,
                   int k) const
  {
    for ( int j = 0; j < adj_.d_; ++j ) {
      if ( sign_vector[k] ) {
        v[j] -= adj_.generators_[k].first[j];
        v[j] += adj_.generators_[k].second[j];
      } else {
        v[j] += adj_.generators_[k].first[j];
        v[j] -= adj_.generators_[k].second[j];
      }
    }
    sign_vector[k] = ( ! sign_vector[k] );
  }
};

template <typename Number_t,
          typename Adjacency_oracle_t,
          typename Flip_functor_t = Sign_flip_functor<Adjacency_oracle_t> >
std::vector<std::vector<Number_t> >
zonotope_vertices_dfs (
    std::vector<Number_t>& current_vertex, // the starting vertex
    std::vector<bool>& sign_vector,
    const Flip_functor_t& flip,
    const Adjacency_oracle_t& is_adjacent  ) // the sign vector of the starting vertex
{
  using std::stack;
  using std::vector;
  using std::pair;
  using std::unordered_set;

  vector<vector<Number_t> > vertices;
  stack<pair<int, int> > sign_flip_stack;
  sign_flip_stack.push( pair<int,int> (-1, -1) );

  const int n = is_adjacent.n_;

  unordered_set<vector<bool> > visited_cells;
  visited_cells.insert(sign_vector);
  vertices.push_back(current_vertex);

  // perform a non-recursive DFS to prevent stack overflow
  while ( ! sign_flip_stack.empty() ) {

    pair<int,int>& current_flip = sign_flip_stack.top();

    int current_flip_index = current_flip.first; // the flip that produced the current sign vector
    int& child_flip_index = current_flip.second; // the last explored child flip from this sign vector

    // search for the next not-yet-visited neighbor
    for (++child_flip_index; child_flip_index != n; ++child_flip_index) {

      if ( is_adjacent(sign_vector, child_flip_index) ) {
        // flipping the sign in child_flip_index results in a neighbor

        flip(current_vertex, sign_vector, child_flip_index);

        if ( visited_cells.find(sign_vector) == visited_cells.end() ) {
          // this is the first time we see current_vertex, so we push to the flip to the stack
          vertices.push_back(current_vertex);
          visited_cells.insert(sign_vector);

          sign_flip_stack.push( pair<int,int> (child_flip_index, -1) );
          break;
        }
        // the vertex has already been visited, so we must undo the flip
        flip(current_vertex, sign_vector, child_flip_index);
      }
    }

    if ( child_flip_index == n) {
      // all neighbors of the current sign vector have been visited
      sign_flip_stack.pop();

      if ( current_flip_index != -1 ) {
        // we're not in the root, so we must undo a flip
        flip(current_vertex, sign_vector, current_flip_index);
      }
    }
  }
  return vertices;
}


template <typename Number_t, typename Adjacency_oracle_t >
std::vector<std::vector<Number_t> >
zonotope_vertices (const std::vector<std::vector<Number_t> >& generators)
{
  using std::vector;
  using std::pair;
  typedef vector<Number_t> Vector_t;
  typedef pair<Vector_t, Vector_t> Segment_t;
  
  const int d = generators[0].size();

  //
  // Preprocess  co-directional generators
  // 
  const int m = generators.size();
  vector<bool> handled (m, false);

  vector<Segment_t> _generators; // the internal representation of the generators
  
  for ( int i = 0; i < m; ++i ) {
    // look for codirectional vectors

    if ( handled[i] ) {
      continue;
    }
    // no multiple of generators[i] exists in _generators yet
    
    const auto& u = generators[i];
    
    pair<vector<Number_t>, vector<Number_t> > _u;
    
    _u.first = u;
    _u.second = vector<Number_t> (d, Number_t(0));
    
    int j;
    
    for ( j = i+1; j < m; ++j ) {
      if ( handled[j] ) {
        // generators[j] is already a multiple of some other generator
        continue;
      }
      
      const auto& v = generators[j];
      
      // add only if v is parallel to u
      int r, s;
      for ( s = 0; s < d; ++s ) {
        if ( u[s]*v[s] != 0 ){
          break;
        }
      }
      for ( r = 0; r < d; ++r ) {
        if ( u[s] * v[r] != u[r] * v[s] ) {
          break;
        }
      }
      if ( r == d ) {
        // u and v are parallel
        for ( r = 0; r < d; ++r ) {
          if ( u[s]*v[s] > 0 ) {
            // u and v are codirectional
            _u.first[r] += v[r];
          } else {
            // u and v have opposite directions
            _u.second[r] += v[r];
          }
        }
        handled[j] = true;
      }
    }

    // push only if _u.first is non-zero
    int r;
    for ( r = 0; r < d; ++r ) {
      if (_u.first[r] != 0 ) {
        break;
      }
    }
    if ( r != d ) {
      _generators.push_back(_u);
    }
    handled[i] = true;
  }

  //
  // Perform arrangement-based adjacency traversal
  // 

  const Adjacency_oracle_t is_adjacent (_generators);
  const Sign_flip_functor<Adjacency_oracle_t> flip (is_adjacent);
  
  const int n = _generators.size();

  // Identify an initial vertex and its sign vector
  vector<bool> sign_vector(n);
  vector<Number_t> current_vertex (d, Number_t(0));
  vector<Number_t> random_point   (d, Number_t(1)); // TODO: Pick a "more random" point

  // TODO: Handle the case when random_point lies on too many hyperplanes at once
  for ( int i = 0; i < n; ++i ) {
    sign_vector[i] = ( 0 < dot<Number_t>(_generators[i].first, random_point) );
    for ( int j = 0; j < d; ++j ) {
      if (sign_vector[i]) {
        current_vertex[j] += _generators[i].first[j];
      } else {
        current_vertex[j] += _generators[i].second[j];
      }
    }
  }

  // enumerate the vertices
  return zonotope_vertices_dfs<Number_t, Adjacency_oracle_t>
      (current_vertex, sign_vector, flip, is_adjacent);
}

} // namespace zonotope

#endif // VERTEX_ENUM_HPP_

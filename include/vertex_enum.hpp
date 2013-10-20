#ifndef _VERTEX_ENUM_HPP__
#define _VERTEX_ENUM_HPP__

// inspired by cell_enum.hpp

#include "linalg.hpp"

// CGAL LP solver
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// STL dependencies
#include <vector>
#include <set>
#include <stack>
#include <utility>
#include <unordered_set>

// For debugging
#include <iostream>

// TODO: Make sure that we have the correct vertices for sum([0,v] : v in V)

using std::vector;
using std::set;
using std::stack;
using std::pair;
using std::unordered_set;

template <typename Input_number_t, typename Exact_number_t> 
struct Zonotope_vertex_enumerator {

  typedef vector<vector<Exact_number_t> > VertexContainer;
  
  VertexContainer vertices;
  
  const vector<vector<Input_number_t> > generators_;
  const int n_;
  const int d_;

  Zonotope_vertex_enumerator( const vector<vector<Input_number_t> >& generators ) : 
    generators_(generators),
    n_(generators.size()),
    d_(generators[0].size())
  {

    // Identify an initial vertex and its sign vector
    vector<bool> sign_vector(n_);
    vector<Exact_number_t> current_vertex (d_, Exact_number_t(0));
    vector<Input_number_t> random_point (d_, Input_number_t(1)); // TODO: Pick a "more random" point

    // TODO: Handle the case when random_point lies on too many hyperplanes at once
    for ( int i = 0; i < n_; ++i ) {
      sign_vector[i] = ( 0 < dot<Input_number_t>(generators_[i], random_point) );
      for ( int j = 0; j < d_; ++j ) {
        current_vertex[j] += sign_vector[i]*generators_[i][j];
      }
    }

    // enumerate the vertices
    enumerateVertices(current_vertex, sign_vector);
  }

  bool isAdjacentCell_CGAL(const vector<bool>& sign_vector, const int k) const {

    typedef CGAL::Quadratic_program<Input_number_t> Program;
    typedef CGAL::Quadratic_program_solution<Exact_number_t> Solution;

    Program lp (CGAL::EQUAL, true, 0, false, 0);

    // Set the objective function
    // lp.set_c(k, 1);

    // Set the constraints
    for ( int i = 0; i < n_; ++i ) {
      if ( i != k ) {
        for ( int j = 0; j < d_; ++j ) {
          if ( sign_vector[i] ) {
            lp.set_a(i, j, generators_[i][j]);
          } else {
            lp.set_a(i, j, -generators_[i][j]);
          }
        } 
        // else {
        //   if (sign_vector[k]) {
        //     lp.set_a(k, j, generators_[k][j]);
        //   } else {
        //     lp.set_a(k, j, -generators_[k][j]);
        //   }
        // }
      }
    }
    for ( int j = 0; j < d_; ++j ) {
      if ( sign_vector[k] ) {
        lp.set_b(j, generators_[k][j]);
      } else {
        lp.set_b(j, -generators_[k][j]);
      }
    }

    // Solve the program
    Solution lp_solution = CGAL::solve_nonnegative_linear_program(lp, Exact_number_t());
    return lp_solution.is_infeasible();
    // return ( lp_solution.objective_value() == 1 );
  }

  bool isAdjacentCell (const vector<bool>& sign_vector, const int k) const {
    return isAdjacentCell_CGAL(sign_vector, k);
  }

  void flip(vector<Exact_number_t>& v, vector<bool>& sign_vector, int k) {
    for ( int j = 0; j < d_; ++j ) {
      if ( sign_vector[k] ) {
        v[j] += generators_[k][j];
      } else {
        v[j] -= generators_[k][j];
      }
    }
    sign_vector[k] = ( ! sign_vector[k] );
  }

  void enumerateVertices (
      vector<Exact_number_t>& current_vertex, // the starting vertex
      vector<bool>& sign_vector) // the sign vector of the starting vertex
  {

    stack<pair<int, int> > sign_flip_stack;
    sign_flip_stack.push( pair<int,int> (-1, -1) );

    unordered_set<vector<bool> > visited_cells; 
    visited_cells.insert(sign_vector);
    vertices.push_back(current_vertex);

    // perform a non-recursive DFS to prevent stack overflow
    while ( ! sign_flip_stack.empty() ) { 

      pair<int,int>& current_flip = sign_flip_stack.top();

      int current_flip_index = current_flip.first; // the flip that produced the current sign vector
      int& child_flip_index = current_flip.second; // the last explored child flip from this sign vector

      // search for the next not-yet-visited neighbor
      for (++child_flip_index; child_flip_index != n_; ++child_flip_index) {

        if ( isAdjacentCell(sign_vector, child_flip_index) ) {
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

      if ( child_flip_index == n_) {
        // all neighbors of the current sign vector have been visited
        sign_flip_stack.pop();

        if ( current_flip_index != -1 ) {
          // we're not in the root, so we must undo a flip
          flip(current_vertex, sign_vector, current_flip_index);
        }
      }
    }
  }
};


#endif

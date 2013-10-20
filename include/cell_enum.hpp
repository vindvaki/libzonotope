#ifndef CELL_ENUM_HPP_
#define CELL_ENUM_HPP_

/**
 * An implementation of the oracle functors for reverse search cell
 * enumeration in arrangements. It can also be use for vertex
 * enumeration in zonotopes by providing the appropriate output
 * functor.
 */

#include "reverse_search.hpp"
#include "linalg.hpp"

#include <vector>

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#include <iostream>

typedef std::vector<int> Cell_t;
const Cell_t NONE_CELL (0);

template <typename Number_t> 
struct Arrangement_adjacency_oracle {
  const int n;
  const int d;

  typedef CGAL::Quadratic_program<Number_t> Program;
  typedef CGAL::Quadratic_program_solution<Number_t> Solution;

  const std::vector<std::vector<Number_t> >& A;
  const std::vector<Number_t>& b;

  Arrangement_adjacency_oracle( const std::vector<std::vector<Number_t> >& A,
                                const std::vector<Number_t>& b ) :
    A(A), 
    b(b), 
    n(A.size()),
    d(A[0].size())
  { }

  Cell_t operator() (const Cell_t& c, const int neighbor_counter) const {

    const int k = neighbor_counter - 1;

    Program lp (CGAL::LARGER, false, 0, false, 0);

    // Set the objective function
    for ( int j = 0; j < d; ++j ) {
      lp.set_c(j, c[k]*A[k][j]);
    }
    lp.set_c0(-c[k]*b[k]);

    // Set the constraints
    for (int i = 0; i < n; ++i) {
      // set c[i]*A[i][j] >= c[i]*b[i]
      if ( i != k ) {
        for ( int j = 0; j < d; ++j ) {
          lp.set_a(j, i, c[i]*A[i][j]);
        }
        lp.set_b(k, c[i]*b[i]);
      }
    }

    // Solve the program
    Solution lp_solution = CGAL::solve_linear_program(lp, Number_t());

    if ( lp_solution.is_infeasible() ) {
      return NONE_CELL;
    } 
    if ( (lp_solution.objective_value() < 0) || lp_solution.is_unbounded() ) {
      // there exists a feasible solution with a strictly negative objective value
      Cell_t c_result = c;
      c_result[k] *= -1;
      return c_result;
    }
    // the objective value is >= 0
    return NONE_CELL;
  }
};

template <typename Number_t> 
struct Arrangement_finite_local_search {
  const Arrangement_adjacency_oracle<Number_t> Adj;

  Arrangement_finite_local_search( const Arrangement_adjacency_oracle<Number_t>& Adj ) :
    Adj(Adj) {}

  Cell_t operator() (const Cell_t& c) const {

    // find the smallest k for which Adj(c, k) != NONE_CELL
    for ( int k = 1; k <= Adj.n; ++k ) {
      if ( c[k-1] < 0 ) {
        Cell_t c_next = Adj(c, k);
        if ( c_next != NONE_CELL ) {
          return c_next;
        }
      }
    }
    // there is no such k, so c must be the sink
    return c;
  }
};

template <typename Output_functor> 
struct Translated_output_functor {
  const Cell_t& translation_sign_vector;
  Output_functor& original_output_functor;
  
  Translated_output_functor( Output_functor& original_output_functor, 
                             const Cell_t& translation_sign_vector ) :
    translation_sign_vector(translation_sign_vector), 
    original_output_functor(original_output_functor)
    { }
  
  void operator() (const Cell_t& c) {
    Cell_t c_translated = c;
    for ( int i = 0; i < c.size(); ++i ) {
      c_translated[i] *= translation_sign_vector[i];
    }
    original_output_functor( c_translated );
  }
};

template <typename Number_t, typename Output_functor>
void enumerate_cells_reverse_search (
    const std::vector<std::vector<Number_t> >& A,
    const std::vector<Number_t>& b,
    Output_functor& output) 
{
  const int n = A.size();
  const int d = A[0].size();

  // TODO: a better way to select the root interior point
  std::vector<Number_t> root_interior_point (d, 1); 
  Cell_t root_sign_vector (n, +1); 
  
  // re-orient the hyperplanes to fit root_sign_vector with respect to
  // root_interior_point
  std::vector<std::vector<Number_t> > A_translated = A;
  std::vector<Number_t> b_translated = b;

  Cell_t translation_sign_vector (n, 1);
  for ( int i = 0; i < n; ++i ) {
    if ( dot<Number_t> (root_interior_point, A[i]) <= b[i] ) {
      translation_sign_vector[i] = -1;
      b_translated[i] *= -1;
      for ( int j = 0; j < d; ++j ) {
        A_translated[i][j] *= -1;
      }
    }
  }

  // init functors
  Arrangement_adjacency_oracle<Number_t> Adj (A_translated, b_translated);
  Arrangement_finite_local_search<Number_t> f (Adj);
  Translated_output_functor<Output_functor>
    translated_output (output, translation_sign_vector);

  // perform the reverse search
  reverse_search<
    Arrangement_adjacency_oracle<Number_t>,
    Arrangement_finite_local_search<Number_t>,
    Cell_t, 
    Translated_output_functor<Output_functor> > 
      ( Adj, Adj.n, root_sign_vector, NONE_CELL, f, translated_output );
}

#endif // CELL_ENUM_HPP_

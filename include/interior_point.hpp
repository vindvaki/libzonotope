#ifndef INTERIOR_POINT_HPP_
#define INTERIOR_POINT_HPP_

#include "hyperplane.hpp"
#include <vector>

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

/**
 * @brief A routine to compute an interior point of a polytope given
 *        as an iterator over hyperplanes.
 *
 * @return An interior point, if it exists, and a vector of size 0
 *         otherwise.
 */
template <typename ET, typename Hyperplane_iterator>
std::vector<CGAL::Quotient<ET> >
interior_point(
  const Hyperplane_iterator& begin,
  const Hyperplane_iterator& end ) 
{

  const int d = (*begin).normal.size();
  
  // In case we get bad input, we return an empty vector
  const std::vector<CGAL::Quotient<ET> > NO_INTERIOR_POINT (0); 

  CGAL::Quadratic_program<ET> lp (CGAL::LARGER, false, 0, false, 0);
  // lp will be the program
  //
  //   minimize    -y
  //   subject to  a*x - y >= -b
  //                     y >= 0
  //
  // which has a feasible solution with y > 0 iff the polytope has an
  // interior point.


  //
  // Set the objective function
  // 
  lp.set_c(d, -1); // minimize -y

  //
  // Set the inequalitie coefficients
  //
  int current_ieq = 0;
  for ( Hyperplane_iterator iter = begin; iter != end; ++iter, ++current_ieq ) {
    const Hyperplane<ET>& h = *iter;

    for ( int j = 0; j < d; ++j ) {
      lp.set_a(j, current_ieq, h.normal[j]); // a*x
    }
    lp.set_a(d, current_ieq, -1);            // a*x - y
    lp.set_b(current_ieq, -h.offset);         // a*x - y >= -b
  }
  lp.set_a(d, current_ieq, 1);               // y >= 0

  //
  // Solve the program and extract the solution
  // 
  CGAL::Quadratic_program_solution<ET> lp_solution = CGAL::solve_linear_program(lp, ET());

  if ( lp_solution.is_infeasible() ) {
    // the intersection is empty
    return NO_INTERIOR_POINT;
  } 
  if ( ( lp_solution.objective_value() < 0 ) || lp_solution.is_unbounded() ) {
    // there exists a feasible solution (x, y) with y > 0
    std::vector<CGAL::Quotient<ET> >
      result ( lp_solution.variable_values_begin(),
               lp_solution.variable_values_end() - 1 );
    // result = the d-dimensional x from the augmented interior point (x, y)
    return result;
  }
  // the intersection is a single point, i.e. the objective value is 0
  return NO_INTERIOR_POINT;
}

#endif

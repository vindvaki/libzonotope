#ifndef ZONOTOPE_VERTEX_ADJACENCY_ORACLE_CGAL_HPP_
#define ZONOTOPE_VERTEX_ADJACENCY_ORACLE_CGAL_HPP_

#include <vector>
#include <utility>

// CGAL LP solver
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

namespace zonotope {

template <typename Input_number_t, typename Exact_number_t>
struct Zonotope_vertex_adjacency_oracle_CGAL {

  typedef std::vector<bool> Cell_t;
  typedef std::vector<Input_number_t> Vector_t;
  typedef std::pair<Vector_t, Vector_t> Segment_t;
  typedef std::vector<Segment_t> Generator_container_t;

  const Generator_container_t& generators_;
  const int d_;
  const int n_;

  Zonotope_vertex_adjacency_oracle_CGAL(const Generator_container_t& generators) :
    generators_(generators),
    d_(generators[0].first.size()),
    n_(generators.size()) {}

  bool operator() (const Cell_t& sign_vector, const int k) const {

    typedef CGAL::Quadratic_program<Input_number_t> Program;
    typedef CGAL::Quadratic_program_solution<Exact_number_t> Solution;

    Program lp (CGAL::EQUAL, true, 0, false, 0);

    // Set the constraints
    for ( int i = 0; i < n_; ++i ) {
      if ( i != k ) {
        for ( int j = 0; j < d_; ++j ) {
          if ( sign_vector[i] ) {
            lp.set_a(i, j, generators_[i].first[j]);
          } else {
            lp.set_a(i, j, -generators_[i].first[j]);
          }
        }
      }
    }
    for ( int j = 0; j < d_; ++j ) {
      if ( sign_vector[k] ) {
        lp.set_b(j, generators_[k].first[j]);
      } else {
        lp.set_b(j, -generators_[k].first[j]);
      }
    }

    // Solve the program
    Solution lp_solution = CGAL::solve_nonnegative_linear_program(lp, Exact_number_t());
    return lp_solution.is_infeasible();
  }
};

} // namespace zonotope

#endif // ZONOTOPE_VERTEX_ADJACENCY_ORACLE_CGAL_HPP_

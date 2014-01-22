
libzonotope README
==================

This is a collection of scripts and algorithms for zonotope
manipulations. The most interesting functions are the generic C++
algorithms implemented in the header files in the `include/`
subdirectory of this project.

A zonotope is a convex polytope that is the Minkowski sum of line
segments. We write

    Z(V) := [0,v_0] + ... + [0,v_{n-1}]

for the Minkowski sum of the line segments `[0,v_i]`, `i=0..n-1`,
where `V=[v_0 .. v_{n-1}]` is the `d×n` matrix with columns `v_i`, and
`0` is the origin in `d` dimensions.


Features
--------

In its current state, the project should be considered no more than a 
proof of concept. For more sophisticated vertex enumeration, we recommend
[MinkSum][1]. The main reason someone might use the library right now
is its small size, and the control it gives over the components used
in its algorithms. We currently provide implementations that depend on the
[linear programming solver][3] from [CGAL](http://cgal.org), but
it's fairly easy to adapt the code to a different solver.


The main files are:

- `include/combination_traversal.cpp`: A generic combination traversal
  algorithm for lexicographical traversal of combinations. Currently
  lacking is a proper definition of the concept
  `Combination_container`, though it's already implemented in
  `Combination_kernel_container<NT>`.

- `include/combination_kernel_container.hpp`: A combination container
  that implements incremental kernel updates.

- `include/zonotope_output_functor.hpp`: An example output functor for
  the combination traversal that constructs the H-representation of
  the zontope. The functor accepts all combinations, but forwards only
  the (d-2)-combinations it accepts to `handle_event_points` which
  generates a batch of halfspaces.
  
- `include/event_point_2.hpp`: Implements the function
  `handle_event_points`, which generates a batch of halfspaces in time
  O(n*log(n)).
  
- `zonotope_halfspaces.hpp`: Outputs the H-represenation by combining
  the appropriate output functor and combination container with the
  combination traversal algorithm.

- `include/linalg.hpp`: Implements in particular the function
  `update_kernel` for efficient kernel updates (used in
  `Combination_kernel_container`).
  
- `include/reverse_search.hpp`: Implements a completely generic
  [reverse search algorithm][2] by Avis and Fukuda,
  using functors for adjacency checks and local search.

- `include/cell_enum.hpp`: Implements the functors for reverse search
  to perform cell enumeration in arrangements, and a wrapper function
  that finds a root cell. By providing an appropriate output functor,
  this can be used for vertex enumeration in zonotopes. This algorithm
  is also from the [original reverse search paper by Avis and Fukuda][2].

- `include/vertex_enum.hpp`: Implements vertex enumeration in
  zonotopes using depth-first-search in the dual arrangement with
  manual stack management, for comparison with reverse search (it can
  be concluded that we need to improve our implementation of reverse
  search).
  
[1]: https://sites.google.com/site/christopheweibel/research/minksum
[2]: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.4487
[3]: http://doc.cgal.org/latest/QP_solver/index.html

TODO
----

- Put the library into a common **namespace**!

- More elaborate **testing** (particularly for highly degenerate
  inputs).

- Special functions for two dimensions (i.e. **zonogons**).

- Implement `O(n^{d-1}\log{n})` general position vertex enumeration.

- Forward a user-facing output functor to `handle_event_points` (by 
  passing it through the zonotope output functor calling 
  `handle_event_points`)

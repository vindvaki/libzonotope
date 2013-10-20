#ifndef _REVERSE_SEARCH_HPP__
#define _REVERSE_SEARCH_HPP__

// TODO: Consider using an output stream instead of a normal functor.

template <typename Adjacency_oracle,
         typename Finite_local_search, 
         typename Vertex, 
         typename Output_functor>
inline void reverse_search(
    const Adjacency_oracle& isAdjacent, 
    const int max_degree, 
    const Vertex& sink, 
    const Vertex& NONE,
    const Finite_local_search& f,
    Output_functor& output) 
{
  Vertex current_vertex = sink;
  int neighbor_counter = 0;
  do {
    while ( neighbor_counter < max_degree ) {
      ++neighbor_counter;
      const Vertex next_vertex = isAdjacent(current_vertex, neighbor_counter);
      if ( ( next_vertex != NONE ) && ( f(next_vertex) == current_vertex ) ) {
        // reverse traverse (with respect to f)
        current_vertex = next_vertex;
        neighbor_counter = 0;
      }
    }
    if ( current_vertex != sink ) {
      output(current_vertex);
      // forward traverse (with respect to f)
      const Vertex prev_vertex = current_vertex;
      current_vertex = f(current_vertex);
      neighbor_counter = 0;
      do {
        // restore neighbor_counter
        ++neighbor_counter;
      } while ( isAdjacent(current_vertex, neighbor_counter) != prev_vertex ); 
    }
  } while ( ( current_vertex != sink ) || ( neighbor_counter != max_degree ) );
  output(sink);
}

#endif

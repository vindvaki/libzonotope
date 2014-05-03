#ifndef CONTAINER_OUTPUT_FUNCTOR_HPP_
#define CONTAINER_OUTPUT_FUNCTOR_HPP_

namespace zonotope {

template <typename Container_t>
struct Container_output_functor {
  
  typedef typename Container_t::value_type value_type;
  
  Container_t& data;

  Container_output_functor(Container_t& data) : data(data) {}

  bool operator() (const value_type& val) {
    data.insert(data.end(), val);
    // iterator position data.end() is necessary for lists and vectors,
    // and optional for sets.
    return true;
  }
};

} // namespace zonotope

#endif // CONTAINER_OUTPUT_FUNCTOR_HPP_

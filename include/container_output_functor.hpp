#ifndef CONTAINER_OUTPUT_FUNCTOR_HPP_
#define CONTAINER_OUTPUT_FUNCTOR_HPP_

namespace zonotope {

template <typename Container_t, typename Input_t, typename Output_t = Input_t>
struct Container_output_functor {
  Container_t& data;
  Type_casting_functor<Input_t, Output_t> Cast_type;

  Container_output_functor(Container_t& data) : data(data) {}

  bool operator() (const Input_t& val) {
    data.insert( Cast_type(val) );
    return true;
  }
};

} // namespace zonotope

#endif // CONTAINER_OUTPUT_FUNCTOR_HPP_

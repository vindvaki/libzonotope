#ifndef OUTPUT_FUNCTOR_BASE_HPP_
#define OUTPUT_FUNCTOR_BASE_HPP_

namespace zonotope {

/**
 * A template for a base output functor for things that rely on a
 * const reference to the generators.
 */
template <typename NT>
struct Output_functor_base {

  typedef std::vector<std::vector<NT> > Generator_container_t;
  
  const Generator_container_t& generators;
  const int n;
  const int d;

  Output_functor_base(const Generator_container_t& generators)
    : generators(generators)
    , n(generators.size())
    , d(generators[0].size())
    { }
};

} // namsepace zonotope

#endif // OUTPUT_FUNCTOR_BASE_HPP_

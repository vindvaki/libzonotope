#ifndef ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_
#define ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_

namespace zonotope {

template <typename Number_t, typename Combination_container>
struct Zonotope_volume_output_functor 
{
  typedef typename Combination_container::Inverse_number_t Internal_number_t;

  const int dimension;
  Number_t  volume;
  Zonotope_volume_output_functor (int dimension)
    : dimension (dimension)
    , volume (0)
  { }

  bool operator() (const Combination_container& c) {
    if ( c.size() == dimension ) {
      // c.scaling == 1 iff all input coefficients are integers.
      // Number_t is only allowed to be integral if c.scaling == 1
      Internal_number_t det = c.determinant;
      if ( det < 0 ) {
        det = -det; 
      }
      volume += static_cast<Number_t>(det) / static_cast<Number_t>(c.scaling);
      return true;
    }
    return false;
  }
};

} // namespace zonotope

#endif // ZONOTOPE_VOLUME_OUTPUT_FUNCTOR_HPP_

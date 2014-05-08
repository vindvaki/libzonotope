// libzonotope dependencies
#include "hyperplane.hpp"
#include "zonotope_halfspaces.hpp"
#include "container_output_functor.hpp"

// testing
#include "test_utils.hpp"
#include <iostream>
#include <chrono>

// containers
#include <vector>
#include <set>
#include <algorithm>


// types and constants
const int Dimension = 6;
typedef zonotope::Zonotope_data<Dimension> Zonotope_data_t;
typedef Zonotope_data_t::ZZ ZZ;
typedef Zonotope_data_t::QQ QQ;
typedef Zonotope_data_t::FF FF;
typedef double Output_number_t;

typedef double Input_number_t;
Input_number_t input_lo = -1000;
Input_number_t input_hi = 1000;
long seed = 0;



// we use this to guarantee combination equality for output testing
// the default comparison operator< is not stable across kernels, or
// even output number types
template <typename H>
struct Compare_hyperplanes_by_combination_first {
  bool operator() (const H& a, const H& b) const {
    return a.compare_by_combination_first(b);
  }
};

using zonotope::Hyperplane;

template <typename Hyperplane_t>
using Cmp_fn_t = Compare_hyperplanes_by_combination_first<Hyperplane_t>;

template <typename H>
using Container_t = std::set< H, Cmp_fn_t<H> >;


template <typename Input_number_t, typename Output_number_t, typename Kernel_number_t>
Container_t< Hyperplane<Dimension, Output_number_t> >
test_zonotope_halfspaces(int num_generators) {
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;
  using std::chrono::duration_cast;

  typedef Hyperplane<Dimension, Output_number_t>  Hyperplane_t;
  typedef Container_t<Hyperplane_t > Halfspaces_container_t;
  typedef zonotope::Container_output_functor<Halfspaces_container_t> Container_outfn_t;

  std::vector<Input_number_t> values = zonotope::random_buffer(Dimension*num_generators, input_lo, input_hi, seed);
  Zonotope_data_t z (num_generators, &values[0]);
  Halfspaces_container_t halfspaces_container;
  Container_outfn_t halfspaces_outfn (halfspaces_container);

  auto t_start = high_resolution_clock::now();

  zonotope::zonotope_halfspaces<Zonotope_data_t, Container_outfn_t, Kernel_number_t>  (z, halfspaces_outfn);

  auto t_end = high_resolution_clock::now();

  auto t_delta = duration_cast<milliseconds>(t_end - t_start).count();

  std::cout << "ieqs = " << halfspaces_container.size() << "\n"
            << "time = " << t_delta << " ms\n";

  return halfspaces_container;
}


int main(int argc, char** argv) {
  using namespace zonotope;

  for ( int argi = 1; argi < argc; ++argi ) {

    const int num_generators = atol(argv[argi]);
    assert( num_generators >= Dimension );

    std::cout << "d = " << Dimension << "\n"
              << "n = " << num_generators << "\n";

    std::cout << "\n";


    std::cout << "ZZ -> double\n";
    auto halfspaces_zz_dd = test_zonotope_halfspaces<Input_number_t, double, ZZ> (num_generators);
    std::cout << "\n";

    std::cout << "FF -> double\n";
    auto halfspaces_ff_dd = test_zonotope_halfspaces<Input_number_t, double, FF> (num_generators);
    std::cout << "\n";

    std::cout << "ZZ -> QQ\n";
    auto halfspaces_zz_qq = test_zonotope_halfspaces<Input_number_t, QQ, ZZ> (num_generators);
    std::cout << "\n";

    // if the halfspace containers don't have equal size, then the rest of the
    // tests are nonsensical
    assert( (halfspaces_zz_dd.size() == halfspaces_ff_dd.size()) &&
            (halfspaces_zz_qq.size() == halfspaces_ff_dd.size()) );

    const int num_halfspaces = halfspaces_ff_dd.size();

    auto iter_zz_dd = halfspaces_zz_dd.begin();
    auto iter_ff_dd = halfspaces_ff_dd.begin();
    auto iter_zz_qq = halfspaces_zz_qq.begin();

    bool combination_equality_zz = true;
    bool combination_equality_ff = true;

    // deviation from exact qq outputs
    double error_tot_zz = 0.0;
    double error_max_zz = 0.0;

    double error_tot_ff = 0.0;
    double error_max_ff = 0.0;

    for ( int m = 0; m < num_halfspaces; ++m, ++iter_zz_dd, ++iter_ff_dd, ++iter_zz_qq ) {

      const auto& h_ff_dd = *iter_zz_dd;
      const auto& h_zz_dd = *iter_ff_dd;
      const auto& h_zz_qq = *iter_zz_qq;

      // if we have NaN values in the outputs then the rest of the tests are
      // nonsensical (and there is a greater error to be addressed)
      assert(!std::isnan(h_ff_dd.offset()));
      assert(!std::isnan(h_zz_dd.offset()));

      combination_equality_zz = combination_equality_zz && (h_zz_dd.combination() == h_zz_qq.combination());
      combination_equality_ff = combination_equality_ff && (h_ff_dd.combination() == h_zz_qq.combination());

      Hyperplane<Dimension, double> h_zz_qq_dd (h_zz_qq.combination());
      for ( int i = 0; i < h_zz_qq_dd.size(); ++i ) {
        h_zz_qq_dd.data(i) = static_cast<double> (h_zz_qq.data(i));
      }

      auto error_ff = (h_zz_dd.data() - h_zz_qq_dd.data()).template lpNorm<1>();
      auto error_zz = (h_ff_dd.data() - h_zz_qq_dd.data()).template lpNorm<1>();
      error_tot_ff += error_ff;
      error_tot_zz += error_zz;

      error_max_zz = std::max(error_zz, error_max_zz);
      error_max_ff = std::max(error_ff, error_max_ff);
    }
    std::cout << "combination_equality_zz = " << combination_equality_zz << "\n"
              << "combination_equality_ff = " << combination_equality_ff << "\n"
              << "error_tot_zz = " << error_tot_zz << "\n"
              << "error_tot_ff = " << error_tot_ff << "\n"
              << "error_max_zz = " << error_max_zz << "\n"
              << "error_max_ff = " << error_max_ff << "\n"
              << "\n\n------------\n\n";
  }

  return 0;
}

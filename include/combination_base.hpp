#ifndef COMBINATION_BASE_HPP_
#define COMBINATION_BASE_HPP_

#include <algorithm>
#include <vector>

namespace zonotope {

struct Combination_base {
  typedef std::vector<int> Element_container;
  typedef Element_container::size_type size_type;

  size_type MAX_SIZE;
  const int MAX_ELEMENT;

  /**
   * A sorted stack of up to MAX_SIZE integers in 0..MAX_ELEMENT-1, that
   * can still be extended to include MAX_SIZE integers.
   */
  Element_container elements;

  /**
   * @brief Construct an empty combination of integers in
   *        0..(MAX_ELEMENT-1) with room for up to MAX_SIZE elements.
   */
  Combination_base( const size_type MAX_SIZE, const int MAX_ELEMENT )
    : MAX_SIZE(MAX_SIZE)
    , MAX_ELEMENT(MAX_ELEMENT)
    { }

  /**
   * @brief The largest value in the combination, or -1 if empty
   */
  int back() const {
    if (!elements.size()) {
      return -1;
    }
    return elements.back();
  }

  size_type size() const {
    return elements.size();
  }

  virtual void extend(const int i) {
    elements.push_back(i);
  }

  /**
   * @brief true iff the combination is "valid".
   *
   * For the base class, we only need the combination to have no more
   * than MAX_SIZE elements. Ever overloaded insetance must enforce
   * this in particular, but they are also free to restrict their
   * definition of "valid".
   */
  virtual bool is_valid() const {
    return size() <= MAX_SIZE;
  }

  /**
   * @brief true iff i is in the combination
   *
   * Performs an O(log(k)) search for i in a combination of size k.
   */
  bool find(int i) const {
    return std::binary_search(elements.begin(), elements.end(), i);
  }

  int next_elements_begin() const {
    return back()+1;
  }

  int next_elements_end() const {
    // At most we can have MAX_SIZE elements.  No element can be larger than
    // MAX_ELEMENTS.  (MAX_SIZE - size()) is the number of elements we must
    // still add.  Thus, MAX_ELEMENT - (MAX_SIZE - size()) is the largest
    // element that we can add while still being able to reach MAX_SIZE
    return MAX_ELEMENT - (MAX_SIZE - size()) + 1;
  }
};

} // namespace zonotope

#endif // COMBINATION_BASE_HPP_

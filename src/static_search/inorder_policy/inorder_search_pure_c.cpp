/********************************************************************
*
* The pure-C implementation the binary search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __INORDER_SEARCH_PURE_C_CPP__
#define __INORDER_SEARCH_PURE_C_CPP__
#include <iterator>

namespace pure_c {
 template <typename RandomIterator, typename T>
 bool inorder_search(RandomIterator begin,
		     RandomIterator beyond,
		     const T& value) {
    typedef typename std::iterator_traits<RandomIterator>::difference_type diff_type;
  initialize:
    diff_type half, len = beyond - begin;
    T middle_value;
    RandomIterator middle;

    goto not_finished;
  descend_right:
    begin = middle;
    begin = begin + 1;
    len = len - half;
    len = len - 1;  
  not_finished:
    if (len == 0) goto return_false;
  node_contains:
    half = len >> 1;
    middle = begin + half;
    middle_value = *middle;
    if (middle_value == value) goto return_true;
    if (middle_value < value) goto descend_right;
  descend_left:
    len = half;
    goto not_finished;

  return_false:
    return false;
  return_true:
    return true;
 }
}
#endif //__INORDER_SEARCH_PURE_C_CPP__

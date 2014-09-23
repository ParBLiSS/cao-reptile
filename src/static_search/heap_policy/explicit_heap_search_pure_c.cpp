/********************************************************************
*
* The pure-C implementation the explicit heap search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __EXPLICIT_HEAP_SEARCH_PURE_C_CPP__
#define __EXPLICIT_HEAP_SEARCH_PURE_C_CPP__
#include <iterator>
#include <algorithm>
#include <vector>
#include <iostream>
#include "lower_bound_pure_c.h"

namespace pure_c {
 template <typename RandomIterator, typename T>
 bool explicit_heap_search(RandomIterator begin,
			   RandomIterator beyond,
			   unsigned int degree_minus_1,
			   const T& value) {
   typedef typename std::iterator_traits<RandomIterator>::difference_type diff_type;
  initialize:
   T x, l_middle_value, *lower_bound, *l_middle;
   diff_type y, l_len, l_half;
   RandomIterator current_node = begin;

   goto not_finished;
  descend_tree:
   y = lower_bound - reinterpret_cast<T*>(current_node);
   current_node = current_node->p[y];
  not_finished:
   if (current_node >= beyond) goto return_false;
  node_contains:
   lower_bound = reinterpret_cast<T*>(current_node);
   LOWER_BOUND(lower_bound, x)
   if (x == value) goto return_true;
   goto descend_tree;

  return_false:
   return false;
  return_true:
   return true;
 }
}
#endif //__EXPLICIT_HEAP_SEARCH_PURE_C_CPP__

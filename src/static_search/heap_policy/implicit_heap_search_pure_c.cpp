/********************************************************************
*
* The pure-C implementation the implicit heap search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __IMPLICIT_HEAP_SEARCH_PURE_C_CPP__
#define __IMPLICIT_HEAP_SEARCH_PURE_C_CPP__
#include <iterator>
#include <algorithm>
#include <vector>
#include "lower_bound_pure_c.h"

namespace pure_c {
 template <typename RandomIterator, typename T>
 bool implicit_heap_search(RandomIterator begin,
			   RandomIterator beyond,
			   int degree,
			   const T& value) {
   typedef typename std::iterator_traits<RandomIterator>::difference_type diff_type;
  initialize:
   RandomIterator lower_bound, l_middle, current_node = begin;
   diff_type l_len, l_half, index = 0, size = beyond-begin;
   T l_middle_value, x;
   int degree_minus_1 = degree - 1;

   goto not_finished;
  descend_tree:
   index = index * degree;
   x = lower_bound-current_node;
   x = x + 1;
   x = degree_minus_1 * x;
   index  = index + x;
   current_node = begin + index;
  not_finished:
   if (index >= size) goto return_false;
  node_contains:
   lower_bound = current_node;
   LOWER_BOUND(lower_bound, x)
   if (x == value) goto return_true;
   goto descend_tree;

  return_false:
   return false;
  return_true:
   return true;
 }
}

#endif //__IMPLICIT_HEAP_SEARCH_PURE_C_CPP__

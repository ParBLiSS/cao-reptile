/********************************************************************
*
* The pure-C implementation the explicit height partitioned
* search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __EXPLICIT_HP_SEARCH_PURE_C_CPP__
#define __EXPLICIT_HP_SEARCH_PURE_C_CPP__
#include <iterator>
#include <vector>
#include <iostream>

namespace pure_c {
 template <typename RandomIterator, typename T>
 bool explicit_hp_search(RandomIterator begin,
			 RandomIterator beyond,
			 const T& value) {
  initialize:
   typedef typename std::iterator_traits<RandomIterator>::difference_type diff_type;
   RandomIterator current_element = begin;
   T x;
   diff_type size = beyond - begin, current_bfs_index = 1;

   goto not_finished;
  descend_right:
   current_element = current_element->left_child;
   current_bfs_index = current_bfs_index*2;
   current_bfs_index++;
  not_finished:
   if (current_bfs_index > size) goto return_false;
  node_contains:
   x = current_element->e;
   if (x == value) goto return_true;
   if (x > value) goto descend_right;
  descend_left:
   current_element = current_element->right_child;
   current_bfs_index = current_bfs_index*2;
   goto not_finished;

  return_false:
   return false;
  return_true:
   return true;
 }
}

#endif //__EXPLICIT_HP_SEARCH_PURE_C_CPP__

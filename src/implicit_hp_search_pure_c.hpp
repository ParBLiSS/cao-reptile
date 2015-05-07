/********************************************************************
*
* The pure-C implementation the implicit height partitioned
* search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __IMPLICIT_HP_SEARCH_PURE_C_CPP__
#define __IMPLICIT_HP_SEARCH_PURE_C_CPP__
#include <iterator>
#include <vector>

namespace pure_c {
  const int max_depth = 40;
    template <typename RandomIterator, typename Tp>
 RandomIterator implicit_hp_search(RandomIterator begin,
			 RandomIterator beyond,
                                   //RandomIterator Pos,
			 RandomIterator T,
			 RandomIterator B,
			 RandomIterator D,
			 const Tp& value) {
   typedef typename std::iterator_traits<RandomIterator>::difference_type diff_type;
   //initialize:
   diff_type current_depth = 0, current_bfs_index = 1, size = beyond - begin;
   Tp x, y, z, w;
   Tp PosArr[max_depth] = {0};
   Tp* Pos = PosArr, *q;
   RandomIterator p;
   *Pos = 0;

   goto not_finished;
  descend_right:
   current_bfs_index = current_bfs_index*2;
   current_bfs_index++;
   current_depth++;
   if (current_bfs_index > size) goto return_false;
   p = D + current_depth;
   x = *p;
   q = Pos + x;
   x = *q;
   p = T + current_depth;
   y = *p;
   z = current_bfs_index & y;
   p = B + current_depth;
   w = *p;
   z = z * w;
   w = x + y;
   w = w + z;
   q = Pos + current_depth;
   *q = w;
  not_finished:
    if (current_bfs_index > size) goto return_false;
  //node_contains:
   q = Pos+current_depth;
   x = *q;
   p = begin+x;
   x = *p;
   if (x == value) goto return_true;
   if (x < value) goto descend_right;
  //descend_left:
   current_bfs_index = current_bfs_index*2;
   current_depth++;
   if (current_bfs_index > size) goto return_false;
   p = D + current_depth;
   x = *p;
   q = Pos + x;
   x = *q;
   p = T + current_depth;
   y = *p;
   z = current_bfs_index & y;
   p = B + current_depth;
   w = *p;
   z = z * w;
   w = x + y;
   w = w + z;
   q = Pos + current_depth;
   *q = w;
   goto not_finished;

  return_false:
   return beyond;
   //return false;
  return_true:
   return p;
   //return true;
 }
}
#endif //__IMPLICIT_HP_SEARCH_PURE_C_CPP__

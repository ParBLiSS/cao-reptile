/***************************************************************
 * The lower_bound implementations of this file are copied from
 * S. Mortensen, Refining the pure-C Cost Model, 
 * M. Sc. Thesis, Department of Computer Science, 
 * University of Copenhagen (2001). 
 ***************************************************************/

#ifndef __LOWER_BOUND_OPTIMIZED_PURE_C_CPP__
#define __LOWER_BOUND_OPTIMIZED_PURE_C_CPP__
#include <iterator>
#include "log2.h"

template<typename RandomIterator, typename T> 
RandomIterator lower_bound_optimized(RandomIterator begin,
				     RandomIterator end,
				     const T& val)
{
  typename std::iterator_traits<RandomIterator>::difference_type n = end - begin;
  if (n == 0) return end;
  ptrdiff_t i = (1 << floor_log2(n)) - 1;
  begin = begin[i] < val ? begin + (n - i) : begin;
  while (i > 0) {
    i = i >> 1;
    begin = begin[i] < val ? begin + i + 1: begin;
  }
  return begin;
}

#endif //__LOWER_BOUND_OPTIMIZED_PURE_C_CPP__

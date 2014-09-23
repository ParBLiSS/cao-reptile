/********************************************************************
*
* pure-C implementation of lower_bound
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __LOWER_BOUND_PURE_C__
#define __LOWER_BOUND_PURE_C__

#define LOWER_BOUND(lower_bound, x) \
   l_len = degree_minus_1; \
   goto l_not_finished; \
  l_right: \
   lower_bound = l_middle; \
   lower_bound = lower_bound + 1; \
   l_len = l_len - l_half; \
   l_len = l_len - 1; \
  l_not_finished: \
   if (l_len == 0) goto l_return; \
  l_found: \
   l_half = l_len >> 1; \
   l_middle = lower_bound + l_half; \
   l_middle_value = *l_middle; \
   if (l_middle_value < value) goto l_right; \
  l_left: \
   l_len = l_half; \
   goto l_not_finished; \
  l_return: \
   x = *lower_bound;

#endif //__LOWER_BOUND_PURE_C__

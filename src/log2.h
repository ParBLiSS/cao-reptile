/*******************************************************
 *
 * log2.h
 *
 * Implementation of ceil_log2 and hyperfloor. 
 * The hyperfloor(N) is defined as the largest power 
 * of 2 <= N.
 *
 * By Frederik Rønn, University of Copenhagen.
 *
 ******************************************************/

#ifndef __LOG2_H__
#define __LOG2_H__
#include <iostream>

inline size_t floor_log2(unsigned long n)
{
  unsigned char log_table[256] = {
  0xff,
  0,
  1,1,
  2,2,2,2,
  3,3,3,3,3,3,3,3,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
  long rv = 0;
  if (n & 0xffff0000) { rv += 16; n >>= 16; }
  if (n & 0xff00) { rv += 8; n >>= 8; }
  return rv + log_table[n];
}

inline size_t ceil_log2(unsigned long n)
{
  size_t i = floor_log2(n);
  return (n-(1<<i)==0)?i:i+1;
}

inline size_t hyperfloor(unsigned long n)
{
  return 1<<floor_log2(n);
}

template<typename SizeType>
inline SizeType int_log(const SizeType& b,
                        const SizeType& x){
    double st = b;
    SizeType vl = 1;
    while(st < x){
        st = st * b;
        vl++;
    }
    return vl;
}

#endif  //__LOG2_H__

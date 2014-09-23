/********************************************************************
*
* Defines a node of the explicit heap search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __ELEMENT_N_H__
#define __ELEMENT_N_H__
#include <vector>

template <typename Element, unsigned int N, size_t CacheLineSize>
struct element_N {
  Element e[N-1];
  element_N<Element, N, CacheLineSize> *p[N];
  char dummy[CacheLineSize-(N*sizeof(element_N<Element, N, CacheLineSize>*)+(N-1)*sizeof(Element))];
};

#endif //__ELEMENT_N_H__

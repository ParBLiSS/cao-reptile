/********************************************************************
*
* Generic Search program
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __GENERIC_SEARCH_CPP__
#define __GENERIC_SEARCH_CPP__

template <typename RandomIterator, typename T, typename LayoutPolicy>
bool generic_search(RandomIterator begin,
		    RandomIterator beyond,
		    LayoutPolicy policy,
		    const T& value) {
  policy.initialize(begin, beyond);
  while(policy.not_finished()) {
    if (policy.node_contains(value)) {
      return true;
    }
    else
      policy.descend_tree(value);
  }
  return false;
}

#endif //__GENERIC_SEARCH_CPP__

/********************************************************************
*
* Empty layout policy. 
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __EMPTY_POLICY_CPP__
#define __EMPTY_POLICY_CPP__

template<typename RandomIterator, typename Value>
class empty_policy {
  public:
    empty_policy() { }
    inline void initialize(RandomIterator begin,
                           RandomIterator beyond) { }
    inline bool not_finished() { }
    inline bool node_contains(const Value& value) { }
    inline void descend_tree(const Value& value) { }
};

#endif //__EMPTY_POLICY_CPP__

/********************************************************************
*
* The layout policy for the explicit heap search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __EXPLICIT_HEAP_LAYOUT_CPP__
#define __EXPLICIT_HEAP_LAYOUT_CPP__

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cmath>

using namespace std;

template <typename RandomIterator, typename Value>
class explicit_heap_policy {
  public:
    explicit_heap_policy(int degree) { m_degree = degree; }
    inline void initialize(RandomIterator begin,
			   RandomIterator beyond) { 
      m_current_node = begin; 
      m_beyond = beyond;
    }
    inline bool not_finished() {
      return (m_current_node < m_beyond);
    }
    inline bool node_contains(const Value &value){
      m_lower_bound = 
	lower_bound(&(m_current_node->e[0]), &(m_current_node->e[m_degree-1]), value);
      return (*m_lower_bound == value);
    }
    inline void descend_tree(const Value& element) {
      m_current_node = 
	m_current_node->p[m_lower_bound-reinterpret_cast<Value*>(m_current_node)];
    }
  private:
    int m_degree;
    RandomIterator m_current_node, m_beyond;
    Value *m_lower_bound;
};

#endif // __EXPLICIT_HEAP_LAYOUT_CPP__

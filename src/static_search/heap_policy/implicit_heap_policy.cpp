/********************************************************************
*
* The layout policy for the implicit heap search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __IMPLICIT_HEAP_POLICY_CPP__
#define __IMPLICIT_HEAP_POLICY_CPP__

#include <iostream>
#include <iterator>

using namespace std;

template <typename RandomIterator, typename Value>
class implicit_heap_policy {
  public:
    implicit_heap_policy(int degree) { m_degree = degree; }
    inline void initialize(RandomIterator begin,
                           RandomIterator beyond) {
      m_index = 0;
      m_size = beyond-begin-m_degree;
      m_current_node = begin;
      m_lower_bound = begin;
      m_begin = begin;
    }
    inline bool not_finished() { return m_index < m_size; }
    inline bool node_contains(const Value &value) {
      m_lower_bound = 
	std::lower_bound(m_current_node, m_current_node+m_degree-1, value);
      return (*m_lower_bound == value);
    }
    inline void descend_tree(const Value& value) {
      m_index = 
	m_index*m_degree + (m_degree-1)*((m_lower_bound-m_current_node)+1);
      m_current_node = m_begin+m_index;
    }
  private:
    int m_degree;
    typename iterator_traits<RandomIterator>::difference_type 
      m_index, m_size;
    RandomIterator m_current_node, m_lower_bound, m_begin;
};

#endif // __IMPLICIT_HEAP_POLICY_CPP__

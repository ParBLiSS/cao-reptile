/********************************************************************
*
* The layout policy for the explicit height partitioned
* search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __EXPLICIT_HEIGHT_PARTITIONING_POLICY_CPP__
#define __EXPLICIT_HEIGHT_PARTITIONING_POLICY_CPP__

#include "hp_precomputed_table.h"
#include "../misc/log2.h"
#include <iterator>

template <typename RandomIterator, typename Value>
class explicit_height_partitioning_policy {
  public:
    explicit_height_partitioning_policy() {}
    inline void initialize(RandomIterator begin,
                           RandomIterator beyond) { 
      m_current_element = begin; 
      m_current_bfs_index = 1;
      m_size = beyond-begin;
    }
    inline bool not_finished() { return (m_current_bfs_index < m_size); }
    inline bool node_contains(const Value& value) { return m_current_element->e == value; }
    inline void descend_tree(const Value& value) {
      m_current_bfs_index = m_current_bfs_index*2;
      if(m_current_element->e > value)
        m_current_element = m_current_element->left_child;
      else {
        m_current_element = m_current_element->right_child;
        m_current_bfs_index++;
      }
    }
  private:
    RandomIterator m_current_element;
    typename iterator_traits<RandomIterator>::difference_type 
      m_size, m_current_bfs_index;
};

#endif //__EXPLICIT_HEIGHT_PARTITIONING_POLICY_CPP__

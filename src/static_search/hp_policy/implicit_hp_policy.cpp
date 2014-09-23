/********************************************************************
*
* The layout policy for the implicit height partitioned
* search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __IMPLICIT_HEIGHT_PARTITIONING_POLICY_CPP__
#define __IMPLICIT_HEIGHT_PARTITIONING_POLICY_CPP__

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include "hp_precomputed_table.h"

template <typename RandomIterator, typename Value>
class implicit_height_partitioning_policy {
  public:
    implicit_height_partitioning_policy(precomputed_table<Value> *precomp_table) {
      m_precomp_table = precomp_table;
    }
    inline void initialize(RandomIterator begin,
                           RandomIterator beyond) {
      m_current_depth     = 0;
      m_begin             = begin;
      m_current_bfs_index = 1;
      m_size              = beyond-begin;
      m_precomp_table->initialise();
    }
    inline bool not_finished() { return m_current_bfs_index <= m_size; }
    inline bool node_contains(const Value& value) {
      return m_begin[m_precomp_table->Pos[m_current_depth]] == value;
    }
    inline void descend_tree(const Value& value) {
      if (m_begin[m_precomp_table->Pos[m_current_depth]] > value)
        m_current_bfs_index = m_current_bfs_index*2;
      else
        m_current_bfs_index = m_current_bfs_index*2+1;
      m_current_depth++;
      m_precomp_table->Pos[m_current_depth] = 
	m_precomp_table->Pos[m_precomp_table->D[m_current_depth]]+
        m_precomp_table->T[m_current_depth]+
        (m_current_bfs_index & m_precomp_table->T[m_current_depth])*
	m_precomp_table->B[m_current_depth];
    }
  private:
    precomputed_table<Value> *m_precomp_table;
    typename iterator_traits<RandomIterator>::difference_type 
      m_current_depth, m_current_bfs_index, m_size;
    RandomIterator m_begin;
};

#endif // __IMPLICIT_HEIGHT_PARTITIONING_POLICY_CPP__

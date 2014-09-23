/********************************************************************
*
* The layout policy for the binary search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#include <iterator>
#include <iostream>

#ifndef __INORDER_SEARCH_CPP__
#define __INORDER_SEARCH_CPP__

template<typename RandomIterator, typename Value>
class inorder_search_policy {
  public:
    inorder_search_policy() { }
    inline void initialize(RandomIterator begin,
                           RandomIterator beyond) {
      m_begin = begin;
      m_beyond = beyond;
      m_len = beyond-begin;
    }
    inline bool not_finished() { return (m_len > 0); }
    inline bool node_contains(const Value& value) {
      m_half = m_len >> 1;
      m_middle = m_begin;
      m_middle = m_middle + m_half;
      return *m_middle == value;
    }
    inline void descend_tree(const Value& value) {
      if (*m_middle < value) {
        m_begin = m_middle;
        ++m_begin;
        m_len = m_len - m_half - 1;
      }
      else
        m_len = m_half;
    }
  private:
    typename std::iterator_traits<RandomIterator>::difference_type m_half;
    typename std::iterator_traits<RandomIterator>::difference_type m_len;
    RandomIterator m_begin, m_beyond, m_middle;
};

#endif //__INORDER_SEARCH_CPP__

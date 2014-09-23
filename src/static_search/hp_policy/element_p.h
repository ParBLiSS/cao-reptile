/********************************************************************
*
* Defines a node of the explicit height partitioned search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __ELEMENT_P_H__
#define __ELEMENT_P_H__

template<typename Element>
struct element_p {
  Element e;
  element_p<Element> *left_child;
  element_p<Element> *right_child;
};

#endif // __ELEMENT_P_H__

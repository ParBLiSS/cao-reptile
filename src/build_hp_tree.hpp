/********************************************************************
*
* Builds a tree in the height-partitioned layout.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/
#ifndef _BUILD_HP_TREE_H_
#define _BUILD_HP_TREE_H_

#include "log2.h"

using namespace std;

template<typename RandomIterator>
void build_hp_tree(RandomIterator begin_in,
                     RandomIterator beyond_in,
                     RandomIterator begin_out,
                     RandomIterator beyond_out,
                     int height,
                     int inc) {
  int bottom_height = ((height==2)?1:hyperfloor(height-1));
  int top_height    = height-bottom_height;
  int bottom_size   = (1<<bottom_height)-1;
  int top_size      = (1<<top_height)-1;

  if (top_height==1 && bottom_height==1) {
    begin_out[1] = begin_in[0];
    begin_out[0] = begin_in[1*inc];
    begin_out[2] = begin_in[2*inc];
    return;
  }

  if (top_height==1) {
    begin_out[0] = begin_in[bottom_size*inc];
  } else {
    build_hp_tree(begin_in+bottom_size*inc,
                  beyond_in,
                  begin_out,
                  beyond_out,
                  top_height,
                  bottom_size*inc+inc);
  }

  for(int i=0;i<=top_size;i++) {
    build_hp_tree(begin_in+(i*bottom_size+i)*inc,
         	  beyond_in+(i+1)*bottom_size*inc+i,
         	  begin_out+top_size+i*bottom_size,
         	  beyond_out,
         	  bottom_height,
          	  inc);
  }
}

template<typename RandomInIterator,
         typename RandomIdIterator,
         typename RandomCountIterator>
void coblivious_layout(RandomInIterator begin_in,
                       RandomInIterator beyond_in,
                       RandomIdIterator begin_id,
                       RandomIdIterator beyond_id,
                       RandomCountIterator begin_count,
                       RandomCountIterator beyond_count,
                       int height, int inc) {
  int bottom_height = ((height==2)?1:hyperfloor(height-1));
  int top_height    = height-bottom_height;
  int bottom_size   = (1<<bottom_height)-1;
  int top_size      = (1<<top_height)-1;

  if (top_height==1 && bottom_height==1) {
      begin_id[1] = begin_in[0].ID;
      begin_id[0] = begin_in[1*inc].ID;
      begin_id[2] = begin_in[2*inc].ID;
      begin_count[1] = begin_in[0].count;
      begin_count[0] = begin_in[1*inc].count;
      begin_count[2] = begin_in[2*inc].count;
    return;
  }

  if (top_height==1) {
    begin_id[0] = begin_in[bottom_size*inc].ID;
    begin_count[0] = begin_in[bottom_size*inc].count;
  } else {
    coblivious_layout(begin_in+bottom_size*inc,
                  beyond_in,
                  begin_id,
                  beyond_id,
                  begin_count,
                  beyond_count,
                  top_height,
                  bottom_size*inc+inc);
  }

  for(int i=0;i<=top_size;i++) {
    coblivious_layout(begin_in+(i*bottom_size+i)*inc,
         	  beyond_in+(i+1)*bottom_size*inc+i,
         	  begin_id+top_size+i*bottom_size,
         	  beyond_id,
         	  begin_count+top_size+i*bottom_size,
         	  beyond_count,
         	  bottom_height,
          	  inc);
  }
}

#endif /* _BUILD_HP_TREE_H_ */

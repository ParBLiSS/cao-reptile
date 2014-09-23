/********************************************************************
*
* Builds a tree in the height-partitioned layout.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

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

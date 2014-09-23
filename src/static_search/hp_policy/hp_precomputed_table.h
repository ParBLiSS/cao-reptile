/********************************************************************
*
* Defines the precomputed table that holds information about the
* layout of a height partitioned search tree.
*
* Frederik Rønn, June 2003, University of Copenhagen.
*
*********************************************************************/

#ifndef __HEIGHT_PARTITIONING_PRECOMPUTED_TABLE__
#define __HEIGHT_PARTITIONING_PRECOMPUTED_TABLE__

#include <vector>
#include "../misc/log2.h"

using namespace std;

template<typename element>
class precomputed_table {
  public:
   precomputed_table(element height) : D(height),
                                       T(height),
                                       B(height),
                                       Pos(height) {
     D.resize(height);
     T.resize(height);
     B.resize(height);
     Pos.resize(height);
     build_precomputed_table(height, 0);
   }

    void initialise() {
      Pos[0] = 0;
    };

    vector<element> D;
    vector<element> T;
    vector<element> B;
    vector<element> Pos;

  private:
    void build_precomputed_table(element height, element depth) {
      element bottom_height = ((height==2)?1:hyperfloor(height-1));
      element top_height    = height-bottom_height;
      
      if (height==1) return;
      
      D[depth+top_height]   = depth;
      T[depth+top_height]   = (1<<top_height)-1;
      B[depth+top_height]   = (1<<bottom_height)-1;
      Pos[depth+top_height] = 0;
      
      build_precomputed_table(top_height, depth);
      build_precomputed_table(bottom_height, depth+top_height);
    }
};

#endif //__HEIGHT_PARTITIONING_PRECOMPUTED_TABLE__

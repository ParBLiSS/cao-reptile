#ifndef _FLAT_LAYOUT_H_
#define _FLAT_LAYOUT_H_
#include <vector>
#include <iterator>
#include <algorithm>
#include "cao_util.h"

template <typename RandomIterator,
          typename IdType,
          typename CountType>
class ECDataFlatLayout{
public:
    typedef typename std::vector<IdType>::iterator id_iterator;
    std::vector<IdType> mIds;
    std::vector<CountType> mCounts;
    ECDataFlatLayout(){};

    void init(RandomIterator intr,
              size_type nTotal, int kSize){
        assert(nTotal > 0);
        mIds.resize(nTotal);
        mCounts.resize(nTotal);
        for(size_type i = 0u; i < nTotal;i++) {
	  mIds[i] = (intr + i)->ID;
          mCounts[i] = (intr + i)->count;
        }
    }

    bool find(const IdType &rID) {
      int final = -1;
      int lb = 0, ub = mIds.size() - 1, mid;
      while (lb <= ub) {
        mid = (lb + ub) / 2;
        if (mIds[mid] == rID) {
            final = mid;
            break;
        }
        else if (mIds[mid] < rID)
            lb = mid + 1;
        else if (mIds[mid] > rID)
            ub = mid - 1;
      }
      return (final != -1);
    };
 
   int getCount(const IdType& rID, CountType& count){
      int final = -1;
      int lb = 0, ub = mIds.size() - 1, mid;
      while (lb <= ub) {
        mid = (lb + ub) / 2;
        if (mIds[mid] == rID) {
            final = mid;
            break;
        }
        else if (mIds[mid] < rID)
            lb = mid + 1;
        else if (mIds[mid] > rID)
            ub = mid - 1;
      }
      if(final != -1){
            count = mCounts[final];
            //std::cout << "count " << (int) count << std::endl;
            return final;
      } else {
            return -1;
      }
    }
};

#endif // 

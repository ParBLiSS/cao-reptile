#ifndef _FLAT_LAYOUT_H_
#define _FLAT_LAYOUT_H_
#include <vector>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <cassert>
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
              size_type nTotal){
        assert(nTotal > 0);
        mIds.resize(nTotal);
        mCounts.resize(nTotal);
        for(size_type i = 0u; i < nTotal;i++) {
            mIds[i] = (intr + i)->ID;
            mCounts[i] = (intr + i)->count;
        }
    }

    bool find(const IdType &rID) const{
      long final = -1;
      long lb = 0, ub = mIds.size() - 1, mid;
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

    long getCount(const IdType& rID, CountType& count) const{
      long final = -1;
      long lb = 0, ub = mIds.size() - 1, mid;
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

    void serialize(std::ofstream& ofs) const{
        ofs.write((char*)&mIds[0], sizeof(IdType)*mIds.size());
        ofs.write((char*)&mCounts[0], sizeof(CountType)*mCounts.size());
    }

};

#endif //

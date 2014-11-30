#ifndef _COBLIVIOUS_LAYOUT_H_
#define _COBLIVIOUS_LAYOUT_H_
#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <stdint.h>
#include "build_hp_tree.hpp"
#include "hp_precomputed_table.h"
#include "implicit_hp_search_pure_c.hpp"
#include "io_util.h"
#include "log2.h"
#include "cao_util.h"

template <typename RandomInputItr,
          typename IdType,
          typename CountType>
class ECDataCOLayout{

    void log2_steps(size_type n, std::vector<IdType>& steps,
                    size_type limit){
        size_type mn = floor_log2(n);
        std::cout << limit << "\t" << mn << "\t";
        while(mn > 2){
            if((1 << mn) & n){
                steps.push_back(mn);
                limit -= 1;
            }
            mn -= 1;
            std::cout << mn << "\t";
            if(limit == 0)
                break;
        }
        std::cout << std::endl;
    }

    void init_layout(RandomInputItr inItrBegin, size_type nTotal){
        RandomInputItr inItr = inItrBegin;
        id_iterator idItr = mIds.begin();
        count_iterator countItr = mCounts.begin();
        for(size_t i = 0u; i < treeSteps.size();i++){
            size_t cstep = 1 << treeSteps[i];
            coblivious_layout(inItr, inItr + cstep,
                              idItr, idItr + cstep,
                              countItr, countItr + cstep,
                              treeSteps[i], 1);
            precomputed_table<IdType> pct(treeSteps[i]);

            stTotal = stTotal + cstep - 1;
            inItr = inItr + cstep - 1;
            idItr = idItr + cstep - 1;
            countItr = countItr + cstep - 1;
            pctSteps.push_back(pct);
        }
        // layout out side of steps
        for(size_t i = stTotal; i < nTotal;i++){
            mIds[i] = (inItrBegin + i)->ID;
            mCounts[i] = (inItrBegin + i)->count;
        }
    }
    void init_steps(RandomInputItr inItrBegin){
        // asssumes that treesteps are initialized
        assert(treeSteps.size() > 0);
        stepIds.resize(treeSteps.size());
        stepPrefixes.resize(treeSteps.size());
        stepPrefixes[0] = (1 << treeSteps[0]) - 1;
        stepIds[0] = (inItrBegin + stepPrefixes[0])->ID;

        // init step mIds
        for(size_t i = 1u; i < treeSteps.size();i++) {
            stepPrefixes[i] = stepPrefixes[i - 1] + (1 << treeSteps[i]) - 1;
            stepIds[i] = (inItrBegin + stepPrefixes[i])->ID;
        }
    }

    void init_steps(RandomInputItr inItrBegin, size_type nTotal,
                    size_type limit){
        stTotal = 0;
        // init steps
        log2_steps(nTotal, treeSteps, limit);
        //
        if(treeSteps.size() > 0) {
            init_steps(inItrBegin);
        }
    }

public:

    typedef typename std::vector<IdType>::iterator id_iterator;
    typedef typename std::vector<CountType>::iterator count_iterator;
    typedef size_type_vector::iterator sizet_iterator;
    std::vector< precomputed_table<IdType> > pctSteps;
    std::vector<IdType> treeSteps;
    std::vector<IdType> stepPrefixes;
    std::vector<IdType> mIds;
    std::vector<IdType> stepIds;
    std::vector<CountType> mCounts;
    size_type stTotal;

    ECDataCOLayout(){ }

    void init(RandomInputItr inItrBegin,
              size_type nTotal, size_type limit){
        assert(limit > 0);
        assert(nTotal > limit);
        mIds.resize(nTotal);
        mCounts.resize(nTotal);
        // init steps
        init_steps(inItrBegin, nTotal, limit);
        // layout of each step
        init_layout(inItrBegin, nTotal);
    }

    ECDataCOLayout(RandomInputItr inItrBegin,
                   size_type nTotal,
                   size_type limit){
        init(inItrBegin, nTotal, limit);
    }

    id_iterator getPosition(const IdType& val){
        id_iterator idtr = std::upper_bound(stepIds.begin(),
                                            stepIds.end(), val);
        if(idtr == stepIds.end()){
            //std::cout << "END" << std::endl;
            // may do linear search, because it is very small
            idtr = std::lower_bound(mIds.begin() + stTotal,
                                    mIds.end(), val);
        } else {
            size_type dist = idtr - stepIds.begin(),
                bdist = 0, edist = stepPrefixes[dist];
            //std::cout << "DIST " <<  dist << " " << edist << std::endl;
            if(dist > 0){
                bdist = stepPrefixes[dist - 1];
            }
            idtr = pure_c::implicit_hp_search(mIds.begin() + bdist,
                                              mIds.begin() + edist,
                                              pctSteps[dist].Pos.begin(),
                                              pctSteps[dist].T.begin(),
                                              pctSteps[dist].B.begin(),
                                              pctSteps[dist].D.begin(), val);
        }
        if(idtr != mIds.end() && ((*idtr) == val))
            return idtr;
        else
            return mIds.end();
    }

    int getCount(const IdType& val, CountType& count){
        id_iterator idtr = getPosition(val);
        if(idtr != mIds.end()){
            int dist = idtr - mIds.begin();
            count = mCounts[dist];
            //std::cout << std::endl << (int)count << std::endl;
            return dist;
        } else {
            count = 0;
            return -1;
        }
    }

    bool find(const IdType& val){
        return (getPosition(val) != mIds.end());
    }

};
int test_oblivious();

#endif /* _COBLIVIOUS_LAYOUT_H_ */

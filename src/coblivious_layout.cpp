#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <stdint.h>
#include "coblivious_layout.hpp"
#include "build_hp_tree.hpp"
#include "hp_precomputed_table.h"
#include "implicit_hp_search_pure_c.hpp"
#include "io_util.h"
#include "log2.h"
#include "cao_util.h"


template <typename RandomInputItr,
          typename IdType,
          typename CountType>
class KmerCOLayout{

    void log2_steps(size_type n, std::vector<IdType>& steps,
                    size_type limit = 10){
        size_type mn = floor_log2(n);
        while(mn > limit){
            if((1 << mn) & n)
                steps.push_back(mn);
            mn = mn - 1;
        }
    }

    void init_layout(RandomInputItr inItrBegin, size_type nTotal,
                     size_type limit){
        RandomInputItr inItr = inItrBegin;
        id_iterator idItr = kmerIds.begin();
        count_iterator countItr = kmerCounts.begin();
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
            kmerIds[i] = (inItrBegin + i)->ID;
            kmerCounts[i] = (inItrBegin + i)->count;
        }
    }

    void init_steps(RandomInputItr inItrBegin, size_type nTotal,
                    size_type limit){
        stTotal = 0;
        // init steps
        log2_steps(nTotal, treeSteps, limit);
        if(treeSteps.size() > 0) {
            stepKmerIds.resize(treeSteps.size());
            stepPrefixes.resize(treeSteps.size());
            stepPrefixes[0] = (1 << treeSteps[0]) - 1;
            stepKmerIds[0] = (inItrBegin + stepPrefixes[0])->ID;
        }
        // init step kmerIds
        for(size_t i = 1u; i < treeSteps.size();i++) {
            stepPrefixes[i] = stepPrefixes[i - 1] + (1 << treeSteps[i]) - 1;
            stepKmerIds[i] = (inItrBegin + stepPrefixes[i])->ID;
        }
    }

public:

    typedef typename std::vector<IdType>::iterator id_iterator;
    typedef typename std::vector<CountType>::iterator count_iterator;
    typedef size_type_vector::iterator sizet_iterator;
    std::vector< precomputed_table<IdType> > pctSteps;
    std::vector<IdType> treeSteps;
    std::vector<IdType> stepPrefixes;
    std::vector<IdType> kmerIds;
    std::vector<IdType> stepKmerIds;
    std::vector<CountType> kmerCounts;
    size_type stTotal;


    KmerCOLayout(RandomInputItr inItrBegin,
                 size_type nTotal,
                 size_type limit){
        assert(limit > 0);
        assert(nTotal > limit);
        kmerIds.resize(nTotal);
        kmerCounts.resize(nTotal);
        // init steps
        init_steps(inItrBegin, nTotal, limit);
        // layout of each step
        init_layout(inItrBegin, nTotal, limit);
    }

    bool find(IdType val){
        // int dist = 0;
        // bool r = pure_c::implicit_hp_search(kmerIds.begin(),
        //                                   kmerIds.begin() + 31,
        //                                   pct.Pos.begin(),
        //                                   pct.T.begin(),
        //                                   pct.B.begin(),
        //                                   pct.D.begin(), val);
        // return r;
        id_iterator idtr = std::upper_bound(stepKmerIds.begin(),
                                            stepKmerIds.end(), val);
        if(idtr == stepKmerIds.end()){
            //std::cout << "END" << std::endl;
            // may do linear search, because it is very small
            return std::binary_search(kmerIds.begin() + stTotal,
                                      kmerIds.end(), val);
        } else {
            size_type dist = idtr - stepKmerIds.begin(),
                bdist = 0, edist = stepPrefixes[dist];
            //std::cout << "DIST " <<  dist << " " << edist << std::endl;
            if(dist > 0){
                bdist = stepPrefixes[dist - 1];
            }
            return pure_c::implicit_hp_search(kmerIds.begin() + bdist,
                                              kmerIds.begin() + edist,
                                              pctSteps[dist].Pos.begin(),
                                              pctSteps[dist].T.begin(),
                                              pctSteps[dist].B.begin(),
                                              pctSteps[dist].D.begin(), val);
        }
    }

};

typedef struct kd_s{
    unsigned int ID;
    unsigned char count;
} kd_t;


int test_oblivious(){
    size_t m = 1, beg = 32,
        end = 65;
    for(size_t i = beg; i  < end; i += m)
    {
        size_t n = i;

        std::vector<unsigned int> idt;
        std::vector<char> ct;
        std::vector<kd_t> tst;

        tst.resize(n);
        idt.resize(n);
        ct.resize(n);

        for(size_t i = 0u; i < n; i++) {
            tst[i].ID = i+1;
            tst[i].count = n - i;
        }


        KmerCOLayout< std::vector<kd_t>::iterator, unsigned int, char >
            colayout(tst.begin(), n, 1);
        write_stvector(colayout.treeSteps, std::cout);
        write_stvector(colayout.stepPrefixes, std::cout);
        write_stvector(colayout.stepKmerIds, std::cout);
        write_stvector(colayout.kmerIds, std::cout);

        std::cout << colayout.find(n + 2) << " "
                  << colayout.find(n + 1) << " "
                  << colayout.find(0) << " "
                  << colayout.find(1) << " "
                  << colayout.find(31) << " "
                  << colayout.find(32) << " "
                  << colayout.find(33) << " "
                  << colayout.find(46) << " "
                  << colayout.find(47) << " "
                  << colayout.find(48) << " "
                  << std::endl;
        // std::cout << colayout.find(32)
        //           << std::endl;
    }
    return 0;

}

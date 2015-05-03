#include "util.h"
#include "WorkDistribution.hpp"
#include "ECImpl.hpp"

#include <fstream>
#include <string>
#include <mpi.h>


// Make ReadStore swappable
static void swap(ReadStore& left, ReadStore& right){
    left.swap(right);
}


// Error correction for a batch of reads
struct BatchEC{
    void operator()(int tid, int rank, const ReadStore& rbatch, ECImpl& ecr)
    {
        int rID = rbatch.readId;
        for(unsigned i = 0; i < rbatch.size(); i++){
            int position = rbatch.readsOffset[i];
            int qposition = rbatch.qualsOffset[i];
            char* addr = const_cast<char*> (&(rbatch.readsString[position]));
            char* qAddr = const_cast<char*> (&(rbatch.qualsString[qposition]));
            ecr.readEC(rID, addr, qAddr);
            rID++;
        }
    // TODO: think about writing output when ecr has a lot of error records
    }
};

// Error correction for a single read within the batch
struct UnitEC{
    void operator()(unsigned ipos, const ReadStore& rbatch, ECImpl& ecr){
        if(ipos > rbatch.size())
            return;
        int position = rbatch.readsOffset[ipos];
        int qposition = rbatch.qualsOffset[ipos];
        char* addr = const_cast<char*> (&(rbatch.readsString[position]));
        char* qAddr = const_cast<char*> (&(rbatch.qualsString[qposition]));
        ecr.readEC(rbatch.readId + ipos, addr, qAddr);
        // TODO: think about writing output when ecr has a lot of error records
    }
};

// Load a batch of reads starting from offset
struct ReadBatchLoader{
    void operator()(const ECImpl& ecr, unsigned long woffStart,
                     unsigned long woffEnd, ReadStore& tmp){
        std::ifstream fin(ecr.getPara().iFaName);
        fin.seekg(woffStart);
        readBatch(&fin, UINT_MAX, woffEnd, tmp);
    }
};

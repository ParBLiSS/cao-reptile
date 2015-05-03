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
    void operator()(const ReadStore& rbatch, ECImpl& ecr, int, int) {
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
    void operator()(const ReadStore& rbatch, ECImpl& ecr, unsigned ipos) {
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

int getTotalWork(const std::string& fileName){
    std::ifstream fin(fileName.c_str());
    fin.seekg(0,std::ios::end);
    return fin.tellg();
}

void ec_wdist(ECData& ecdata){
    Para& params = ecdata.getParams();
    unsigned nThreads = (unsigned) params.numThreads;
    std::vector<ECImpl> threadEC(nThreads + 1,
                                 ECImpl(ecdata, params));
    unsigned tWork = getTotalWork(params.iFaName);
    unsigned nWorkers = params.m_size * params.numThreads;
    unsigned tChunk = tWork / (nWorkers * params.workFactor);
    WorkDistribution<ReadStore, unsigned long, ECImpl,
                     ReadBatchLoader, BatchEC, UnitEC> wdist(tWork, threadEC,
                                                             params.numThreads,
                                                             tChunk);
    if(params.m_rank == 0){
        wdist.masterMain();
    } else{
        wdist.slaveMain();
    }
}

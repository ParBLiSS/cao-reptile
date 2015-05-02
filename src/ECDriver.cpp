/***
 **
 *  File: ECDriver.cpp
 *  Created: Dec 12, 2009 4:05 PM
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *  Copyright 2011-12 Ankit Shah, Sriram P C, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1)
 *
 *  Reptile is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Reptile is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <set>
#include <cmath>
#include <algorithm>
#include <thread>

#include "util.h"
#include "ECImpl.hpp"
#include "ECDriver.hpp"

void print(int i) {
    std::cout << " " << i;
}

void printHex(int i) {
    std::cout << " " << std::hex << i;
}

void ECDriver::ec() {
    if(inPara_.writeOutput != 0)
        oHandle_.open(outFname_.c_str());

    if(inPara_.storeReads) {
        // reads are already obtained from input file and stored in ECData
        processBatch(ecdata_.getReads(),ecdata_.getQuals(),
                     ecdata_.getReadsOffsets(),ecdata_.getQualsOffsets());
    } else {
        processReadsFromFile();
    }
    if(oHandle_.good())
        oHandle_.close();
}

void ECDriver::processBatch(const cvec_t &ReadsString,const cvec_t &QualsString,
                            const ivec_t &ReadsOffset,const ivec_t &QualsOffset) {
  if(inPara_.numThreads > 1){
    processBatchMT(ReadsString, QualsString, ReadsOffset, QualsOffset);
  } else {
    processBatchST(ReadsString, QualsString, ReadsOffset, QualsOffset);
  }
}

void ECDriver::processBatchST(const cvec_t &ReadsString,const cvec_t &QualsString,
                              const ivec_t &ReadsOffset,const ivec_t &QualsOffset) {
    ECImpl ecr(ecdata_, inPara_);
#ifdef DEBUG
    std::stringstream out3 ;
#endif
    for(unsigned long i = 0; i < ReadsOffset.size();i++) {
        int position = ReadsOffset[i];
        int qposition = QualsOffset[i];
        char* addr = const_cast<char*> (&ReadsString[position]);
        char* qAddr = const_cast<char*> (&QualsString[qposition]);
#ifdef DEBUG
        out3 << addr << std::endl;
#endif
        ecr.readEC(readID_, addr, qAddr);
        readID_++;
    }

#ifdef DEBUG
    for(int i = 0; i< size;i++){
        if(rank == i){
            std::cout << out3.str();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    if(oHandle_.good())
        ecr.writeErrors(oHandle_);
}

void ec_thread(int tid, int nthreads, int grid, ECImpl& ecr,
               const cvec_t &ReadsString,const cvec_t &QualsString,
               const ivec_t &ReadsOffset,const ivec_t &QualsOffset)
{
  //std::cout << "Launched by thread" << tid << std::endl;
  // block decomposition of work.
  unsigned long ntotal =  ReadsOffset.size();
  unsigned long nbegin = BLOCK_LOW(tid, nthreads, ntotal);
  unsigned long nend = BLOCK_HIGH(tid, nthreads, ntotal);
  int rID = grid + nbegin; // global read id + the id we start with
  for(unsigned long i = nbegin; i <= nend;i++) {
    int position = ReadsOffset[i];
    int qposition = QualsOffset[i];
    char* addr = const_cast<char*> (&ReadsString[position]);
    char* qAddr = const_cast<char*> (&QualsString[qposition]);
    ecr.readEC(rID, addr, qAddr);
    rID++;
  }

}

void ECDriver::processBatchMT(const cvec_t &ReadsString,const cvec_t &QualsString,
                              const ivec_t &ReadsOffset,const ivec_t &QualsOffset) {
    std::vector<std::thread> tvx(inPara_.numThreads);
    std::vector<ECImpl> threadEC(inPara_.numThreads,
                                 ECImpl(ecdata_, inPara_));

    //Launch a group of threads
    for (int i = 0; i < inPara_.numThreads; ++i) {
        tvx[i] = std::thread(ec_thread, i, readID_, inPara_.numThreads,
                             std::ref(threadEC[i]),
                             std::cref(ReadsString), std::cref(QualsString),
                             std::cref(ReadsOffset), std::cref(QualsOffset));
    }

    std::cout << "Launched from the main" << std::endl;

    //Join the threads with the main thread
    for (int i = 0; i < inPara_.numThreads; ++i)
        tvx[i].join();

    // output is written serially
    for (int i = 0; i < inPara_.numThreads; ++i)
        if(oHandle_.good())
            threadEC[i].writeErrors(oHandle_);
}

void ECDriver::processReadsFromFile() {
    std::ifstream read_stream(inPara_.iFaName.c_str());
    assert(read_stream.good() == true);

    read_stream.seekg(inPara_.offsetStart, std::ios::beg);

    cvec_t ReadsString;   // full string
    cvec_t QualsString;   // full quality score
    ivec_t ReadsOffset;
    ivec_t QualsOffset;

    while(1){
        bool lastRead = readBatch(&read_stream, inPara_.batchSize,
                                  inPara_.offsetEnd, ReadsString, ReadsOffset,
                                  QualsString, QualsOffset, readID_);
#ifdef DEBUG
        {
            std::stringstream out ;
            out << "PROC : " << inPara_.mpi_env->rank()  << " "
                << ReadsOffset.size() << std::endl;
            std::cout << out.str();
        }
#endif
        assert(ReadsOffset.size() == QualsOffset.size());

        processBatch(ReadsString,QualsString,ReadsOffset,QualsOffset);

       // std::cout << out3.str();
        ReadsString.resize(0);
        QualsString.resize(0);
        ReadsOffset.resize(0);
        QualsOffset.resize(0);
        if(lastRead) break;
    }

    return;
}

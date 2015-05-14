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
#include "ECWDist.hpp"

void print(int i) {
    std::cout << " " << i;
}

void printHex(int i) {
    std::cout << " " << std::hex << i;
}

void ECDriver::ecDynamic(){
    ec_wdist(ecdata_);
}

void ECDriver::ecStatic() {
    if(inPara_.writeOutput != 0)
        oHandle_.open(inPara_.outputFilename.c_str());

    if(inPara_.storeReads) {
        // reads are already obtained from input file and stored in ECData
        correctBatch(ecdata_.getReadStore());
    } else {
        correctReadsFromFile();
    }

    if(oHandle_.good())
        oHandle_.close();
}

void ECDriver::ec(){
    if(inPara_.dynamicWorkDist == 0){
        ecStatic();
    } else {
        ecDynamic();
    }
}

void ECDriver::correctBatch(const ReadStore& rbatch) {
  if(inPara_.numThreads > 1){
      correctBatchMT(rbatch);
  } else {
      correctBatchST(rbatch);
  }
}

void ECDriver::correctBatchST(const ReadStore& rbatch) {
    readID_ = rbatch.readId;
    ECImpl ecr(ecdata_, inPara_);
#ifdef DEBUG
    std::stringstream out3 ;
#endif
    for(unsigned long i = 0; i < rbatch.readsOffset.size();i++) {
        int position = rbatch.readsOffset[i];
        int qposition = rbatch.qualsOffset[i];
        char* addr = const_cast<char*> (&(rbatch.readsString[position]));
        char* qAddr = const_cast<char*> (&(rbatch.qualsString[qposition]));
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

void ec_thread(int tid, int nthreads, const ReadStore& rbatch, ECImpl& ecr)
{
  //std::cout << "Launched by thread" << tid << std::endl;
  // block decomposition of work.
  unsigned long ntotal =  rbatch.readsOffset.size();
  unsigned long nbegin = BLOCK_LOW(tid, nthreads, ntotal);
  unsigned long nend = BLOCK_HIGH(tid, nthreads, ntotal);
  int rID = rbatch.readId + nbegin; // global read id + the id we start with
  for(unsigned long i = nbegin; i <= nend;i++) {
    int position = rbatch.readsOffset[i];
    int qposition = rbatch.qualsOffset[i];
    char* addr = const_cast<char*> (&(rbatch.readsString[position]));
    char* qAddr = const_cast<char*> (&(rbatch.qualsString[qposition]));
    ecr.readEC(rID, addr, qAddr);
    rID++;
  }
}

void ECDriver::startThreads(std::vector<std::thread>& tvx,
                            std::vector<ECImpl>& threadEC,
                            const ReadStore& rbatch){
    //Launch a group of threads
    for (int i = 0; i < inPara_.numThreads; ++i) {
        tvx[i] = std::thread(ec_thread, i, inPara_.numThreads,
                             std::cref(rbatch),
                             std::ref(threadEC[i]));
    }
}

void ECDriver::joinThreads(std::vector<std::thread>& tvx){
    for (int i = 0; i < inPara_.numThreads; ++i)
        tvx[i].join();
}

void ECDriver::writeErrors(std::vector<ECImpl>& threadEC){
    for (int i = 0; i < inPara_.numThreads; ++i)
        if(oHandle_.good())
            threadEC[i].writeErrors(oHandle_);
}

void ECDriver::correctBatchMT(const ReadStore& rbatch) {
    std::vector<std::thread> tvx(inPara_.numThreads);
    std::vector<ECImpl> threadEC(inPara_.numThreads,
                                 ECImpl(ecdata_, inPara_));

    startThreads(tvx, threadEC, rbatch);
    //std::cout << "Launched from the main" << std::endl;

    //Join the threads with the main thread
    joinThreads(tvx);
    // output is written serially
    writeErrors(threadEC);
}

void ECDriver::correctReadsFromFile() {
  if(inPara_.numThreads > 1){
      correctReadsFromFileMT();
  } else {
      correctReadsFromFileST();
  }
}

void ECDriver::correctReadsFromFileMT(){
    std::ifstream read_stream(inPara_.iFaName.c_str());
    std::vector<std::thread> tvx(inPara_.numThreads);
    std::vector<ECImpl> threadEC(inPara_.numThreads,
                                 ECImpl(ecdata_, inPara_));

    ReadStore rbatch[2];
    int cbatch = 0;
    bool lastRead = false;
    assert(read_stream.good() == true);
    read_stream.seekg(inPara_.offsetStart, std::ios::beg);
    // load first batch
    lastRead = readBatch(&read_stream, inPara_.batchSize,
                         inPara_.offsetEnd, rbatch[0]);
    while(lastRead == false){
        // start threads
        startThreads(tvx, threadEC, rbatch[cbatch]);
        // load the next set of reads ?
        int nbatch = (cbatch + 1) % 2;

        rbatch[nbatch].reset();
        lastRead = readBatch(&read_stream, inPara_.batchSize,
                             inPara_.offsetEnd, rbatch[nbatch]);
        // join threads
        joinThreads(tvx);
        cbatch = nbatch;
    }
    // last batch
    startThreads(tvx, threadEC, rbatch[cbatch]);
    joinThreads(tvx);
    // output is written serially
    writeErrors(threadEC);
}

void ECDriver::correctReadsFromFileST() {
    std::ifstream read_stream(inPara_.iFaName.c_str());
    assert(read_stream.good() == true);

    read_stream.seekg(inPara_.offsetStart, std::ios::beg);

    ReadStore rbatch;
    while(1) {
        bool lastRead = readBatch(&read_stream, inPara_.batchSize,
                                  inPara_.offsetEnd, rbatch);
#ifdef DEBUG
        {
            std::stringstream out ;
            out << "PROC : " << inPara_.m_rank  << " "
                << rbatch.readsOffset.size() << std::endl;
            std::cout << out.str();
        }
#endif
        assert(rbatch.readsOffset.size() == rbatch.qualsOffset.size());

        correctBatch(rbatch);

        rbatch.reset();

        if(lastRead) break;
    }

    return;
}

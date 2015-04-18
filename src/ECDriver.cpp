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
    readID_  = inPara_.startFromLineNo;
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
    processBatchST(ReadsString, QualsString, ReadsOffset, QualsOffset);
  } else {
    processBatchMT(ReadsString, QualsString, ReadsOffset, QualsOffset);
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

void ec_thread(int tid, int nthreads, ECImpl& ecr,
               const cvec_t &ReadsString,const cvec_t &QualsString,
               const ivec_t &ReadsOffset,const ivec_t &QualsOffset)
{
  std::cout << "Launched by thread" << tid << std::endl;
  // block decomposition of work.
  unsigned long ntotal =  ReadsOffset.size();
  unsigned long nbegin = BLOCK_LOW(tid, nthreads, ntotal);
  unsigned long nend = BLOCK_HIGH(tid, nthreads, ntotal);
  int rID = nbegin; // rID is NOT the same as nbegin, but for our expts, it is OK
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
  std::thread t[inPara_.numThreads];
  std::vector<ECImpl> threadEC(inPara_.numThreads, 
                               ECImpl(ecdata_, inPara_));

  //Launch a group of threads
  for (int i = 0; i < inPara_.numThreads; ++i) {
    t[i] = std::thread(ec_thread, i, inPara_.numThreads, std::ref(threadEC[i]),
                       std::cref(ReadsString), std::cref(QualsString),
                       std::cref(ReadsOffset), std::cref(QualsOffset));
  }

  std::cout << "Launched from the main" << std::endl;

  //Join the threads with the main thread
  for (int i = 0; i < inPara_.numThreads; ++i) 
      t[i].join();

  // output is written serially
  for (int i = 0; i < inPara_.numThreads; ++i) 
      if(oHandle_.good())
        threadEC[i].writeErrors(oHandle_);

}

void ECDriver::processReadsFromFile() {
    std::ifstream read_stream(inPara_.iFaName.c_str());
    assert(read_stream.good() == true);

    std::ifstream qual_stream(inPara_.iQName.c_str());
    assert(qual_stream.good() == true);

    read_stream.seekg(inPara_.offsetStart,std::ios::beg);
    qual_stream.seekg(inPara_.qOffsetStart,std::ios::beg);
    bIO::FASTA_input fasta(read_stream);
    bIO::FASTA_input qual(qual_stream);
    ++fasta;++qual;

    cvec_t ReadsString;   // full string
    cvec_t QualsString;   // full quality score
    ivec_t ReadsOffset;
    ivec_t QualsOffset;

    while(1){
        bool lastRead = readBatch( fasta,qual,ReadsString,ReadsOffset,
                                   QualsString,QualsOffset,inPara_);
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


void ECDriver::tableMaker() {

    // std::cout << "Constructing Facility Tables ... \n";
    /*
     * 1. split K positions into random chunks based
     *    on inPara_.k and inPara_.eSearch
     */
    ivec_t indices(inPara_.K, 0);
    for (int i = 0; i < inPara_.K; ++i) indices[i] = i;
    std::random_shuffle(indices.begin(), indices.end());
    //std::for_each(indices.begin(), indices.end(), print);  std::cout <<"\n";

    int unitSize = (log2(inPara_.eSearch) / 2);
    unsigned num = inPara_.K / unitSize;

    for (unsigned i = 0; i < num; ++i) {
        ivec_t tmp(unitSize, 0);
        std::copy(indices.begin() + i*unitSize, indices.begin() + (i + 1)*unitSize, tmp.begin());
        //tmp.insert(tmp.end(), indices.begin() + i*unitSize,
        //        indices.begin() + (i + 1) * unitSize);
        maskIdx_.push_back(tmp);
    }

    int lastUnitSize = inPara_.K % unitSize;
    if (lastUnitSize != 0) {
        ivec_t tmp;
        tmp.insert(tmp.end(), indices.begin() + num*unitSize, indices.end());
        maskIdx_.push_back(tmp);
        num++;
    }

    //std::cout << "\t" << num << " Tables to be created\n";

    /*
     * 2. Create Masks for tables
     */
    masks_.resize(num);
    for (unsigned i = 0; i < num; masks_[i] = 0xFFFFFFFF, i++);
    for (unsigned i = 0; i < num; ++i) {
        for (unsigned j = 0; j < maskIdx_[i].size(); ++j) {
            masks_[i] &= ~(1 << (2 * maskIdx_[i][j]));
            masks_[i] &= ~(1 << (2 * maskIdx_[i][j] + 1));
        }
    }

    //std::for_each(masks_.begin(), masks_.end(), printHex);
    //std::cout << "\n";

    /*
     * 3. Create Tables Then sort according to masks for each table
     */
    for (unsigned i = 0; i < num; ++ i){
        for(int j = 0; j < ecdata_.getKmerCount();j++)
            table_.insert(table_.end(), ecdata_.getKmerAt(j).ID);
    }
    /*
     * sort tables
     */
    int tableSize = ecdata_.getKmerCount();
    for (unsigned i = 0; i < num; ++i) {
        std::sort(table_.begin() + i*tableSize,
                table_.begin() + (i + 1) * tableSize, TComp(masks_[i]));
    }

    // std::cout << "\tdone !\n\n";

}

bool checkPoint(const ivec_t& indices, int dPoint) {
    for (unsigned i = 0; i < indices.size(); ++i) {
        if (indices[i] >= dPoint)
            return true;
    }
    return false;
}

//
// This unit neighbor function uses the table of masked lists
//
void ECDriver::tableUnitNeighbor(uvec_t& myNB, int dPoint) {

#ifdef DEBUG
    std::cout << "DPOINT : " <<  dPoint
              << " INPUT : " << myNB.size() << " : " ;
    for(int i = 0; i < myNB.size(); i++ )
        std::cout << myNB[i] << " " ;
    std::cout << std::endl;
#endif
    //
    if (myNB.size() == 0){
#ifdef DEBUG
        std::cout << "OUTPUT : NONE" << std::endl ;
#endif
        return;
    }

    uvec_t inputIDs(myNB.size(), 0);
    for (unsigned i = 0; i < myNB.size(); inputIDs.push_back(myNB[i]), ++ i);
    std::sort(inputIDs.begin(), inputIDs.end());

    uvec_t nbIDs; // neigbhoring IDs

    for (unsigned i = 0; i < myNB.size(); ++i) {

        uint32_t myID = myNB[i];

        uvec_t myNeighbor; //store ID's neighbors by searching in all tables

        /*
         *  search in every table
         */

        int tableSize = table_.size() / masks_.size();
        for (unsigned j = 0; j < masks_.size(); ++j) {

            if (!checkPoint(maskIdx_[j], dPoint)) continue; //a

            std::pair<uvec_t::iterator, uvec_t::iterator> bounds;

            bounds = std::equal_range(table_.begin() + j*tableSize,
                    table_.begin() + (j + 1) * tableSize, myID, TComp(masks_[j]));

            // identify if HD = 1
            for (uvec_t::iterator it = bounds.first; it != bounds.second; ++it) {
                // check if differed position satisfies dPoint
                int idx = inPara_.K - 1;
                bool dPointFlag = true;

                //if (*it == myID) continue;
                if (std::binary_search(inputIDs.begin(), inputIDs.end(), *it)) continue;

                uint32_t tmp = (*it) ^ myID;
                int cnt = 0;
                while (tmp) {

                    if (idx < dPoint) {
                        dPointFlag = false;
                        break;
                    }

                    if (tmp & 0x3) cnt++;
                    tmp >>= 2;
                }
                if (cnt == 1 && dPointFlag) {
                    myNeighbor.push_back(*it);
                }
            }
        }

        nbIDs.insert(nbIDs.begin(), myNeighbor.begin(), myNeighbor.end());
    }

    /*
     * remove duplicates
     */
    std::sort(nbIDs.begin(), nbIDs.end());
    uvec_t::iterator it
            = std::unique_copy(nbIDs.begin(), nbIDs.end(), nbIDs.begin());
    nbIDs.resize(it - nbIDs.begin());

    /*
     * now have all nbIDs from table, then search in [snKArray_]
     * to get complete multiplicity and quality score information
     */
#ifdef DEBUG
    uvec_t outnbs;
    std::cout << "OUTPUT : " ;
#endif
    for (unsigned i = 0; i < nbIDs.size(); ++i) {

        if (ecdata_.findKmer(nbIDs[i])){
#ifdef DEBUG
            outnbs.push_back(nbIDs[i]);
#endif
            myNB.push_back(nbIDs[i]);
        }
    }
#ifdef DEBUG
    for (int i = 0; i < outnbs.size(); ++i)
            std::cout << outnbs[i] << " " ;
    std::cout << std::endl;
#endif
}

/*
 * Keep elements of N only if the ID equals to the last 2*len bits in tiles
 */
void updateNodes(kcvec_t& N, kcvec_t& N_rv, const kcvec_t& tiles, int len) {

    std::set<uint64_t> IDs, IDs_rv;
    kcvec_t tmpVec;
    for (unsigned i = 0; i < tiles.size(); ++i) {
        uint64_t last2k = tiles[i].ID & ((0x1 << 2 * len) - 1);
        IDs.insert(last2k);
        IDs_rv.insert(reverse_complementary <uint32_t, uint32_t > (last2k, len));
    }
    for (unsigned i = 0; i < N.size(); ++i) {
        if (IDs.count(N[i].ID)) tmpVec.push_back(N[i]);
    }
    N = tmpVec;
    tmpVec.clear();
    for (unsigned i = 0; i < N_rv.size(); ++i) {
        if (IDs_rv.count(N_rv[i].ID)) tmpVec.push_back(N_rv[i]);
    }
    N_rv = tmpVec;
}

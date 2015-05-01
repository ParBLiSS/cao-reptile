/***
 *
 *  Author: Ankit Shah <shah.29ankit@gmail.com>
 *          Sriram P C <sriram.pc@gmail.com>
 *
 *  Copyright 2011-2012 Ankit Shah, Sriram P C, Srinivas Aluru
 *
 *  This file is part of Parallel Reptile (version 1.1)
 *
 *  PReptile is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PReptile is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <stdint.h>
#include <limits.h>

#include "util.h"
#include "ECData.hpp"


static bool goodQuality(char* qAddr, int kvalue, const Para& myPara){
    if (!qAddr) return false;
    int badpos = 0;
    for (int i = 0; i < kvalue; ++ i){
        if ((int) qAddr[i] < myPara.qualThreshold){
            badpos++;
        }
    }
    if (badpos > myPara.maxBadQPerKmer) return false;
    return true;
}


template <typename KeyIDType>
void add_kmers(char *line,char *qAddr,
               bool QFlag,int read_length,int kLength,ECData& ecdata)
{
    const Para& params = ecdata.getParams();
    int i = 0,failidx = 0;
    KeyIDType ID;
    // Insert the first kmer
    while( i < (read_length - kLength) &&
           !toID(ID, failidx, line + i, kLength)) i += failidx + 1;

    if(i < (read_length - kLength)) {
        if(QFlag == false ||
           (QFlag && goodQuality(qAddr+i,kLength,params))) {
            ecdata.addToArray(ID,1);
        }
    }

    while(i < (read_length - kLength) ) {
        KeyIDType begin = (KeyIDType) char_to_bits(line[i]);
        int end =  char_to_bits(line[i+kLength]);
        bool addKmer = false;

        // if the next char is valid, calculate using previous one
        if(end != -1) {
            ID = ((ID - (begin << (2*kLength -2))) << 2 ) + end;
            addKmer = true;
            i++;
        } else {
            // otherwise, shift to the next valid kmer
            failidx = -1; i += kLength + 1;
            do {
                i += failidx + 1;
                addKmer = toID(ID, failidx, line + i, kLength);
            } while (!addKmer && i < (read_length - kLength));
        }

        if( addKmer ) {
            if(QFlag == false ||
               (QFlag && goodQuality(qAddr+i,kLength,params)))
            {
                ecdata.addToArray(ID,1);
            }
        }
    }
}

void processBatch(cvec_t &ReadsString,cvec_t &QualsString,
                  ivec_t &ReadsOffset,ivec_t &QualsOffset,
                  ECData& ecdata)
{
    static int _batch = 0;
    const Para& params = ecdata.getParams();
    int kLength  = params.K,
        tileLength = kLength + params.step;

#ifdef DEBUG
    double tBatchStart = MPI_Wtime();
#endif
    ecdata.setBatchStart();
    //std::cout << " PROCESS: " << ecdata->m_params.mpi_env->rank()
    //         << " LOADING  BATCH " << _batch << std::endl;
    for(unsigned long i = 0; i < ReadsOffset.size();i++) {
        int position = ReadsOffset[i];
        int qposition = QualsOffset[i];
        char* addr = const_cast<char*> (&ReadsString[position]);
        char* qAddr = const_cast<char*> (&QualsString[qposition]);
        int read_length = strlen(addr);
        add_kmers<kmer_id_t>(addr,qAddr,false,read_length,kLength,ecdata);
        add_kmers<tile_id_t>(addr,qAddr,true,read_length,tileLength,ecdata);
    }
    ecdata.mergeBatchKmers();
#ifdef DEBUG
	std::stringstream out;
    	out << " PROCESS: " << ecdata.m_params.mpi_env->rank()
        << " BATCH " << _batch
        << " LOAD TIME " << MPI_Wtime()-tBatchStart << std::endl;
        std::cout << out.str();
#endif
    _batch++;

}

void processReadsFromFile(ECData& ecdata){
    // get the parameters
    const Para& params = ecdata.getParams();

    std::ifstream read_stream(params.iFaName.c_str());
    assert(read_stream.good() == true);

    read_stream.seekg(params.offsetStart,std::ios::beg);

    cvec_t ReadsString;   // full string
    cvec_t QualsString;   // full quality score
    ivec_t ReadsOffset;
    ivec_t QualsOffset;

    while(1){
        int readid = 0;
        bool lastRead = readBatch(&read_stream, params.batchSize,
                                  params.offsetEnd, ReadsString, ReadsOffset,
                                  QualsString, QualsOffset, readid);
        std::stringstream out ;
#ifdef DEBUG
        out << "PROC : " << params.mpi_env->rank()  << " "
            << ReadsOffset.size() << " " << lastRead << " QS: "
            << QualsOffset.size() << std::endl;
        std::cout << out.str();
#endif
        assert(ReadsOffset.size() == QualsOffset.size());

        processBatch(ReadsString,QualsString,
                     ReadsOffset,QualsOffset,ecdata);

        ReadsString.resize(0);
        QualsString.resize(0);
        ReadsOffset.resize(0);
        QualsOffset.resize(0);
        if(lastRead) break;
    }
}

// ------------------------------------------------------
// For map reduce :
//    map function to emit the reads
// In general :
//    emit the key-value pairs
// ------------------------------------------------------
void count_kmers(ECData& ecdata){
    const Para& params = ecdata.getParams();

    if(params.storeReads) {
        // reads are already collected and stored
        processBatch(ecdata.m_ReadsString,ecdata.m_QualsString,
                     ecdata.m_ReadsOffset,ecdata.m_QualsOffset,
                     ecdata);
    } else {
        // process reads from input file
        processReadsFromFile(ecdata);
    }
}

// Get the reads corresponding to this processor and
// store it in ECData object
bool getReadsFromFile(ECData& ecdata){

    Para& params = ecdata.getParams();

    std::ifstream read_stream(params.iFaName.c_str());
    if(!read_stream.good()) {
        std::cout << "open " << params.iFaName << "failed :|\n";
        exit(1);
    }
    read_stream.seekg(params.offsetStart,std::ios::beg);
    params.batchSize = INT_MAX; // Get all reads of this processor
    int readid = 0;
    bool lastRead = readBatch(&read_stream, params.batchSize,
                              params.offsetEnd, ecdata.m_ReadsString,
                              ecdata.m_ReadsOffset,ecdata.m_QualsString,
                              ecdata.m_QualsOffset, readid);
    return lastRead;
}

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

void processBatch(ReadStore& rbatch, ECData& ecdata)
{
    static int _batch = 0;
    const Para& params = ecdata.getParams();
    int kLength  = params.K,
        tileLength = params.tileLength;

#ifdef DEBUG
    double tBatchStart = MPI_Wtime();
#endif
    ecdata.setBatchStart();
    //std::cout << " PROCESS: " << ecdata->m_params.mpi_env->rank()
    //         << " LOADING  BATCH " << _batch << std::endl;
    for(unsigned long i = 0; i < rbatch.readsOffset.size();i++) {
        int position = rbatch.readsOffset[i];
        int qposition = rbatch.qualsOffset[i];
        char* addr = const_cast<char*> (&(rbatch.readsString[position]));
        char* qAddr = const_cast<char*> (&(rbatch.qualsString[qposition]));
        int read_length = strlen(addr);
        add_kmers<kmer_id_t>(addr,qAddr,false,read_length,kLength,ecdata);
        add_kmers<tile_id_t>(addr,qAddr,true,read_length,tileLength,ecdata);
    }
    ecdata.mergeBatch();
#ifdef DEBUG
	std::stringstream out;
    	out << " PROCESS: " << ecdata.m_params.m_rank
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

    ReadStore rbatch;
    while(1){
        bool lastRead = readBatch(&read_stream, params.batchSize,
                                  params.offsetEnd, rbatch);
#ifdef DEBUG
        std::stringstream out ;
        out << "PROC : " << params.m_rank  << " "
            << rbatch.readsOffset.size() << " " << lastRead << " QS: "
            << rbatch.qualsOffset.size() << std::endl;
        std::cout << out.str();
#endif
        assert(rbatch.readsOffset.size() == rbatch.qualsOffset.size());

        processBatch(rbatch, ecdata);

        rbatch.reset();
        if(lastRead) break;
    }
}

void count_kmers(ECData& ecdata){
    const Para& params = ecdata.getParams();

    if(params.storeReads) {
        // reads are already collected and stored
        processBatch(ecdata.m_reads, ecdata);
    } else {
        // process reads from input file
        processReadsFromFile(ecdata);
    }
}


template <typename KmerType>
void kmerSpectrumBatch(ReadStore& rbatch, ECData& ecdata,
                      int kLength, bool QFlag)
{
    static int _batch = 0;

#ifdef DEBUG
    const Para& params = ecdata.getParams();
    double tBatchStart = MPI_Wtime();
#endif
    KmerType tmp = 0;
    ecdata.setBatchStart(tmp);
    assert(tmp == 1);
    //std::cout << " PROCESS: " << ecdata->m_params.mpi_env->rank()
    //         << " LOADING  BATCH " << _batch << std::endl;
    for(unsigned long i = 0; i < rbatch.readsOffset.size();i++) {
        int position = rbatch.readsOffset[i];
        int qposition = rbatch.qualsOffset[i];
        char* addr = const_cast<char*> (&(rbatch.readsString[position]));
        char* qAddr = const_cast<char*> (&(rbatch.qualsString[qposition]));
        int read_length = strlen(addr);
        add_kmers<KmerType>(addr,qAddr,QFlag,read_length,kLength,ecdata);
    }
    tmp = 0;
    ecdata.mergeBatch(tmp);
    assert(tmp == 1);
#ifdef DEBUG
	std::stringstream out;
    	out << " PROCESS: " << ecdata.m_params.m_rank
        << " BATCH " << _batch
        << " LOAD TIME " << MPI_Wtime()-tBatchStart << std::endl;
        std::cout << out.str();
#endif
    _batch++;

}

template <typename KmerType>
void kmerSpectrumFile(ECData& ecdata, int kLength, bool qFlag){
    // get the parameters
    const Para& params = ecdata.getParams();

    std::ifstream read_stream(params.iFaName.c_str());
    assert(read_stream.good() == true);

    read_stream.seekg(params.offsetStart,std::ios::beg);

    ReadStore rbatch;
    while(1){
        bool lastRead = readBatch(&read_stream, params.batchSize,
                                  params.offsetEnd, rbatch);
#ifdef DEBUG
        std::stringstream out ;
        out << "PROC : " << params.m_rank  << " "
            << rbatch.readsOffset.size() << " " << lastRead << " QS: "
            << rbatch.qualsOffset.size() << std::endl;
        std::cout << out.str();
#endif
        assert(rbatch.readsOffset.size() == rbatch.qualsOffset.size());

        kmerSpectrumBatch<KmerType>(rbatch, ecdata, kLength, qFlag);

        rbatch.reset();
        if(lastRead) break;
    }
}

void local_kmer_spectrum(ECData& ecdata){
    const Para& params = ecdata.getParams();

    if(params.storeReads) {
        // reads are already collected and stored
        kmerSpectrumBatch<kmer_id_t>(ecdata.m_reads, ecdata,
                                    params.K, false);
    } else {
        // process reads from input file
        kmerSpectrumFile<kmer_id_t>(ecdata, params.K, false);
    }
}

void local_tile_spectrum(ECData& ecdata){
    const Para& params = ecdata.getParams();

    if(params.storeReads) {
        // reads are already collected and stored
        kmerSpectrumBatch<tile_id_t>(ecdata.m_reads, ecdata,
                                    params.tileLength, true);
    } else {
        // process reads from input file
        kmerSpectrumFile<tile_id_t>(ecdata, params.tileLength, true);
    }
}

/***
 **
 *  File: util.h
 *  Created: Dec 12, 2009 4:05 PM
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *  Copyright 2011-12 Ankit Shah, Sriram P C, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1).
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

#ifndef _UTIL_H
#define	_UTIL_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <complex> //pow
#include <ctype.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <stdint.h>

#include <sys/time.h>
#include <MPI_env.hpp>
#include <cassert>
#include "fasta_file.hpp"

typedef std::vector<int> ivec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::vector<uvec_t> uuvec_t;
typedef std::vector<char> cvec_t;
typedef std::vector<bool> bvec_t;
typedef std::vector<std::string> strvec_t;
typedef std::map<int, int> imap_t;
typedef std::set<int> iset_t;
typedef std::pair<int, int> ipair_t;
typedef std::pair<uint32_t, uint32_t> upair_t;
typedef uint32_t kmer_id_t;
typedef uint64_t tile_id_t;
//typedef MPI_UNSIGNED_LONG_LONG mpi_kmer_id_t;
#define mpi_kmer_id_t MPI_UNSIGNED // MPI_UNSIGNED_LONG_LONG
#define mpi_tile_id_t MPI_UNSIGNED_LONG_LONG // MPI_UNSIGNED_LONG_LONG

class Para {
public:
    std::string iFaName;
    int QFlag;
    std::string iQName;
    int batchSize;  //specify number of reads to be loaded
    int K;
    std::string oErrName;
    int step;
    int qualThreshold;
    int Qlb;
    int maxBadQPerKmer;
    int eSearch;
    int tGoodTile;
    int tCard;
    //int tConst;
    int hdMax;
    double tRatio;
    int storeReads;
    bool absentKmers;
    // cache additions
    int kmerCacheSize;
    int tileCacheSize;
    int cacheOptimizedSearch;

    Para(const char *configFile) {
        setPara(configFile);
        load_parallel_params();
    };

    virtual ~Para(){
        delete mpi_env;
    };

private:

    void setPara(const char *configFile) {
        cacheOptimizedSearch = 0;
        absentKmers = false;
        std::string line, s1;
        std::ifstream input(configFile);
        std::istringstream buf;
        storeReads = 0;
        kmerCacheSize = tileCacheSize = 4;
        while(getline(input, line)){
            buf.clear();
            buf.str(line);

            if (buf >> s1){
                if (s1 == "InFaFile")
                    buf >> iFaName;
                else if (s1 == "QFlag")
                    buf >> QFlag;
                else if (s1 == "IQFile")
                    buf >> iQName;
                else if (s1 == "OErrFile")
                    buf >> oErrName;
                else if (s1 == "BatchSize")
                    buf >> batchSize;
                else if (s1 == "KmerLen")
                    buf >> K;
                else if (s1 == "hd_max")
                    buf >> hdMax;
                else if (s1 == "Step")
                    buf >> step;
                else if (s1 == "ExpectSearch")
                    buf >> eSearch;
                else if (s1 == "T_ratio")
                    buf >> tRatio;
                else if (s1 == "QThreshold")
                    buf >> qualThreshold;
                else if (s1 == "MaxBadQPerKmer")
                    buf >> maxBadQPerKmer;
                else if (s1 == "Qlb")
                    buf >> Qlb;
                else if (s1 == "T_expGoodCnt")
                    buf >> tGoodTile;
                else if (s1 == "T_card")
                    buf >> tCard;
                else if (s1 == "StoreReads")
                    buf >> storeReads;
                else if (s1 == "CacheOptimizedSearch")
                    buf >> cacheOptimizedSearch;
                else if (s1 == "KmerCacheSize")
                    buf >> kmerCacheSize;
                else if (s1 == "TileCacheSize")
                    buf >> tileCacheSize;
            }
        }

    }

  public:
    empi::MPI_env *mpi_env;
    long readsPerProcessor;
    unsigned long startFromLineNo;
    unsigned long endTillLineNo;
    unsigned long offsetStart;
    unsigned long qOffsetStart;
    long inLastBatch;

    unsigned long long offsetEnd;
    unsigned long long qOffsetEnd;
    unsigned long long fileSize;

    void load_parallel_params() {
        mpi_env = new empi::MPI_env();
        compute_parallelio_params();
    };

    bool validate() {
        if (iFaName.empty()) {
            if(mpi_env->rank() == 0)
                std::cout << "Err: InFaFile is not specified!\n";
            return false;
        }
        if (oErrName.empty()) {
            if(mpi_env->rank() == 0)
                std::cout << "Err: Output file is not specified!\n";
            return false;
        }
        if (K > 16 || (K + step) > 32) {
            if(mpi_env->rank() == 0)
                std::cout <<
                    "Set K in the range of (0, 16] and K+step in the range of (2, 32]\n";
            return false;
        }
        std::ifstream read_stream(iFaName.c_str());
        if(!read_stream.good()) {
            std::cout << "open " << iFaName << "failed :|\n";
            return false;
        }
        std::ifstream qual_stream(iQName.c_str());
        if(!qual_stream.good()) {
            std::cout << "open " << iQName << "failed :|\n";
            return false;
        }

        if(mpi_env->rank() == 0) {
            std::cout << "Input Parameters:\n----------------------------------\n";
            std::cout << "Short Reads file: " << iFaName << "\n";
            if (QFlag){
                std::cout << "I/QualFile = " << iQName << "\n"
                          << "Qthreshold = " << qualThreshold << "\n";
            }
            std::cout << "O/ErrFile = " << oErrName << "\n";
            std::cout << "(K, step, tile) = " << "(" << K << "," << step << ","
                      << K + step << ")\n"
                      << "BatchSize = " << batchSize << "\n"
                      << "Max Hamming Distance Allowed = " << hdMax << "\n"
                      << "ExpectSearch = " << eSearch << "\n"
                      << "T_ratio = " << tRatio << "\n"
                      << "QThreshold = " << qualThreshold << "\n"
                      << "MaxBadQPerKmer = " << maxBadQPerKmer << "\n"
                      << "Qlb = " << Qlb << "\n"
                      << "T_expGoodCnt = " << tGoodTile << "\n"
                      << "T_card = " << tCard << "\n"
                      << "StoreReads = " << storeReads << "\n"
                      << "CacheOptimizedSearch = " << cacheOptimizedSearch << "\n"
                      << "Kmer Cache = " << kmerCacheSize << "\n"
                      << "Tile Cache = " << tileCacheSize << "\n"
                      << "----------------------------------\n";
        }
        return true;
    }
/*
    void compute_readcount(){
        int rank = mpi_env->rank(),
            size = mpi_env->size();

        // TODO: Should do this more elegantly
        if(rank == 0){
            std::stringstream out;
            std::string tmpfile = "./nolines";
            // generate at tem file
            out << "wc -l " << this->iFaName
                << " | cut -d' ' -f1 > " << tmpfile;
            system(out.str().c_str());
            std::ifstream fin(tmpfile.c_str());
            std::string line;
            std::istringstream buf;

            std::getline(fin,line);
            buf.clear(); buf.str(line);
            buf >> readCount;
            assert(readCount > (size-1));
    #ifdef DEBUG
            std::cout << "NO OF LINES IS " << readCount << std::endl;
    #endif
            readCount /= 2;
            fin.close();
        }
        MPI_Bcast (&readCount, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    }
*/
   void compute_parallelio_params(){
       compute_offsets();
    }

   void compute_offsets(){
       std::ifstream fin(iFaName.c_str());
       fin.seekg(0,std::ios::end);
       fileSize = fin.tellg();
       offsetStart = (fileSize/mpi_env->size()) * mpi_env->rank();
       offsetEnd = (fileSize/mpi_env->size()) * (mpi_env->rank()+1);
       fin.seekg(offsetStart, std::ios::beg);


       bIO::FASTA_input fasta_sr(fin);
       ++fasta_sr;
       std::string lineNo = (*fasta_sr).first;

       fin.seekg(offsetEnd, std::ios::beg);
       bIO::FASTA_input fastaend(fin);
       ++fastaend;
       std::string endLineNo = (*fastaend).first;

       startFromLineNo = strtoul (lineNo.c_str(),NULL,0);
       endTillLineNo = strtoul (endLineNo.c_str(),NULL,0) - 1;
       qOffsetStart = offsetStart + startFromLineNo ;
       qOffsetEnd = offsetEnd + endTillLineNo ;

       // Print the qual size for debug purposes
       std::ifstream qfin(iQName.c_str());
       qfin.seekg(0,std::ios::end);
       long long qfileSize = qfin.tellg();

       std::string curline,nextline;
       qfin.seekg(qOffsetStart,std::ios::beg);
       if(qfin.good()) {
           std::getline(qfin,curline);
           if(qfin.good()) {
               std::getline(qfin,nextline);
               if( curline[0] == '>' && nextline[0] == '>') {
                   qOffsetStart += curline.size();
               }
           }
       }

       qfin.seekg(qOffsetStart,std::ios::beg);
       bIO::FASTA_input qfasta_sr(qfin);
       ++qfasta_sr;
       std::string qlineNo = (*qfasta_sr).first;


       qfin.seekg(qOffsetEnd,std::ios::beg);
       bIO::FASTA_input qfastaend(qfin);
       ++qfastaend;
       std::string qendLineNo = (*qfastaend).first;
       long quendTillLineNo = strtoul (qendLineNo.c_str(),NULL,0) - 1;

       if(readsPerProcessor % batchSize == 0)
           inLastBatch  = batchSize;
       else
           inLastBatch = readsPerProcessor % batchSize;
#ifdef DEBUG
       std::stringstream out;
       std::string testStr;
       unsigned long qualstartFromLineNo = strtoul (qlineNo.c_str(),NULL,0);

       out << "PROC : " << mpi_env->rank() << " FILE SIZE : "
           <<  fileSize <<   " START FROM : "
           << startFromLineNo << " START OFFSET : " << offsetStart
           << " UP TO : " << endTillLineNo << std::endl;
       out << "PROC : " << mpi_env->rank() << " QUAL FILE SIZE : "
           <<  qfileSize << " QUAL START FROM : "
           << qlineNo  << " START OFFSET : " << qOffsetStart
           << " UP TO : " << quendTillLineNo
           << (( qualstartFromLineNo ==  startFromLineNo ) ? " TRUE" : " FALSE")
           << std::endl;
       std::cout << out.str();
#endif
   }

/*
    void compute_parallelio_params(int fileRecordLength,
                    long& perprocessor,long& startfrom,long &offset, long& inLastBatch) {
        int rank = mpi_env->rank(),
            size = mpi_env->size();
        unsigned long n = readCount;

        perprocessor = n/size;
        short temp = n%size;
        startfrom = rank * perprocessor;

        if(rank < temp){
            startfrom += rank;
            perprocessor += 1;
        }
        else
            startfrom += temp;
        if(perprocessor % batchSize == 0)
            inLastBatch  = batchSize;
        else
            inLastBatch = perprocessor % batchSize;

        offset = (4+fileRecordLength) * startfrom + 3;
        long added = 9;
        temp = startfrom - added;

        while(temp > 0 ){
            offset += temp;
            added *= 10;
            temp -= added;
        }
    }
*/

};



inline int char_to_bits(char c) {
    static bool flag = false;
    int cvalue = -1;

    switch (c) {
        case 'A':
        case 'a':
            cvalue = 0;
            break;
        case 'C':
        case 'c':
            cvalue = 1;
            break;
        case 'G':
        case 'g':
            cvalue = 2;
            break;
        case 'T':
        case 't':
            cvalue = 3;
            break;
    }

    if (cvalue != -1) {
        return cvalue;
    } else {
        if(flag == false){
                //std::cout << "Note: non-ACGT character found: " << c
                //          << " , all will be ignored\n";
                flag = true;
        }
        return -1;
    }
}
inline char bits_to_char(int value) {

    switch (value) {
        case 0:
            return 'a';
        case 1:
            return 'c';
        case 2:
            return 'g';
        case 3:
            return 't';
        default:
            return 'n';
    }
}

/*
 *  reverse_complementary test code:
 *  uint32_t test = 314324;
 *  std::cout << toString(test, 10) << "\n";
 *  uint32_t result = reverse_complementary<uint32_t, uint32_t> (test, 10);
 *  std::cout << toString(result, 10) << "\n";
 *
*/
template <typename Tin, typename Tout>
Tout reverse_complementary (Tin num, int len){
    Tout rv = 0;
    for (int i = 0; i < len; ++ i){
        int cvalue = num & 0x3;
        switch (cvalue){
            case 0:
                rv = (rv << 2) | 0x3;
                break;
            case 1:
                rv = (rv << 2) | 0x2;
                break;
            case 2:
                rv = (rv << 2) | 0x1;
                break;
            case 3:
                rv = rv << 2;
                break;
        }
        num >>= 2;
    }
    return rv;
}

template <typename T>
inline bool toID(T& ID, char* addr, int len) {
    ID = 0;
    for (int i = 0; i < len; ++ i){
        int c = char_to_bits(addr[i]);
        if (c == -1) return false;
	ID  = (ID << 2 | c);
    }
    return true;
}

template <typename T>
inline bool toID(T& ID, int &failidx, char* addr, int len) {
    ID = 0; failidx = 0;
    for (int i = 0; i < len; ++ i){
        int c = char_to_bits(addr[i]);
        if (c == -1) {
            failidx = i;
            return false;
        }
        ID  = (ID << 2 | c);
    }
    return true;
}


inline std::string toString(kmer_id_t ID, int len){

    std::string kmer = "";
    for (int i = 0; i < len; ++ i){
        int last = (ID & 0x3);
        char c;
        switch (last){
            case 0: c = 'a';
            break;
            case 1: c = 'c';
            break;
            case 2: c = 'g';
            break;
            case 3: c = 't';
            break;
        }
        kmer += c;
        ID = ID >> 2;
    }
    std::reverse(kmer.begin(), kmer.end());

    return kmer;
}

inline double get_time() {
      timeval t;
      gettimeofday(&t, 0);
      return t.tv_sec + (0.000001 * t.tv_usec);
} // get_time

inline void print_time (const std::string& msg, double& timing){
    double cur_time = get_time();
    std::cout << msg << "(" << cur_time - timing << " secs)\n\n";
    timing = cur_time;
}

//inline void readBatch(FILE *filep,FILE *fileqp,int rounds, int currentRound,
//               cvec_t &ReadsString,cvec_t &QualsString, const Para &myPara)
inline bool readBatch(bIO::FASTA_input& fasta,bIO::FASTA_input& qual,
                      cvec_t &ReadsString,ivec_t &ReadsOffset,
                      cvec_t &QualsString,ivec_t &QualsOffset,const Para &myPara)
{
    unsigned long position = 0;
    typedef bIO::FASTA_input::value_type value_type;


    bool lastRead = false;
    unsigned long curLine = 0;
    for(long j=0; j < myPara.batchSize ; j++){

        const value_type& v = *fasta;

        std::string posstr = v.first;
        std::string read_str = v.second;
        curLine = strtoul(v.first.c_str(),NULL,0);

        if( myPara.mpi_env->rank() < myPara.mpi_env->size() - 1 &&
            curLine > myPara.endTillLineNo)  {
            lastRead = true;
            break;
        }
        int read_length = read_str.length();

        ReadsString.resize(ReadsString.size() + read_length + 1);
        memcpy(&ReadsString[position], read_str.c_str(), read_length + 1);
        ReadsOffset.push_back(position);

        position += read_length + 1;
        if(++fasta == false){
            lastRead = true;
            break;
        }
    }
#ifdef DEBUG
    std::stringstream out;
    out << "READ PROC : " << myPara.mpi_env->rank() << " " << curLine << std::endl;
    std::cout << out.str();
#endif
    position = 0;
    for(long j = 0; j < myPara.batchSize; j++){
        const value_type& v = *qual;

        std::string posstr = v.first;
        std::string quals_str = v.second;
        quals_str.erase(0,1);
        int read_length = quals_str.length();

        unsigned long curLine = strtoul(v.first.c_str(),NULL,0);
        if( myPara.mpi_env->rank() < myPara.mpi_env->size() - 1 &&
            curLine > myPara.endTillLineNo)
            break;

        QualsString.resize(QualsString.size() + read_length + 1);
        memcpy(&QualsString[position], quals_str.c_str(), read_length + 1);
        QualsOffset.push_back(position);

        position += read_length + 1;
        if(++qual == false) break;
    }
    return lastRead;
}

#endif	/* _UTIL_H */

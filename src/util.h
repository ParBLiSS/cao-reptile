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
#include <vector>
#include <string>
#include <utility>
#include <cstdint>

typedef std::vector<int> ivec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::vector<uvec_t> uuvec_t;
typedef std::vector<char> cvec_t;
typedef std::vector<bool> bvec_t;
typedef std::vector<std::string> strvec_t;

typedef std::pair<int, int> ipair_t;
typedef std::pair<uint32_t, uint32_t> upair_t;
typedef uint32_t kmer_id_t;
typedef uint64_t tile_id_t;

#define mpi_kmer_id_t MPI_UNSIGNED // MPI_UNSIGNED_LONG_LONG
#define mpi_tile_id_t MPI_UNSIGNED_LONG_LONG // MPI_UNSIGNED_LONG_LONG
/// macros for block decomposition
#define BLOCK_LOW(i,p,n) ( (i*n) / p)
#define BLOCK_HIGH(i,p,n) ( (((i+1)*n)/p) - 1)
#define BLOCK_SIZE(i,p,n) (BLOCK_LOW((i+1),p,n) - BLOCK_LOW(i,p,n))
#define BLOCK_OWNER(j,p,n) (((p) * ((j)+1)-1)/(n))


class Para {
public:
    std::string iFaName;
    std::string oErrName;
    std::string outputFilename;
    double tRatio;
    bool QFlag;
    bool absentKmers;
    int batchSize;  //specify number of reads to be loaded
    int K;
    int step;
    int qualThreshold;
    int Qlb;
    int maxBadQPerKmer;
    int eSearch;
    int tGoodTile;
    int tCard;
    int hdMax;
    int storeReads;

    // cache optimizations
    int kmerCacheSize;
    int tileCacheSize;
    int cacheOptimizedSearch;
    int writeOutput;

    int writeSpectrum;
    int numThreads;
    std::string kmerSpectrumOutFile;
    std::string tileSpectrumOutFile;

    // size and rank of the MPI world
    int m_size;
    int m_rank;

    // input file offsets
    unsigned long offsetStart;
    unsigned long offsetEnd;
    unsigned long fileSize;

    Para(const char *configFile);
    bool validate();

 private:

    void setPara(const char *configFile);
    void load_parallel_params();

    void compute_offsets();
};

struct ReadStore {
    cvec_t readsString;
    ivec_t readsOffset;
    cvec_t qualsString;
    ivec_t qualsOffset;
    int readId;

    void reset(){
        readsString.resize(0);
        readsOffset.resize(0);
        qualsString.resize(0);
        qualsOffset.resize(0);
    }

    void swap(ReadStore& other){
        std::swap(readsString, other.readsString);
        std::swap(qualsString, other.qualsString);
        std::swap(readsOffset, other.readsOffset);
        std::swap(qualsOffset, other.qualsOffset);
        std::swap(readId, other.readId);
    }
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

// UTILITY FUNCTIONS -------------------------------------------------
// trim taken from stack overflow
// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
            s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

std::string toString(kmer_id_t ID, int len);
double get_time();
void print_time (const std::string& msg, double& timing);
bool readBatch(std::ifstream* fqfs, const long& batchSize,
               const unsigned long& offsetEnd,
               ReadStore& rbatch);
bool goodQuals(const char* qAddr, int len, int threshold);
#endif	/* _UTIL_H */

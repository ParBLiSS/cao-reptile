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

#ifndef _ECDATA_H
#define _ECDATA_H

#include <vector>
#include <string>
#include <cstdint>
#include <climits>
#include <mpi.h>

#include "util.h"
#include "flat_layout.hpp"
#include "caware_layout.hpp"
#include "coblivious_layout.hpp"

// Structure from Reptile
typedef struct KC{
    uint64_t ID;
    int goodCnt;
    int cnt;
    KC(){};
    KC(uint64_t myID, int c1, int c2): ID(myID), goodCnt(c1), cnt(c2){};
}kc_t;

typedef std::vector<kc_t> kcvec_t;

struct Knumcomp {
    bool operator() (const kc_t& e1, const kc_t& e2) const {
        return  (e1.ID < e2.ID);
    }
};


typedef struct kmer_s{
    kmer_id_t ID;
    unsigned char count;
} kmer_t;


typedef struct tile_s{
    tile_id_t ID;
    unsigned char count;
} tile_t;


template <typename StructDataType>
struct StructComp {
    bool operator()(const StructDataType& k1, const StructDataType& k2) const {
        return k1.ID < k2.ID;
    }
};

template <typename KeyDataType>
struct KeyIDComp {
    bool operator()(const KeyDataType& k1, const KeyDataType& k2) const {
        return k1 < k2;
    }
};

typedef StructComp<tile_t> TileComp;
typedef KeyIDComp<tile_id_t> TileIDComp;
typedef StructComp<kmer_t> KmerComp;
typedef KeyIDComp<kmer_id_t> KmerIDComp;
typedef std::vector<kmer_id_t> kmer_id_vector;
typedef std::vector<tile_id_t> tile_id_vector;
typedef std::vector<unsigned char> kcount_vector;
static const unsigned MAX_LEVELS = 64;

class ECData {
  private:
    void registerKmerTypes();
    ECDataCOLayout<kmer_t*, kmer_id_t,
                   unsigned char> m_kmerCOLayout;
    ECDataCOLayout<tile_t*, tile_id_t,
                   unsigned char> m_tileCOLayout;

    ECDataCALayout<kmer_t*, kmer_id_t,
                   unsigned char> m_kmerCALayout;
    ECDataCALayout<tile_t*, tile_id_t,
                   unsigned char> m_tileCALayout;

    ECDataFlatLayout<kmer_t*, kmer_id_t,
                     unsigned char> m_kmerFlatLayout;
    ECDataFlatLayout<tile_t*, tile_id_t,
                     unsigned char> m_tileFlatLayout;

    ParECDataCALayout<kmer_t*, kmer_id_t,
                      unsigned char> m_kmerParCALayout;
    ParECDataCALayout<tile_t*, tile_id_t,
                      unsigned char> m_tileParCALayout;

    // I own karray, ksize and kcount. So, Please be nice to them!
    kmer_t *m_karray;
    long m_ksize, m_kcount;
    // I own tilearray, tilesize and tilecount. So, Please be nice to them!
    tile_t *m_tilearray;
    long m_tilesize, m_tilecount;

    // I am only pointing to this parameter object. I don't own it!
    Para& m_params;

    MPI_Datatype m_mpi_kmer_t;
    MPI_Datatype m_mpi_tile_t;
    mutable uint64_t m_kmerQueries;
    mutable uint64_t m_kmerQueryFails;
    mutable uint64_t m_tileQueries;
    mutable uint64_t m_tileQueryFails;
    mutable unsigned m_kmerLevels[MAX_LEVELS];
    mutable unsigned m_tileLevels[MAX_LEVELS];

    // Store Reads if reqd.
    ReadStore m_reads;
    std::vector<std::string> m_fullReads;
    std::vector<std::string> m_fullQuals;

    std::vector<kmer_id_t> m_byte_kref[3];
    unsigned m_byte_kcount[3];
    std::vector<tile_id_t> m_byte_tref[7];
    unsigned m_byte_tcount[7];

    int currentKmerBatchStart;
    int currentTileBatchStart;

    void padKmerArray(unsigned kSize);
    void padTileArray(unsigned kSize);
    void estimateKmerByteCounters();
    void estimateTileByteCounters();

    void buildCacheAwareLayout(const unsigned& kmerCacheSize,
                               const unsigned& tileCacheSize);
    void buildParCacheAwareLayout(const unsigned& kmerCacheSize,
                                  const unsigned& tileCacheSize);
    void buildCacheObliviousLayout();
    void buildFlatLayout();

    bool findKmerDefault(const kmer_id_t &kmerID) const;
    bool findKmerFlat(const kmer_id_t &kmerID) const;
    bool findKmerCacheAware(const kmer_id_t &kmerID) const;
    bool findKmerCacheOblivious(const kmer_id_t &kmerID) const;
    bool findKmerParCacheAware(const kmer_id_t &kmerID) const;

    long findTileFlat(const tile_id_t &tileID,kc_t& output) const;
    long findTileDefault(const tile_id_t &tileID,kc_t& output) const;
    long findTileCacheAware(const tile_id_t &tileID,kc_t& output) const;
    long findTileCacheOblivious(const tile_id_t &tileID,kc_t& output) const;
    long findTileParCacheAware(const tile_id_t &tileID,kc_t& output) const;

  public:
    friend void sort_kmers(ECData& ecdata);
    friend void count_kmers(ECData& ecdata);
    friend void gather_spectrum(ECData& ecdata);
    friend void local_tile_spectrum(ECData& ecdata, bool qFlag);
    friend void local_kmer_spectrum(ECData& ecdata, bool qFlag);
    friend void dist_kmer_spectrum(ECData& ecdata);
    friend void dist_tile_spectrum(ECData& ecdata);

    Para& getParams() const{return m_params;}

    long getKmerCount() const{return m_kcount;}
    long getTileCount() const{return m_tilecount;}
    const kmer_t& getKmerAt(int j) const{return m_karray[j];}
    const cvec_t& getReads() const{return m_reads.readsString;}
    const cvec_t& getQuals() const{return m_reads.qualsString;}
    const ivec_t& getReadsOffsets() const{return m_reads.readsOffset;}
    const ivec_t& getQualsOffsets() const{return m_reads.qualsOffset;}
    const ReadStore& getReadStore() const {return m_reads;}

    uint64_t getKmerQueries() const {return m_kmerQueries;}
    uint64_t getKmerQueryFails() const {return m_kmerQueryFails;}
    uint64_t getTileQueries() const {return m_tileQueries;}
    uint64_t getTileQueryFails() const {return m_tileQueryFails;}
    unsigned getKmerLevels(const int& j) const {
        return j < (int)MAX_LEVELS ? m_kmerLevels[j] : 0;
    }
    unsigned getTileLevels(const int& j) const {
        return j < (int)MAX_LEVELS ? m_tileLevels[j]: 0;
    }

    bool getReadsFromFile(unsigned& nreads);
    bool getAllReads(long nReads);
    bool loadReads(unsigned long woffStart, unsigned long woffEnd, ReadStore& rbatch) const;
    long getFullReadSize() const{
      return m_fullReads.size();
    }
    bool addToArray(kmer_id_t &ID,int count);
    bool addToArray(kmer_id_t &ID,unsigned char count);
    bool addToArray(tile_id_t &ID,int count);

    void replaceKArray(kmer_t* allData,int allDataCount,int allDataSize);
    void replaceTileArray(tile_t *newTileArray,int newTileCount,
                          int newTileSize);

    void mergeBatch();
    void mergeBatch(kmer_id_t &tmp);
    void mergeBatch(tile_id_t &tmp);
    void setBatchStart();
    void setBatchStart(kmer_id_t &tmp);
    void setBatchStart(tile_id_t &tmp);
    void buildCacheOptimizedLayout();
    void estimateByteCounters();
    void runCAStats(std::ostream& ots);

    bool findKmer(const kmer_id_t &kmerID) const;
    long findTile(const tile_id_t &tileID,kc_t& output) const;

    void printByteCounters(std::ostream& ots) const;
    void printKArray() const;
    void writeQueryStats(const std::string& filename) const;
    void writeSpectrum() const;

    void loadSpectrum();
    void writeKmerDistSpectrum() const;
    void writeTileDistSpectrum() const;
    void writeDistSpectrum() const;

    void makeTileHistogram(long tileHist[UCHAR_MAX + 1]);

    ECData(Para& p);
    virtual ~ECData();
};

#endif

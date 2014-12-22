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
#define MAX_LEVELS 64
class ECData {
  private:
    void registerKmerTypes();
  public:
    // I own karray, ksize and kcount. So, Please be nice to them!
    kmer_t *m_karray;
    int m_ksize, m_kcount;
    // I own tilearray, tilesize and tilecount. So, Please be nice to them!
    tile_t *m_tilearray;
    int m_tilesize, m_tilecount;

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

    MPI_Datatype m_mpi_kmer_t;
    MPI_Datatype m_mpi_tile_t;
    uint64_t m_kmerQueries;
    uint64_t m_kmerQueryFails;
    uint64_t m_tileQueries;
    uint64_t m_tileQueryFails;
    unsigned m_kmerLevels[MAX_LEVELS];
    unsigned m_tileLevels[MAX_LEVELS];


    // Store Reads if reqd.
    cvec_t m_ReadsString;
    cvec_t m_QualsString;
    ivec_t m_ReadsOffset;
    ivec_t m_QualsOffset;

    int currentKmerBatchStart;
    int currentTileBatchStart;

    std::vector<kmer_id_t> m_byte_kref[3];
    unsigned m_byte_kcount[3];
    std::vector<tile_id_t> m_byte_tref[3];
    unsigned m_byte_tcount[3];

    void padKmerArray(unsigned kSize);
    void padTileArray(unsigned kSize);
    void estimateKmerByteCounters();
    void estimateTileByteCounters();
    void estimateByteCounters();
    void printByteCounters(std::ostream& ots);
    void runCAStats(std::ostream& ots);

  public:
    // I am only pointing to this parameter object. I don't own it!
    Para *m_params;

    bool addToArray(kmer_id_t &ID,int count);
    bool addToArray(kmer_id_t &ID,unsigned char count);
    bool addToArray(tile_id_t &ID,int count);

    void printKArray();
    void replaceKArray(kmer_t* allData,int allDataCount,int allDataSize);
    void replaceTileArray(tile_t *newTileArray,int newTileCount,
                          int newTileSize);
    bool findKmer(const kmer_id_t &kmerID);
    bool findKmerDefault(const kmer_id_t &kmerID);
    bool findKmerFlat(const kmer_id_t &kmerID);
    bool findKmerCacheAware(const kmer_id_t &kmerID);
    bool findKmerCacheOblivious(const kmer_id_t &kmerID);

    int findTile(const tile_id_t &tileID,kc_t& output);
    int findTileFlat(const tile_id_t &tileID,kc_t& output);
    int findTileDefault(const tile_id_t &tileID,kc_t& output);
    int findTileCacheAware(const tile_id_t &tileID,kc_t& output);
    int findTileCacheOblivious(const tile_id_t &tileID,kc_t& output);

    void mergeBatchKmers();
    void setBatchStart();
    void buildCacheOptimizedLayout();
    void buildCacheAwareLayout(const unsigned& kmerCacheSize,
                               const unsigned& tileCacheSize);

    void buildCacheObliviousLayout();
    void buildFlatLayout();
    void output(const std::string& filename);
    void writeSpectrum();

    ECData(Para *p);
    virtual ~ECData();
};

#endif

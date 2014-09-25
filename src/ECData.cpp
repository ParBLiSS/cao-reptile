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

#include "util.h"
#include "ECData.hpp"
#include "sort_kmers.hpp"
#include "caware_layout.hpp"
#include "implicit_heap_search_c.h"
#include <cassert>
#include <climits>

ECData::ECData(Para *p) : m_karray(0),m_ksize(0),m_kcount(0),
        m_tilearray(0),m_tilesize(0),m_tilecount(0),m_params(p){
    registerKmerTypes();
}

void ECData::registerKmerTypes(){
    MPI_Datatype type[2] = { mpi_kmer_id_t, MPI_UNSIGNED_CHAR };
    int blocklen[2];
    MPI_Aint disp[2];

    blocklen[0] = 1;
    blocklen[1] = 1;
    disp[0] = 0;
    disp[1] = sizeof(kmer_id_t);

    MPI_Type_create_struct(2, blocklen, disp, type, &m_mpi_kmer_t);
    MPI_Type_commit(&m_mpi_kmer_t);

    type[0] =  mpi_tile_id_t;
    type[1] = MPI_UNSIGNED_CHAR;
    blocklen[0] = 1;
    blocklen[1] = 1;
    disp[0] = 0;
    disp[1] = sizeof(tile_id_t);

    MPI_Type_create_struct(2, blocklen, disp, type, &m_mpi_tile_t);
    MPI_Type_commit(&m_mpi_tile_t);


}

bool ECData::findKmer(const kmer_id_t &kmerID) {
   if(m_karray == 0){
       return false;
   }
   kmer_t searchKmer;
   searchKmer.ID = kmerID;
   bool result = std::binary_search(m_karray,m_karray + m_kcount,
                                    searchKmer, KmerComp());
   if(m_params->absentKmers == true)
       return !result;
   return result;
}

bool ECData::findKmerCacheAware(const kmer_id_t &kmerID) {
    kmer_id_vector::iterator result =
        pure_c::implicit_heap_search(m_kmer_ID.begin(),
                                     m_kmer_ID.end(),
                                     m_kdegree,
                                     kmerID);

    bool final = false;

   if(m_params->absentKmers == true)
       final = (result == m_kmer_ID.end());
   else
       final = (result != m_kmer_ID.end());
    // std::cout << "Find " << kmerID << " : "
    //           << (m_params->absentKmers) << " : "
    //           << final << std::endl;
    return final;
}

int ECData::findTile(const tile_id_t &tileID,kc_t& output) {
    int lb = 0, ub = m_tilecount - 1, mid;
    while (lb <= ub) {
        mid = (lb + ub) / 2;
        if (m_tilearray[mid].ID == tileID) {
            output.ID = m_tilearray[mid].ID;
            output.goodCnt = m_tilearray[mid].count;
            output.cnt = m_tilearray[mid].count;
            return mid;
        }
        else if (m_tilearray[mid].ID < tileID)
            lb = mid + 1;
        else if (m_tilearray[mid].ID > tileID)
            ub = mid - 1;
    }
    return -1;
}

int ECData::findTileCacheAware(const tile_id_t &tileID,kc_t& output){
    tile_id_vector::iterator result =
        pure_c::implicit_heap_search(m_tile_ID.begin(),
                                     m_tile_ID.end(),
                                     m_tdegree,
                                     tileID);
    int final;
    if( result != m_tile_ID.end()) {
        int mid = result - m_tile_ID.begin();
        output.ID = tileID;
        output.goodCnt = output.cnt = m_tile_count[mid];
        final =  mid;
    } else{
        final = -1;
    }
    // std::cout << "Find " << tileID << " : "
    //           << ((final >= 0) ? m_tile_count[final] : 0) <<  std::endl;
    return final;
}


int ECData::getKmerCount(kmer_id_t &ID){
    // TODO: do a binary search in kmer array
    return -1;
}


bool ECData::addToArray(kmer_id_t &ID,int count){
    if(m_karray == 0){
        m_karray = (kmer_t*) malloc(sizeof(kmer_t)*2);
        m_ksize = 2;
    }
    else if(m_kcount+1 > m_ksize){
        m_karray = (kmer_t*) realloc(m_karray, sizeof(kmer_t) * (2* m_ksize));
        m_ksize *= 2;
    }
    assert(m_karray != NULL);

    m_karray[m_kcount].ID = ID;
    m_karray[m_kcount].count = (count < UCHAR_MAX) ? count : UCHAR_MAX;
    m_kcount++;
    return true;
}

bool ECData::addToArray(kmer_id_t &ID,unsigned char count){
    if(m_karray == 0){
        m_karray = (kmer_t*) malloc(sizeof(kmer_t)*2);
        m_ksize = 2;
    }
    else if(m_kcount+1 > m_ksize){
        m_karray = (kmer_t*) realloc(m_karray, sizeof(kmer_t) * (2* m_ksize));
        m_ksize *= 2;
    }
    assert(m_karray != NULL);

    m_karray[m_kcount].ID = ID;
    m_karray[m_kcount].count = (count < UCHAR_MAX) ? count : UCHAR_MAX;
    m_kcount++;
    return true;
}

bool ECData::addToArray(tile_id_t &ID,int count){
    if(m_tilearray == 0){
        m_tilearray = (tile_t*) malloc(sizeof(tile_t)*2);
        m_tilesize = 2;
    }
    else if(m_tilecount+1 > m_tilesize){
        m_tilearray = (tile_t*) realloc(m_tilearray, sizeof(tile_t) * (2* m_tilesize));
        m_tilesize *= 2;
    }
    assert(m_tilearray != NULL);

    m_tilearray[m_tilecount].ID = ID;
    m_tilearray[m_tilecount].count =  (count < UCHAR_MAX) ? count : UCHAR_MAX;
    m_tilecount++;
    return true;
}

void ECData::printKArray(){

    std::cout << "Karray " <<
        m_params->mpi_env->rank() << m_kcount << std::endl;

    for(int i = 0; i < m_kcount;i++){
        // just to print
        std::stringstream out;

        out << m_params->mpi_env->rank() << " " << m_karray[i].count << " "
            << m_karray[i].ID << std::endl;
        std::cout << out.str();
    }
}

void ECData::replaceKArray(kmer_t* allData,int allDataCount,int allDataSize){
    free(m_karray);
    m_karray = allData; m_kcount = allDataCount; m_ksize = allDataSize;
}

void ECData::replaceTileArray(tile_t *newTileArray,int newTileCount,
                              int newTileSize) {
    free(m_tilearray);
    m_tilearray = newTileArray; m_tilecount = newTileCount;
    m_tilesize = newTileSize;
}

ECData::~ECData(){
    free(m_karray);
    free(m_tilearray);
}

void ECData::setBatchStart(){
    currentKmerBatchStart = m_kcount;
    currentTileBatchStart = m_tilecount;
}

void ECData::mergeBatchKmers(){

    // Sort the new ones
    std::sort(m_karray + currentKmerBatchStart,
              m_karray + m_kcount, KmerComp());

    // Merge 2 arrays
    if(currentKmerBatchStart > 0) {
        kmer_t *newarray = (kmer_t*) malloc(sizeof(kmer_t) * m_ksize);
        kmer_t *mergeLast = std::merge(m_karray, m_karray + currentKmerBatchStart,
                                       m_karray + currentKmerBatchStart, m_karray + m_kcount,
                                       newarray, KmerComp());
        // for(int i = 0; i < m_kcount;i++) {
        //     std::cout << " " << newarray[i].ID << std::endl;
        // }
        int newcount = (mergeLast - newarray);
        free(m_karray);
        m_kcount = newcount;
        m_karray = newarray;
    }

    // Eliminiate dupes.
    eliminate_dupes(m_karray,m_kcount);

    // Update tiles
    std::sort(m_tilearray + currentTileBatchStart,
              m_tilearray + m_tilecount, TileComp());

    if(currentTileBatchStart > 0) {
        tile_t *newarray = (tile_t*)malloc(sizeof(tile_t) * m_tilesize);
        tile_t *merge_tile_lst = std::merge(m_tilearray, m_tilearray + currentTileBatchStart,
                                            m_tilearray + currentTileBatchStart,
                                            m_tilearray + m_tilecount,
                                            newarray, TileComp());
        int newtilecount = (merge_tile_lst - newarray);
        free(m_tilearray);
        m_tilecount = newtilecount;
        m_tilearray = newarray;
    }


    eliminate_dupes(m_tilearray,m_tilecount);
    // reset this
    currentKmerBatchStart = 0; currentTileBatchStart = 0;
}

void ECData::buildCacheAwareStructure(const unsigned& kmerCacheSize,
                                      const unsigned& tileCacheSize){
    //
    m_kdegree = kmerCacheSize + 1;
    unsigned rSize = (m_kcount % kmerCacheSize);
    unsigned fillIn = kmerCacheSize - rSize;
    if(rSize > 0){
        kmer_t& lastKmer = m_karray[m_kcount - 1];
        for(unsigned i = 0; i < fillIn; i++){
            addToArray(lastKmer.ID, lastKmer.count);
        }
    }
    m_kmer_ID.resize(m_kcount);
    m_kmer_count.resize(m_kcount);

    caware_layout(m_karray, m_karray + m_kcount,
                  kmerCacheSize, m_kmer_ID, m_kmer_count);
    free(m_karray);
    m_karray = 0;
    m_kcount = 0;

    m_tdegree = tileCacheSize + 1;
    rSize = (m_tilecount % tileCacheSize);
    fillIn = tileCacheSize - rSize;
    if(rSize > 0){
        tile_t& lastTile = m_tilearray[m_tilecount - 1];
        for(unsigned i = 0; i < fillIn; i++){
            addToArray(lastTile.ID, lastTile.count);
        }
    }
    m_tile_ID.resize(m_tilecount);
    m_tile_count.resize(m_tilecount);

    caware_layout(m_tilearray, m_tilearray + m_tilecount,
                  tileCacheSize, m_tile_ID, m_tile_count);
    free(m_tilearray);
    m_tilearray = 0;
    m_tilecount = 0;
}

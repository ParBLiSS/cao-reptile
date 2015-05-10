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
#include <fstream>

template<typename T>
int count_bytes(T val){
    int n = 0;
    if(val == 0)
        return 1;
    while (val != 0) {
        val >>= 8;
        n ++;
    }
    return n;
}

ECData::ECData(Para& p) : m_karray(0),m_ksize(0),m_kcount(0),
        m_tilearray(0),m_tilesize(0),m_tilecount(0),m_params(p){
    m_kmerQueries = 0;
    m_kmerQueryFails = 0;
    m_tileQueries = 0;
    m_tileQueryFails = 0;
    for(unsigned i = 0; i < MAX_LEVELS; i++)
        m_tileLevels[i] = 0;
    for(unsigned i = 0; i < MAX_LEVELS; i++)
        m_kmerLevels[i] = 0;
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

// Get the reads corresponding to this processor and
// store it in ECData object
bool ECData::getReadsFromFile(){
    std::ifstream read_stream(m_params.iFaName.c_str());
    if(!read_stream.good()) {
        std::cout << "open " << m_params.iFaName << "failed :|\n";
        exit(1);
    }
    read_stream.seekg(m_params.offsetStart, std::ios::beg);
    m_params.batchSize = INT_MAX; // Get all reads of this processor

    bool lastRead = readBatch(&read_stream, m_params.batchSize,
                              m_params.offsetEnd, m_reads);
    unsigned nsize = m_reads.size();
    unsigned tmp = nsize;
    MPI_Reduce( &tmp, &nsize, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
    if(m_params.m_rank == 0){
      std::cout << "nreads\t" << nsize << std::endl;
    }
    return lastRead;
}


bool ECData::findKmerDefault(const kmer_id_t &kmerID) const{
  if(m_karray == 0){
    return false;
  }
  bool final = false;
  int lb = 0, ub = m_kcount  - 1, mid;
#ifdef QUERY_COUNTS
  int nlevels = 0;
#endif

  while (lb <= ub) {
      mid = (lb + ub) / 2;
      if (m_karray[mid].ID == kmerID) {
          final = true;
          break;
      }
      else if (m_karray[mid].ID < kmerID)
          lb = mid + 1;
      else if (m_karray[mid].ID > kmerID)
          ub = mid - 1;
#ifdef QUERY_COUNTS
      nlevels += 1;
#endif
  }
  //kmer_t searchKmer;
  //searchKmer.ID = kmerID;
  //bool final = std::binary_search(m_karray,m_karray + m_kcount,
  //                                searchKmer, KmerComp());

#ifdef QUERY_COUNTS
  if(final){
      assert(nlevels < (int)MAX_LEVELS);
      m_kmerLevels[nlevels] += 1;
  }
  m_kmerQueries += 1;
  if(!final)
      m_kmerQueryFails += 1;
#endif


  if(m_params.absentKmers == true)
    final = !final;

  // std::cout << "Find " << kmerID << " : "
  //           << (m_params.absentKmers) << " : "
  //           << final << std::endl;

  return final;
}

bool ECData::findKmerFlat(const kmer_id_t &kmerID) const {
   bool final = m_kmerFlatLayout.find(kmerID);
   if(m_params.absentKmers == true)
       final = !final;
   // std::cout << "Find " << kmerID << " : "
   //           << (m_params.absentKmers) << " : "
   //           << final << std::endl;

   return final;
}

bool ECData::findKmerCacheAware(const kmer_id_t &kmerID) const{

    bool final = m_kmerCALayout.find(kmerID);

   if(m_params.absentKmers == true)
       final = !final;

    // std::cout << "Find " << kmerID << " : "
    //            << (m_params.absentKmers) << " : "
    //            << final << std::endl;
    return final;
}

bool ECData::findKmerCacheOblivious(const kmer_id_t &kmerID) const{
    bool final = m_kmerCOLayout.find(kmerID);
   if(m_params.absentKmers == true)
       final = !final;

    // std::cout << "Find " << kmerID << " : "
    //            << (m_params.absentKmers) << " : "
    //            << final << std::endl;
    return final;
}

bool ECData::findKmer(const kmer_id_t &kmerID) const{
    switch(m_params.cacheOptimizedSearch){
        case 1:
            return findKmerCacheAware(kmerID);
        case 2:
            return findKmerCacheOblivious(kmerID);
        case 3:
            return findKmerFlat(kmerID);
        default:
            return findKmerDefault(kmerID);
    }
}

int ECData::findTileDefault(const tile_id_t &tileID,kc_t& output) const{
  int lb = 0, ub = m_tilecount - 1, mid;
  int final = -1;
#ifdef QUERY_COUNTS
  int nlevels = 0;
#endif
  while (lb <= ub) {
    mid = (lb + ub) / 2;
    if (m_tilearray[mid].ID == tileID) {
      output.ID = m_tilearray[mid].ID;
      output.goodCnt = m_tilearray[mid].count;
      output.cnt = m_tilearray[mid].count;
      final = mid;
      break;
    }
    else if (m_tilearray[mid].ID < tileID)
      lb = mid + 1;
    else if (m_tilearray[mid].ID > tileID)
      ub = mid - 1;
#ifdef QUERY_COUNTS
     nlevels += 1;
#endif
  }
#ifdef QUERY_COUNTS
  if(final){
      assert(nlevels < (int)MAX_LEVELS);
      m_tileLevels[nlevels] += 1;
  }
  m_tileQueries += 1;
  if(final == -1)
      m_tileQueryFails += 1;
#endif
  // std::cout << "Find Tile " << tileID << " : "
  //           << ((final >= 0) ? output.cnt : 0) <<  std::endl;
  return final;
}


int ECData::findTileFlat(const tile_id_t &tileID,kc_t& output) const{
    unsigned char count = 0;
    int final = m_tileFlatLayout.getCount(tileID, count);
    if( final != -1 ) {
        output.ID = tileID;
        output.goodCnt = output.cnt = count;
    }

    return final;
}

int ECData::findTileCacheAware(const tile_id_t &tileID,kc_t& output) const{
    unsigned char count = 0;
    int final = m_tileCALayout.getCount(tileID, count);

    if( final != -1 ) {
        output.ID = tileID;
        output.goodCnt = output.cnt = count;
    }
    // std::cout << "Find Tile " << tileID << " : "
    //           << ((final >= 0) ? count : 0) <<  std::endl;
    return final;
}

int ECData::findTileCacheOblivious(const tile_id_t &tileID,kc_t& output) const{
    unsigned char count = 0;
    int final = m_tileCOLayout.getCount(tileID, count);
    if( final != -1 ) {
        output.ID = tileID;
        output.goodCnt = output.cnt = count;
    }
    // std::cout << "Find Tile " << tileID << " : "
    //           << ((final >= 0) ? count : 0) <<  std::endl;
    return final;
}

int ECData::findTile(const tile_id_t &tileID,kc_t& output) const{
   switch(m_params.cacheOptimizedSearch){
        case 1:
            return findTileCacheAware(tileID, output);
        case 2:
            return findTileCacheOblivious(tileID, output);
        case 3:
            return findTileFlat(tileID, output);
        default:
            return findTileDefault(tileID, output);
    }
}

static const long size_factor = 10;
static const long size_threshold = 536870912;
long get_size(const long& cursize){
  if(cursize < size_threshold)
    return 2 * cursize;
  return (cursize) + (cursize / size_factor);
}

bool ECData::addToArray(kmer_id_t &ID,int count){
    if(m_karray == 0){
        m_karray = (kmer_t*) malloc(sizeof(kmer_t)*2);
        m_ksize = 2;
    }
    else if(m_kcount+1 > m_ksize){
        long tl = get_size(m_ksize);
        m_karray = (kmer_t*) realloc(m_karray, sizeof(kmer_t) * tl);
        m_ksize = tl;
    }
    if(m_karray == NULL){
      std::stringstream out1;
      out1 << m_params.m_rank << " " 
           << m_ksize << " " << m_kcount << " "
           << m_tilesize << " " << m_tilecount << " " << std::endl;
      std::cout << out1.str();
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
        long tl = get_size(m_ksize);
        m_karray = (kmer_t*) realloc(m_karray, sizeof(kmer_t) * tl);
        m_ksize = tl;
    }
    if(m_karray == NULL){
      std::stringstream out1;
      out1 << m_params.m_rank << " " 
           << m_ksize << " " << m_kcount << " "
           << m_tilesize << " " << m_tilecount << " " << std::endl;
      std::cout << out1.str();
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
        long tl = get_size(m_ksize);
        m_tilearray = (tile_t*) realloc(m_tilearray, sizeof(tile_t) * tl);
        m_tilesize = tl;
    }
    if(m_tilearray == NULL){
      std::stringstream out1;
      out1 << m_params.m_rank << " " 
           << m_ksize << " " << m_kcount << " "
           << m_tilesize << " " << m_tilecount << " " << std::endl;
      std::cout << out1.str();
    }
    assert(m_tilearray != NULL);

    m_tilearray[m_tilecount].ID = ID;
    m_tilearray[m_tilecount].count =  (count < UCHAR_MAX) ? count : UCHAR_MAX;
    m_tilecount++;
    return true;
}

void ECData::printKArray() const{

    std::cout << "Karray " <<
        m_params.m_rank << m_kcount << std::endl;

    for(int i = 0; i < m_kcount;i++){
        // just to print
        std::stringstream out;

        out << m_params.m_rank << " " << m_karray[i].count << " "
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

void ECData::writeDistSpectrum() const{
    if(m_kcount == 0 || m_tilecount == 0)
        return;
    std::ofstream kFile (m_params.kmerSpectrumOutFile.c_str(),
                         std::ofstream::binary);
    std::ofstream tFile (m_params.tileSpectrumOutFile.c_str(),
                         std::ofstream::binary);
    kFile.write((char*)m_karray, sizeof(kmer_t)*m_kcount);
    tFile.write((char*)m_tilearray, sizeof(tile_t)*m_tilecount);

    kFile.close();
    tFile.close();
}

void ECData::loadSpectrum(){
    std::ifstream kFile(m_params.kmerSpectrumInFile.c_str(),
                        std::ifstream::binary);
    std::ifstream tileFile(m_params.tileSpectrumInFile.c_str(),
                           std::ifstream::binary);
    if(kFile){
        kFile.seekg(0, kFile.end);
        unsigned fLength = kFile.tellg();
        kFile.seekg(0, kFile.beg);
        if(m_kcount > 0)
            free(m_karray);

        m_kcount = fLength / sizeof(kmer_t);
        m_karray = (kmer_t*) malloc(sizeof(kmer_t) * m_kcount);
        kFile.read((char*) m_karray,  sizeof(kmer_t) * m_kcount);
        m_ksize = m_kcount;
    }
    if(tileFile){
        tileFile.seekg(0, tileFile.end);
        unsigned fLength = tileFile.tellg();
        tileFile.seekg(0, tileFile.beg);
        if(m_tilecount > 0)
            free(m_tilearray);
        m_tilecount = fLength / sizeof(tile_t);
        m_tilearray = (tile_t*) malloc(sizeof(tile_t) * m_tilecount);
        tileFile.read((char*) m_tilearray,  sizeof(tile_t) * m_tilecount);
        m_tilesize = m_tilecount;
    }
}

void ECData::writeSpectrum() const{
    if(m_params.runType != 0){ // don't write when run type is otherwise
        return;
    }
    if(m_params.writeSpectrum > 0 && (m_params.m_rank == 0)){
        assert(m_params.kmerSpectrumOutFile.length() > 0);
        assert(m_params.tileSpectrumOutFile.length() > 0);

        std::ofstream kFile (m_params.kmerSpectrumOutFile.c_str(),
                             std::ofstream::binary);
        std::ofstream tFile (m_params.tileSpectrumOutFile.c_str(), std::ofstream::binary);
        switch(m_params.cacheOptimizedSearch){
        case 1:
            m_kmerCALayout.serialize(kFile);
            m_tileCALayout.serialize(tFile);
            break;
        case 2:
            m_kmerCOLayout.serialize(kFile);
            m_tileCOLayout.serialize(tFile);
            break;
        case 3:
            m_kmerFlatLayout.serialize(kFile);
            m_tileFlatLayout.serialize(tFile);
            break;
        default:
            kFile.write((char*)m_karray, sizeof(kmer_t)*m_kcount);
            tFile.write((char*)m_tilearray, sizeof(tile_t)*m_tilecount);
            break;
        }
        kFile.close();
        tFile.close();
    }
}

ECData::~ECData(){
    free(m_karray);
    free(m_tilearray);
}


void ECData::setBatchStart(tile_id_t &tmp){
    currentKmerBatchStart = m_kcount;
    tmp = 1;
}

void ECData::setBatchStart(kmer_id_t &tmp){
    currentTileBatchStart = m_tilecount;
    tmp = 1;
}

void ECData::setBatchStart(){
    kmer_id_t tk;
    tile_id_t tl;
    setBatchStart(tk);
    setBatchStart(tl);
}

void ECData::mergeBatch(){
    kmer_id_t tk;
    tile_id_t tl;
    mergeBatch(tk);
    mergeBatch(tl);
}

void ECData::mergeBatch(tile_id_t &tmp){
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
    currentTileBatchStart = 0;
    tmp = 1;
}

void ECData::mergeBatch(kmer_id_t &tmp){

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

    // reset this batch
    currentKmerBatchStart = 0;
    tmp = 1;
}

void ECData::padKmerArray(unsigned kSize){
    unsigned rSize = (m_kcount % kSize);
    unsigned fillIn = kSize - rSize;
    if(rSize > 0){
        if(m_params.m_rank == 0) {
            std::stringstream out;
            out << "Padding kmer array : " << fillIn << std::endl;
            std::cout << out.str();
        }
        kmer_t lastKmer = m_karray[m_kcount - 1];
        for(unsigned i = 0; i < fillIn; i++){
            addToArray(lastKmer.ID, lastKmer.count);
        }
    }
}

void ECData::padTileArray(unsigned kSize){
    unsigned rSize = (m_tilecount % kSize);
    unsigned fillIn = kSize - rSize;
    if(rSize > 0){
        if(m_params.m_rank == 0) {
            std::stringstream out;
            out << "Padding tile array : " << fillIn << std::endl;
            std::cout << out.str();
        }
        tile_t lastTile = m_tilearray[m_tilecount - 1];
        for(unsigned i = 0; i < fillIn; i++){
            addToArray(lastTile.ID, lastTile.count);
        }
    }
}

void ECData::buildCacheAwareLayout(const unsigned& kmerCacheSize,
                                   const unsigned& tileCacheSize){
    int rank = m_params.m_rank;
    if(rank == 0) {
       std::stringstream out;
       out << "Build Kmer Cache Aware Layout : " << m_kcount
                 << " Kmer Cache : " << kmerCacheSize << std::endl;
       std::cout << out.str();
    }

    padKmerArray(kmerCacheSize);
    m_kmerCALayout.init(m_karray, m_kcount, kmerCacheSize);

    free(m_karray); m_karray = 0; m_kcount = 0;

    if(rank == 0){
       std::stringstream out;
       out << "Build Tile Cache Aware Layout : " << m_tilecount
                 << " Tile Cache : " << tileCacheSize << std::endl;
       std::cout << out.str();
    }

    padTileArray(tileCacheSize);
    m_tileCALayout.init(m_tilearray, m_tilecount, tileCacheSize);

    free(m_tilearray); m_tilearray = 0; m_tilecount = 0;
}

void ECData::buildCacheObliviousLayout(){
    int rank = m_params.m_rank;

    padKmerArray(1024);
    if(rank == 0){
       std::stringstream out;
       out << "Build Kmer Cache Oblivious Layout : "
                 << m_kcount << std::endl;
       std::cout << out.str();
    }
    m_kmerCOLayout.init(m_karray, m_kcount, 20);

    padTileArray(1024);
    if(rank == 0) {
       std::stringstream out;
       out << "Build Tile Cache Oblivious Layout : "
                 << m_tilecount << std::endl;
       std::cout << out.str();
    }
    m_tileCOLayout.init(m_tilearray, m_tilecount, 20);
}

void ECData::buildFlatLayout(){
    int rank = m_params.m_rank;
    if(rank == 0) {
       std::stringstream out;
       out << "Build Kmer Flat Layout : "
                 << m_kcount << std::endl;
       std::cout << out.str();
    }
    m_kmerFlatLayout.init(m_karray, m_kcount);
    free(m_karray);
    m_karray = 0;
    m_kcount = 0;
    if(rank == 0){
       std::stringstream out;
       out << "Build Tile Flat Layout : "
                 << m_tilecount << std::endl;
       std::cout << out.str();
    }
    m_tileFlatLayout.init(m_tilearray, m_tilecount);
    free(m_tilearray);
    m_tilearray = 0;
    m_tilecount = 0;
}

void ECData::buildCacheOptimizedLayout(){

   switch(m_params.cacheOptimizedSearch){
       case 1:
           buildCacheAwareLayout(m_params.kmerCacheSize,
                                 m_params.tileCacheSize);
           break;
       case 2:
           buildCacheObliviousLayout();
           break;
       case 3:
           buildFlatLayout();
           break;
       default:
           break;
   }
}

void ECData::writeQueryStats(const std::string& filename) const{
    std::ofstream oHandle(filename.c_str(), ios::out | ios::app );
    if (!oHandle.good()) {
        std::cout << "open " << filename << " failed, correct path?\n";
    return;
    }
    oHandle << m_kmerQueryFails << "/"
        << m_kmerQueries << std::endl;
    oHandle << m_tileQueryFails << "/"
        << m_tileQueries << std::endl;
    oHandle.close();
}

void ECData::estimateKmerByteCounters(){
    for(int j = 0; j < 3; j++){
        m_byte_kref[j].push_back(m_karray[0].ID);
        m_byte_kcount[j] = 0;
    }

    for(int i = 1; i < m_kcount;i++){
        std::stringstream out;
        out << i << " " << m_karray[i].ID << " ";
        for(int j = 0; j < 3; j++){
            kmer_id_t& last_ref = m_byte_kref[j].back();
            int n = count_bytes<kmer_id_t>(m_karray[i].ID - last_ref);
            if(n > j + 1) {
                m_byte_kref[j].push_back(m_karray[i].ID);
            } else {
                m_byte_kcount[j] += 1;
            }
            out << j + 1 << " " << last_ref << " " << (m_karray[i].ID - last_ref)
                << " " << n << " ";
        }

        out << std::endl;
        //std::cout << out.str();
    }
}

void ECData::estimateTileByteCounters(){
    for(int j = 0; j < 7; j++) {
        m_byte_tref[j].push_back(m_tilearray[0].ID);
        m_byte_tcount[j] = 0;
    }

    for(int i = 1; i < m_tilecount;i++){
        for(int j = 0; j < 7; j++){
            tile_id_t& last_ref = m_byte_tref[j].back();
            int n = count_bytes<tile_id_t>(m_tilearray[i].ID - last_ref);
            if(n > j + 1) {
                m_byte_tref[j].push_back(m_tilearray[i].ID);
            } else {
              m_byte_tcount[j] += 1;
            }
        }
    }
}

void ECData::estimateByteCounters(){
   if(m_params.cacheOptimizedSearch != 0)
       return;
    estimateKmerByteCounters();
}

void ECData::printByteCounters(std::ostream& ots) const{
   if(m_params.cacheOptimizedSearch != 0)
       return;

    std::stringstream oss;
    oss << "type" << "\t" << "nbyte" << "\t"
        << "nref" << "\t" << "ncount" << std::endl;
    ots << oss.str();
    ots.flush();

    for(int j = 0; j < 3; j++){
        std::stringstream oss2;
        oss2 << "kmer" << "\t" << (j+1) << "\t"
             << m_byte_kref[j].size() << "\t"
             << m_byte_kcount[j] << std::endl;
        ots << oss2.str();
        ots.flush();
    }

    for(int j = 0; j < 7; j++){
        std::stringstream oss2;
        oss2 << "tile" << "\t" << (j+1) << "\t"
              << m_byte_tref[j].size() << "\t"
             << m_byte_tcount[j] << std::endl;
        ots << oss2.str();
        ots.flush();
    }
}

void ECData::runCAStats(std::ostream& ots){
   if(m_params.cacheOptimizedSearch != 1)
       return;

    double kavg, kmax, kmin;
    double tavg, tmax, tmin;
    m_kmerCALayout.compression_stats(kavg, kmin, kmax);
    m_tileCALayout.compression_stats(tavg, tmin, tmax);

    std::stringstream oss;
    oss << "type" << "\t" << "avg" << "\t"
        << "max" << "\t" << "min" << std::endl;

    oss << "kmer" << "\t" << kavg << "\t"
        << kmax << "\t" << kmin << std::endl;

    oss << "tile" << "\t" << tavg << "\t"
        << tmax << "\t" << tmin << std::endl;
    ots << oss.str();
    ots.flush();

}

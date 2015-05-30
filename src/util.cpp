#include "util.h"

#include <cassert>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <iomanip>

#include <sys/time.h>
#include <mpi.h>

Para::Para(const char *configFile){
    setPara(configFile);
    load_parallel_params();
}

void Para::setPara(const char *configFile) {
    cacheOptimizedSearch = 0;
    absentKmers = false;
    std::string line, s1;
    std::ifstream input(configFile);
    std::istringstream buf;
    storeReads = 0;
    kmerCacheSize = tileCacheSize = 4;
    writeOutput = 1;
    writeSpectrum = 0;
    numThreads = 1;
    dynamicWorkDist = 0;
    runType = 0;
    workFactor = 2;
    QFlag = true;
    kinAbsent = false;
    tileLength = 0;
    while(getline(input, line)){
        buf.clear();
        buf.str(line);

        if (buf >> s1){
            if (s1 == "InFaFile")
                buf >> iFaName;
            else if (s1 == "QFlag")
                buf >> QFlag;
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
            else if (s1 == "WriteOutput")
                buf >> writeOutput;
            else if (s1 == "WriteSpectrum")
                buf >> writeSpectrum;
            else if (s1 == "KmerSpectrumOutFile")
                buf >> kmerSpectrumOutFile;
            else if (s1 == "TileSpectrumOutFile")
                buf >> tileSpectrumOutFile;
            else if (s1 == "Threads")
                buf >> numThreads;
            else if (s1 == "WorkFactor")
                buf >> workFactor;
            else if (s1 == "WorkDistribution")
                buf >> dynamicWorkDist;
            else if (s1 == "KmerSpectrumInFile")
                buf >> kmerSpectrumInFile;
            else if (s1 == "TileSpectrumInFile")
                buf >> tileSpectrumInFile;
            else if (s1 == "RunType")
                buf >> runType;
            else if (s1 == "KmerSpectrumInAbsent")
                buf >> kinAbsent;
        }
    }

    tileLength = K + step;
}

bool Para::validate() {
    if (iFaName.empty()) {
        if(m_rank == 0)
            std::cout << "Err: InFaFile is not specified!\n";
        return false;
    }
    if (oErrName.empty()) {
        if(m_rank == 0)
            std::cout << "Err: Output file is not specified!\n";
        return false;
    }
    if(writeOutput != 0 && (runType == 0 || runType == 2)){
        //  each process validates its ouput path
        std::stringstream outs;
        outs << oErrName << "-" << std::setfill('0') << std::setw(5) << m_rank ;
        outputFilename = outs.str();
        std::ofstream oHandle(outputFilename.c_str());
        if (!oHandle.good()) {
            std::cout << "err output : open " << outputFilename
                      << " failed, correct path?\n";
            return false;
        }
        oHandle.close();
    }

    if ((K < 0) || (K > 16) || (step < 1) ||
        (tileLength < K) || (tileLength > 32)) {
        if(m_rank == 0)
            std::cout <<
                "Set K in the range of (0, 16] and K+step in the range of (K, 32]\n";
        return false;
    }
    std::ifstream read_stream(iFaName.c_str());
    if(!read_stream.good()) {
        std::cout << "open " << iFaName << "failed :|\n";
        return false;
    }
    if(runType == 0){
        if(writeSpectrum != 0 && m_rank == 0)  // validate spectrum output file
            if(!validateSpectrumOutput())
                return false;
    }
    if(runType == 1){ // validate spectrum output file
        if(!validateSpectrumOutput())
            return false;
    }
    if(runType == 2){
        if(kmerSpectrumInFile.length() == 0 || tileSpectrumInFile.length() == 0){
            std::cout << "kmer sectrum input or tile spectrum input is "
                      << "mandatory for EC only option" << std::endl;
            return false;
        }
        std::ifstream kFile(kmerSpectrumInFile.c_str(),
                            std::ifstream::binary);
        if(!kFile.good()){
            std::cout << "kmer sectrum input : open " << kmerSpectrumInFile
                      << " failed, correct path?\n" << std::endl;
            return false;
        }
        kFile.close();
        std::ifstream tileFile(tileSpectrumInFile.c_str(),
                               std::ifstream::binary);
        if(!tileFile.good()){
            std::cout << "tile sectrum input : open " << kmerSpectrumInFile
                      << " failed, correct path?\n" << std::endl;
            return false;
        }
        tileFile.close();
        absentKmers = kinAbsent;
    }
    if(m_rank == 0) {
        print();
    }

    return true;
}

void Para::print(){
        std::stringstream oss;
        oss << "--" << std::endl
            << "parameter" << "\t" << "value" << std::endl;
        oss << "short reads file\t" << iFaName << "\n";
        oss << "O/ErrFile\t" << oErrName << "\n";
        oss << "(K, step, tile)\t"
            << "(" << K << "," << step << "," << tileLength << ")\n"
            << "BatchSize\t" << batchSize << "\n"
            << "Max Hamming Dist\t" << hdMax << "\n"
            << "ExpectSearch\t" << eSearch << "\n"
            << "T_ratio\t" << tRatio << "\n"
            << "QThreshold\t" << qualThreshold << "\n"
            << "MaxBadQPerKmer\t" << maxBadQPerKmer << "\n"
            << "Qlb\t" << Qlb << "\n"
            << "T_expGoodCnt\t" << tGoodTile << "\n"
            << "T_card\t" << tCard << "\n"
            << "StoreReads\t" << storeReads << "\n"
            << "CacheOptimizedSearch\t" << cacheOptimizedSearch << "\n"
            << "Kmer Cache\t" << kmerCacheSize << "\n"
            << "Tile Cache\t" << tileCacheSize << "\n"
            << "Threads \t" << numThreads << "\n"
            << "Dynamic Distribution \t" << dynamicWorkDist << "\n"
            << "Work Factor \t" << workFactor << "\n"
            << "Write Output \t" << writeOutput << "\n"
            << "Write Spectrum \t" << writeSpectrum << "\n"
            << "Kmer Spectrum Output \t" << kmerSpectrumOutFile << "\n"
            << "Tile Spectrum Output \t" << tileSpectrumOutFile << "\n"
            << "Kmer Spectrum Input \t" << kmerSpectrumInFile << "\n"
            << "Tile Spectrum Input \t" << tileSpectrumInFile << "\n"
            << "Kmer Input Absent ? \t" << kinAbsent << "\n"
            << "--\n";
        std::cout << oss.str();
        std::cout.flush();
}

void Para::load_parallel_params(){
    m_size = MPI::COMM_WORLD.Get_size();
    m_rank = MPI::COMM_WORLD.Get_rank();
    compute_offsets();
}

void Para::compute_offsets(){
    // offset computation for static work load
    std::ifstream fin(iFaName.c_str());
    fin.seekg(0,std::ios::end);
    fileSize = fin.tellg();
    offsetStart = (fileSize/m_size) * m_rank;
    offsetEnd = (fileSize/m_size) * (m_rank + 1);
}

bool Para::validateSpectrumOutput(){
    if(kmerSpectrumOutFile.length() == 0)
        kmerSpectrumOutFile = "kmer-spectrum";
    if(tileSpectrumOutFile.length() == 0)
        tileSpectrumOutFile = "tile-spectrum";
    //  each process validates its ouput path
    std::stringstream outs, outs2;
    outs << kmerSpectrumOutFile << "-" << std::setfill('0')
         << std::setw(5) << m_rank << ".bin" ;
    kmerSpectrumOutFile = outs.str();
    outs2 << tileSpectrumOutFile << "-" << std::setfill('0')
          << std::setw(5) << m_rank << ".bin" ;
    tileSpectrumOutFile = outs2.str();
    std::ofstream oHandle(kmerSpectrumOutFile.c_str());
    if (!oHandle.good()) {
        std::cout << "kmer spectrum input : open " << kmerSpectrumOutFile
                  << " failed, correct path?\n";
        return false;
    }
    oHandle.close();
    std::ofstream oHandle2(tileSpectrumOutFile.c_str());
    if (!oHandle2.good()) {
        std::cout << "tile spectrum input : open " << tileSpectrumOutFile
                  << " failed, correct path?\n";
        return false;
    }
    oHandle2.close();
    return true;
}

std::string toString(kmer_id_t ID, int len){

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

bool readFastqRecord(std::ifstream* fqfs,
                     std::string fRecord[4])
{
  for(int i = 0; i < 4; i++) {
    fRecord[i] = "";
  }
  for(int i = 0; i < 4; i++) {
    if(!fqfs->good())
      return false;
    std::getline(*fqfs, fRecord[i]);
    fRecord[i] = trim(fRecord[i]);
    if(fRecord[i].length() == 0)
      return false;
  }

  return true;
}

bool isAlphabet(const char& x){
  return (x == 'A' || x == 'a' ||
          x == 'G' || x == 'g' ||
          x == 'C' || x == 'c' ||
          x == 'T' || x == 't' ||
          x == 'N' || x == 'n');
}

bool readFirstFastqRecord(std::ifstream* fqfs,
                          std::string fRecord[4])
{
  // read strings
  for(int i = 0; i < 4; i++)
    fRecord[i] = "";

  for(int i = 0; i < 4; i++) {
    if(!fqfs->good())
      return false;
    std::getline(*fqfs, fRecord[i]);
    fRecord[i] = trim(fRecord[i]);
  }
  for(int i = 1; i < 4; i++)
    if(fRecord[i].length() == 0)
      return false;

  int sdelta = 4, sread = 4;
  // unless, the record is good, we can't turst the first line's first char
  if(fRecord[2][0] == '+' && isAlphabet(fRecord[1][0]) && fRecord[0][0] == '@'){
    return true; // record is good!
  } else if(fRecord[1][0] == '@' && isAlphabet(fRecord[2][0]) && fRecord[3][0] == '+'){
    sdelta = 1; sread = 3; // have three lines of a 'good' record
  } else if(fRecord[2][0] == '@' && isAlphabet(fRecord[3][0])) {
    sdelta = 2; sread = 2; // have two lines of a 'good' record
  } else if (fRecord[1][0] == '+' && fRecord[3][0] == '@'){
    sdelta = 3; sread = 1; // have only one line of a 'good' record
  } else if(fRecord[2][0] == '+' && isAlphabet(fRecord[1][0])){
    sdelta = 4; sread = 0; // have zero lines of a 'good' record
  } else {
    return false; // bad record!
  }
  // bubble up by swapping
  assert(sdelta >= 0);
  for(int i = 0;(i + sdelta) < 4; i++)
    std::swap(fRecord[i], fRecord[i + sdelta]);
  // read extra lines to fill up the record
  for(int i = sread; i < 4; i++) {
    fRecord[i] = "";
    if(!fqfs->good())
      return false;
    std::getline(*fqfs, fRecord[i]);
    fRecord[i] = trim(fRecord[i]);
    if(fRecord[i].length() == 0)
      return false;
  }
  return true;
}

void updateStrStore(const std::string& in_str,
                    cvec_t &StrStore,ivec_t &Offset,
                    unsigned long& position){
    StrStore.resize(StrStore.size() + in_str.length() + 1);
    memcpy(&StrStore[position], in_str.c_str(), in_str.length() + 1);
    Offset.push_back(position);
    position += in_str.length() + 1;
}

bool readBatch(std::ifstream* fqfs, const long& batchSize,
               const unsigned long& offsetEnd, ReadStore& rbatch)
{
    unsigned long rpos = 0, qpos = 0;
    bool lastRead = false;
    std::string fqRecord[4] = {std::string(""), std::string(""),
                               std::string(""), std::string("")};
    if(batchSize == 0)
        return true;

    // 1. read first record : for handling special case
    if(!readFirstFastqRecord(fqfs, fqRecord))
        return false;

    updateStrStore(fqRecord[1], rbatch.readsString, rbatch.readsOffset, rpos);
    updateStrStore(fqRecord[3], rbatch.qualsString, rbatch.qualsOffset, qpos);
    rbatch.readId = std::stoi(fqRecord[0].substr(1).c_str());

    // 2. read as much as batch size
    for(long j=1; j < batchSize ; j++){
        for(int i = 0; i < 4; i++)
            fqRecord[i] = "";

        if(!readFastqRecord(fqfs, fqRecord)){
            lastRead = true;
            break;
        }

        updateStrStore(fqRecord[1], rbatch.readsString, rbatch.readsOffset, rpos);
        updateStrStore(fqRecord[3], rbatch.qualsString, rbatch.qualsOffset, qpos);

        if(!fqfs->good() || fqfs->tellg() >= (std::streamoff) offsetEnd){
            lastRead = true;
            break;
        }
    }

    return lastRead;
}


void print_time (const std::string& msg, double& timing){
    double cur_time = get_time();
    std::cout << msg << "(" << cur_time - timing << " secs)\n\n";
    timing = cur_time;
}

double get_time() {
      timeval t;
      gettimeofday(&t, 0);
      return t.tv_sec + (0.000001 * t.tv_usec);
} // get_time


timespec local_time(){
  struct timespec tstart;
#ifdef __MACH__
  struct timeval now;
  int rv = gettimeofday(&now, NULL);
  if (!rv) {
      tstart.tv_sec  = now.tv_sec;
      tstart.tv_nsec = now.tv_usec * 1000;
  }
#else
  clock_gettime(CLOCK_REALTIME, &tstart);
#endif
  return tstart;
}


double elapsed(const clock_t& end, const clock_t& start){
  return (double (end - start))/ ((double) CLOCKS_PER_SEC);
}

double elapsed_local(const timespec& finish, const timespec& start){
  double tdiff;
  tdiff = (finish.tv_sec - start.tv_sec);
  tdiff += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  return tdiff;
}

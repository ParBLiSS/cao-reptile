#include "util.h"

void Para::setPara(const char *configFile) {
    cacheOptimizedSearch = 0;
    absentKmers = false;
    std::string line, s1;
    std::ifstream input(configFile);
    std::istringstream buf;
    storeReads = 0;
    kmerCacheSize = tileCacheSize = 4;
    writeOutput = 1;
    writeSpectrum = -1;
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
            else if (s1 == "WriteOutput")
                buf >> writeOutput;
            else if (s1 == "KmerSpectrumOutFile")
                buf >> writeSpectrum;
            else if (s1 == "KmerSpectrumOutFile")
                buf >> kmerSpectrumOutFile;
            else if (s1 == "TileSpectrumOutFile")
                buf >> tileSpectrumOutFile;
        }
    }
}

bool Para::validate() {
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

void Para::compute_offsets(){
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
#ifdef DEBUG
    long long qfileSize = qfin.tellg();
#endif

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
#ifdef DEBUG
    long quendTillLineNo = strtoul (qendLineNo.c_str(),NULL,0) - 1;
#endif

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
void Para::compute_readcount(){
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

void Para::compute_parallelio_params(int fileRecordLength,
                       long& perprocessor,long& startfrom,
                       long &offset, long& inLastBatch) {
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

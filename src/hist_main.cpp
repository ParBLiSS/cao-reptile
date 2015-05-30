#include <mpi.h>
#include <sstream>
#include "ECRunStats.hpp"
#include "ECData.hpp"

extern "C" void hist_reduce(void* in, void* inout, int* len, MPI_Datatype*) {
    int n = *len;
    long* srcHist = static_cast<long*>(in);
    long* destHist = static_cast<long*>(inout);
    for(int i = 0; i < n; i++){
        destHist[i] += srcHist[i];
    }
}

void load_reads(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart_read_p = local_time();
    // If we have to store the reads, we read and store the reads
    if(ecdata.getParams().storeReads) {
        unsigned nReads;
        ecdata.getReadsFromFile(nReads);
    }

    ecstx.tstop_read_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateFileReadTime(ofs);
    }
}

void hist_dist_spectrum(ECData& ecdata, ECRunStats& ecstx,
                        std::ostream& ofs, bool qFlag){
    ecstx.tstart = MPI_Wtime();  ecstx.tstart_kmer_p = local_time();

    local_tile_spectrum(ecdata, qFlag);
    dist_tile_spectrum(ecdata);

    ecstx.tstop_kmer_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();

    ecstx.updateDistSpectrumTime(ecdata, ofs);
}

void make_hist(ECData& ecdata, std::ostream& ofs,
               std::ostream& histfs){
    const int histLength = UCHAR_MAX + 1;
    double tstart = MPI_Wtime(), tstop;

    long kmerHist[histLength], localHist[histLength];
    for(int i = 0; i < histLength; i++)
        localHist[i] = kmerHist[i] = 0;

    // compute histogram
    ecdata.makeTileHistogram(localHist);

    // add-up histogram values
    MPI_Op opSum;
    MPI_Op_create(hist_reduce, 1, &opSum);
    MPI_Reduce(&localHist[0], &kmerHist[0], histLength, MPI_LONG,
               opSum, 0, MPI_COMM_WORLD);

    MPI_Op_free(&opSum);

    // write histogram to output file
    if(ecdata.getParams().m_rank == 0){
        for(int i = 0; i < histLength; i ++)
            histfs << i << "\t" << kmerHist[i] << std::endl;
        histfs << std::endl;
        tstop = MPI_Wtime();
        ofs << "hist construction\t" << tstop - tstart << std::endl;
    }
}

void hist_spectrum(Para& params){
    std::ostream& ofs = std::cout;
    ECData ecd(params);
    ECRunStats ecstx;
    std::string oHistName = params.oErrName + ".hist";
    std::ofstream oHandle;
    if(params.m_rank == 0)
        oHandle.open(oHistName.c_str());

    load_reads(ecd, ecstx, ofs);

    hist_dist_spectrum(ecd, ecstx, ofs, false);

    if(params.m_rank == 0)
        oHandle << "Tile Spectrum Histogram" << std::endl;
    make_hist(ecd, ofs, oHandle);

    // clear current array
    ecd.replaceKArray(0, 0, 0);
    ecd.replaceTileArray(0, 0, 0);

    hist_dist_spectrum(ecd, ecstx, ofs, true);

    if(params.m_rank == 0)
        oHandle << "Tile Qual Spectrum Histogram" << std::endl;
    make_hist(ecd, ofs, oHandle);
}

bool validate_params(Para& params){
    params.tCard = 0;
    if (params.iFaName.empty()) {
        if(params.m_rank == 0)
            std::cout << "Err: InFaFile is not specified!\n";
        return false;
    }
    if (params.oErrName.empty()) {
        if(params.m_rank == 0)
            std::cout << "Err: Output file is not specified!\n";
        return false;
    }
    if(params.m_rank == 0){
        std::string oHistName = params.oErrName + ".hist";
        std::ofstream oHandle(oHistName.c_str());
        if (!oHandle.good()) {
            std::cout << "err output : open " << oHistName
                      << " failed, correct path?\n";
            return false;
        }
        oHandle.close();
    }
    if ((params.K < 0) || (params.K > 16) || (params.step < 1) ||
        (params.tileLength < params.K) || (params.tileLength > 32)) {
        if(params.m_rank == 0)
            std::cout <<
                "Set K in the range of (0, 16] and K+step in the range of (K, 32]\n";
        return false;
    }
    std::ifstream read_stream(params.iFaName.c_str());
    if(!read_stream.good()) {
        std::cout << "open " << params.iFaName << "failed :|\n";
        return false;
    }
    if(params.m_rank == 0) {
        std::stringstream oss;
        oss << "--" << std::endl
            << "parameter" << "\t" << "value" << std::endl;
        oss << "short reads file\t" << params.iFaName << "\n";
        oss << "O/ErrFile\t" << params.oErrName << "\n";
        oss << "(K, step, tile)\t"
            << "(" << params.K << "," << params.step << ","
            << params.tileLength << ")\n"
            << "BatchSize\t" << params.batchSize << "\n"
            << "QThreshold\t" << params.qualThreshold << "\n"
            << "MaxBadQPerKmer\t" << params.maxBadQPerKmer << "\n"
            << "T_card\t" << params.tCard << "\n"
            << "StoreReads\t" << params.storeReads << "\n"
            << "--\n";
        std::cout << oss.str();
        std::cout.flush();
    }
    return true;
}

int main(int argc,char *argv[]){

    try{
        MPI::Init(argc, argv);
    } catch(MPI::Exception& ex) {
        if(ex.Get_error_code() != MPI_SUCCESS) {
            std::cout << "Error starting MPI program. Terminating.\n" ;
            MPI::COMM_WORLD.Abort(ex.Get_error_code());
        }
    }
    // Basic Validations
    if (argc <= 1) {
        if(MPI::COMM_WORLD.Get_rank()  == 0) {
            std::cout << "Syntax:preptile /path/to/config-file" << std::endl;
        }
        MPI::COMM_WORLD.Abort(1);
    }

    std::ifstream input(argv[1]);
    if(!input) {
        if(MPI::COMM_WORLD.Get_rank() == 0){
            std::cout << "ERROR:Can not open Input File!" << std::endl;
            std::cout << "Syntax:preptile /path/to/config-file" << std::endl;
        }
        MPI::COMM_WORLD.Abort(1);
    }
    input.close();
    Para params(argv[1]);
    // Validations on the parameters given in config file
    if(validate_params(params) == false) {
        if(params.m_rank == 0)
            std::cout << "Validation Failed" << std::endl;
        MPI::COMM_WORLD.Abort(1);
    }

    hist_spectrum(params);

    MPI::Finalize();
    return 0;
}

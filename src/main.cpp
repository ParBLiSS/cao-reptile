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

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "ECData.hpp"
#include "ECDriver.hpp"
#include "find_neighbors.h"
#include "count_kmers.hpp"
#include "sort_kmers.hpp"
#include <time.h>

double elapsed(clock_t& end, clock_t& start){
  return (double (end - start))/ ((double) CLOCKS_PER_SEC);
}

struct ECStats{
    double tstartInit, tstart, tstop;
    double read_sync_start, read_sync_stop,
        kmer_sync_start, kmer_sync_stop,
        ec_sync_start, ec_sync_stop;
    clock_t tstart_read_p, tstop_read_p,
        tstart_kmer_p, tstop_kmer_p,
        tstart_ec_p, tstop_ec_p;
    // update global timings
    void updateFileReadTime(std::ostream& ofs);
    void updateSpectrumTime(ECData& ecd, std::ostream& ofs);
    void updateECTime(std::ostream& ofs);
    // report timigs
    void reportTimings(Para& params, std::ostream& ofs);
    void reportQueryCounts(ECData& ecd, std::ostream& ofs);
};

void run_reptile(ECData& ecdata,Para& params){
    std::stringstream out;
    out << params.oErrName << params.mpi_env->rank() ;
    std::string filename = out.str();
    if(params.writeOutput != 0){
        std::ofstream oHandle(filename.c_str());
        if (!oHandle.good()) {
            std::cout << "open " << filename << " failed, correct path?\n";
            return;
        }
        oHandle.close();
    }

    // Run reptile
    ECDriver ecdr(ecdata, filename, params);
    // Commented since it is no longer used
    // if(params.useMaskedLists) {
    //     ecdr.tableMaker(*params);

    //     tstop = MPI_Wtime();
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (params.mpi_env->rank() == 0) {
    //         std::cout << "TIME TO BUILD TABLE " << tstop-tstart
    //                   << " (secs)" << std::endl;
    //     }
    //     tstart = tstop;
    // }
    ecdr.ec();
    return;
}

int parallelEC( char *inputFile){
    Para params(inputFile);
    empi::MPI_env *mpi_env = params.mpi_env;
    std::ostream& ofs = std::cout;
    // Validations on the parameters given in config file
    if(params.validate() == false) {
        if(mpi_env->rank() == 0)
            std::cout << "Validation Failed" << std::endl;
        MPI::COMM_WORLD.Abort(1);
    }

    // Object to encapsulate error-correction data
    ECData ecdata(params);
    ECStats ecstx;
    ecstx.tstart = ecstx.tstartInit = MPI_Wtime(); ecstx.tstart_read_p = clock();
    // If we have to store the reads, we read and store the reads
    if(params.storeReads) {
        getReadsFromFile(ecdata);
    }
    ecstx.tstop_read_p = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        ecstx.updateFileReadTime(ofs);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    ecstx.tstart = MPI_Wtime();  ecstx.tstart_kmer_p = clock();
    // counts the k-mers and loads them in the ECData object
    count_kmers(ecdata);
    // sort kmers and tiles
    sort_kmers(ecdata);
    ecstx.tstop_kmer_p = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        ecstx.updateSpectrumTime(ecdata, ofs);
    }

    ecstx.tstart = MPI_Wtime(); ecstx.tstart_ec_p = clock();
    // Cache Optimized layout construction
    ecdata.buildCacheOptimizedLayout();
    // run reptile
    run_reptile(ecdata, params);
    ecstx.tstop_ec_p = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        ecstx.updateECTime(ofs);
    }

    ecstx.reportTimings(params, ofs);
#ifdef QUERY_COUNTS
    ecstx.reportQueryCounts(ecdata, ofs);
#endif
    ecdata.writeSpectrum();
    return 0;
}

void ECStats::reportTimings(Para& params, std::ostream& ofs){
    // Output for counting the number of failures and success
    //std::stringstream out;
    //out << params.oErrName << params.mpi_env->rank() ;
    //ecdata.output(out.str());
    empi::MPI_env *mpi_env = params.mpi_env;
    int p = params.mpi_env->size();
    if(mpi_env->rank() == 0){
        std::stringstream oss;
        oss << "--" << std::endl
            << "nproc" << "\t" << "phase" << "\t"
            << "start" << "\t" << "stop" << "\t"
            << "duration" << std::endl;
        oss << p << "\t" << "read global" << "\t"
            << read_sync_start << "\t"
            << read_sync_stop << "\t"
            << read_sync_stop - read_sync_start
            << std::endl;
        oss << p << "\t" << "kmer global" << "\t"
            << kmer_sync_start << "\t"
            << kmer_sync_stop << "\t"
            << kmer_sync_stop - kmer_sync_start
            << std::endl;
        oss << p << "\t" << "ec global" << "\t"
            << ec_sync_start << "\t"
            << ec_sync_stop << "\t"
            << ec_sync_stop - ec_sync_start
            << std::endl;
        oss << p << "\t" << "final global" << "\t"
            << tstartInit << "\t" << tstop << "\t"
            << (tstop - tstartInit) << std::endl;
        oss << "--" << std::endl;
        oss << "rank" << "\t" << "phase" << "\t"
            << "start" << "\t" << "stop" << "\t"
            << "duration"
            << std::endl;
        ofs << oss.str();
        ofs.flush();
    }
    for(int i = 0; i < p; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        if(i == mpi_env->rank()){
            std::stringstream oss;
            oss << i << "\t" << "read local" << "\t"
                << elapsed(tstart_read_p, tstart_read_p) << "\t"
                << elapsed(tstop_read_p, tstart_read_p) << "\t"
                << elapsed(tstop_read_p, tstart_read_p)
                << std::endl;
            oss << i << "\t" << "kmer local" << "\t"
                << elapsed(tstart_kmer_p, tstart_read_p) << "\t"
                << elapsed(tstop_kmer_p, tstart_read_p) << "\t"
                << elapsed(tstop_kmer_p, tstart_kmer_p)
                << std::endl;
            oss << i << "\t" << "ec local" << "\t"
                << elapsed(tstart_ec_p, tstart_read_p) << "\t"
                << elapsed(tstop_ec_p, tstart_read_p) << "\t"
                << elapsed(tstop_ec_p, tstart_ec_p)
                << std::endl;
            ofs << oss.str();
            ofs.flush();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void ECStats::reportQueryCounts(ECData& ecdata, std::ostream& ofs){
    empi::MPI_env *mpi_env = ecdata.getParams().mpi_env;
    if(mpi_env->rank() == 0) {
      std::stringstream oss;
      oss << "--" << std::endl;
      oss << "proc" << "\t" << "type" << "\t" << "query counts"  << "\t"
          << "query fails" << "\t" << "query success" << std::endl;
      ofs << oss.str();
      ofs.flush();
    }
    int p = mpi_env->size();
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < p; i++){
      MPI_Barrier(MPI_COMM_WORLD);
      if(i == mpi_env->rank()){
        std::stringstream oss;
        oss << i << "\t" << "kmer" << "\t" << ecdata.getKmerQueries()
            << "\t" << ecdata.getKmerQueryFails() << "\t"
            << (ecdata.getKmerQueries()) - (ecdata.getKmerQueryFails())
            << std::endl;
        oss << i << "\t" << "tile" <<  "\t" << ecdata.getTileQueries()
            << "\t" << ecdata.getTileQueryFails() << "\t"
            << (ecdata.getTileQueries()) - (ecdata.getTileQueryFails())
            << std::endl;
        ofs << oss.str();
        ofs.flush();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_env->rank() == 0) {
        std::stringstream oss;
        oss << "proc" << "\t" << "type" << "\t" ;
        for(unsigned j = 0; j < MAX_LEVELS; j++)
            oss << "L" << j << "\t" ;
        oss <<  std::endl;
        ofs << oss.str();
        ofs.flush();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < p; i++){
      MPI_Barrier(MPI_COMM_WORLD);
      if(i == mpi_env->rank()){
        std::stringstream oss;
        oss << i << "\t" << "kmer";
        for(unsigned j = 0; j < MAX_LEVELS; j++)
          oss << "\t" << ecdata.getKmerLevels()[j];
        oss << std::endl;
        oss << i << "\t" << "tile";
        for(unsigned j = 0; j < MAX_LEVELS; j++)
          oss << "\t" << ecdata.getTileLevels()[j];
        oss << std::endl;
        ofs << oss.str();
        ofs.flush();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void ECStats::updateFileReadTime(std::ostream&){
    read_sync_start = tstart;
    read_sync_stop = tstop;
}

void ECStats::updateSpectrumTime(ECData& ecdata, std::ostream& ofs){
    std::stringstream oss;
    oss << "kmer count\t" << ecdata.getKmerCount() << std::endl;
    oss << "tile count\t" << ecdata.getTileCount() << std::endl;
    std::cout << oss.str();
    std::cout.flush();
    kmer_sync_start = tstart;
    kmer_sync_stop = tstop;
    oss << "K-SPECTRUM CONSTRUCTION TIME " << tstop-tstart
        << " (secs)" << std::endl;
    ofs << oss.str();
    ofs.flush();
}

void ECStats::updateECTime(std::ostream& ofs){
    std::stringstream oss;
    ec_sync_start = tstart;
    ec_sync_stop = tstop;
    oss << "ERR CORRECTION TIME " << tstop-tstart
        << " (secs)" << std::endl;
    oss << "TOTAL TIME " << tstop-tstartInit
        << " (secs)" << std::endl;
    ofs << oss.str();
    ofs.flush();
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
    /*
    if(MPI::COMM_WORLD.Get_size() < 2) {
        std::cout << "ERROR:Requires at least two process" << std::endl;
        MPI::COMM_WORLD.Abort(1);
    }
    */
    std::ifstream input(argv[1]);
    if(!input) {
        if(MPI::COMM_WORLD.Get_rank() == 0){
            std::cout << "ERROR:Can not open Input File!" << std::endl;
            std::cout << "Syntax:preptile /path/to/config-file" << std::endl;
        }
        MPI::COMM_WORLD.Abort(1);
    }
    input.close();

    parallelEC(argv[1]);

    MPI::Finalize();
    return 0;
}

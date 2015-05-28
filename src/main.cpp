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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <mpi.h>

#include "ECData.hpp"
#include "ECDriver.hpp"
#include "count_kmers.hpp"
#include "sort_kmers.hpp"
#include "ECRunStats.hpp"

timespec local_time(){
  struct timespec tstart;
  clock_gettime(CLOCK_REALTIME, &tstart);
  return tstart;
}

void construct_dist_spectrum(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart = MPI_Wtime();  ecstx.tstart_kmer_p = local_time();

    // counts the k-mers and loads them in the ECData object
    local_kmer_spectrum(ecdata);
    // sort kmers
    dist_kmer_spectrum(ecdata);

    local_tile_spectrum(ecdata);
    dist_tile_spectrum(ecdata);

    ecstx.tstop_kmer_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();

    ecstx.updateDistSpectrumTime(ecdata, ofs);

}

void load_spectrum(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart = MPI_Wtime();  ecstx.tstart_kmer_p = local_time();

    ecdata.loadSpectrum();

    ecstx.tstop_kmer_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateSpectrumTime(ecdata, ofs);
    }
}

void load_reads(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart_read_p = local_time();
    // If we have to store the reads, we read and store the reads
    if(ecdata.getParams().storeReads) {
        ecdata.getReadsFromFile();
    }

    ecstx.tstop_read_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateFileReadTime(ofs);
    }
}

void construct_spectrum(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart = MPI_Wtime();  ecstx.tstart_kmer_p = local_time();

    // counts the k-mers and loads them in the ECData object
    count_kmers(ecdata);
    // sort kmers and tiles
    sort_kmers(ecdata);
    // gather kemers and tiles
    gather_spectrum(ecdata);

    ecstx.tstop_kmer_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateSpectrumTime(ecdata, ofs);
    }
}

void run_reptile(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart = MPI_Wtime(); ecstx.tstart_ec_p = local_time();

    // Cache Optimized layout construction
    ecdata.buildCacheOptimizedLayout();

    // Run reptile
    ECDriver ecdr(ecdata, ecdata.getParams());
    ecdr.ec();

    ecstx.tstop_ec_p = local_time();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateECTime(ofs);
    }
}

void dist_spectrum(Para& params){
    std::ostream& ofs = std::cout;
    ECData ecdata(params);
    ECRunStats ecstx;

    load_reads(ecdata, ecstx, ofs);

    construct_dist_spectrum(ecdata, ecstx, ofs);

    ecdata.writeDistSpectrum();

    ecstx.reportTimings(params, ofs);
}

void run_ec_only(Para& params){
    std::ostream& ofs = std::cout;
    ECData ecdata(params);
    ECRunStats ecstx;

    load_reads(ecdata, ecstx, ofs);

    load_spectrum(ecdata, ecstx, ofs);

    run_reptile(ecdata, ecstx, ofs);

    ecstx.reportTimings(params, ofs);

#ifdef QUERY_COUNTS
    ecstx.reportQueryCounts(ecdata, ofs);
#endif

    ecdata.writeSpectrum();
}

int parallelEC(Para& params){
    std::ostream& ofs = std::cout;

    ECData ecdata(params);
    ECRunStats ecstx;

    load_reads(ecdata, ecstx, ofs);

    construct_spectrum(ecdata, ecstx, ofs);

    run_reptile(ecdata, ecstx, ofs);

    ecstx.reportTimings(params, ofs);

#ifdef QUERY_COUNTS
    ecstx.reportQueryCounts(ecdata, ofs);
#endif

    ecdata.writeSpectrum();
    return 0;
}

int main(int argc,char *argv[]){
    int provided;
    try{
       //MPI::Init(argc, argv);

       provided = MPI::Init_thread(argc, argv, MPI_THREAD_FUNNELED);
       assert(provided == MPI_THREAD_FUNNELED);
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
    if(params.validate() == false) {
        if(params.m_rank == 0)
            std::cout << "Validation Failed" << std::endl;
        MPI::COMM_WORLD.Abort(1);
    }

    if(params.runType == 0)
        parallelEC(params);
    else if(params.runType == 1)
        dist_spectrum(params);
    else if(params.runType == 2)
        run_ec_only(params);

    MPI::Finalize();
    return 0;
}

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
#include "kmer_count.hpp"
#include "find_neighbors.h"
#include "sort_kmers.hpp"
#include "ECData.hpp"
#include "Parser.h"
#include <time.h>

void run_reptile(ECData *ecdata,Para *params){

    // Run reptile
    Parser myParser(ecdata);
    // Commented since it is no longer used
    // if(params->useMaskedLists) {
    //     myParser.tableMaker(*params);

    //     tstop = MPI_Wtime();
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (params->mpi_env->rank() == 0) {
    //         std::cout << "TIME TO BUILD TABLE " << tstop-tstart
    //                   << " (secs)" << std::endl;
    //     }
    //     tstart = tstop;
    // }

    myParser.ec(*params);
    std::stringstream out;
    out << params->oErrName << params->mpi_env->rank() ;
    if(params->writeOutput != 0)
      myParser.output(out.str());
    return;
}

int parallelEC( char *inputFile){
    Para *params = new Para(inputFile);
    empi::MPI_env *mpi_env = params->mpi_env;

    // Validations on the parameters given in config file
    if(params->validate() == false) {
        if(mpi_env->rank() == 0)
            std::cout << "Validation Failed" << std::endl;
        MPI::COMM_WORLD.Abort(1);
    }

    // Object to encapsulate error-correction data
    ECData *ecdata = new ECData(params);
    double tstartInit = MPI_Wtime(),
        tstart = tstartInit;
    double read_sync_start, read_sync_stop,
        kmer_sync_start, kmer_sync_stop,
        ec_sync_start, ec_sync_stop;
    time_t tstart_read_p, tstop_read_p,
        tstart_kmer_p, tstop_kmer_p,
        tstart_ec_p, tstop_ec_p;
    // If we have to store the reads, we read and store the reads
    time(&tstart_read_p);
    if(params->storeReads) {
        getReadsFromFile(ecdata);
    }
    time(&tstop_read_p);
    MPI_Barrier(MPI_COMM_WORLD);
    double tstop = MPI_Wtime();

    if (mpi_env->rank() == 0) {
        read_sync_start = tstart;
        read_sync_stop = tstop;
        std::cout << "READING FILE " << tstop-tstart
                  << " (secs)" << std::endl;
    }

     MPI_Barrier(MPI_COMM_WORLD);

    tstart = MPI_Wtime();
    time(&tstart_kmer_p);
    // counts the k-mers and loads them in the ECData object
    kmer_count(ecdata);
    // sort kmers and tiles
    kmer_sort(ecdata);
    time(&tstop_kmer_p);

    MPI_Barrier(MPI_COMM_WORLD);
    tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        kmer_sync_start = tstart;
        kmer_sync_stop = tstop;
        std::cout << "K-Spectrum COnstruction Time " << tstop-tstart
                  << " (secs)" << std::endl;
    }
    tstart = MPI_Wtime();
    time(&tstart_ec_p);
    // Cache Optimized layout construction
    ecdata->buildCacheOptimizedLayout();
    // run reptile
    run_reptile(ecdata, params);
    time(&tstop_ec_p);
    MPI_Barrier(MPI_COMM_WORLD);
    tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        ec_sync_start = tstart;
        ec_sync_stop = tstop;
        std::cout << "ERR CORRECTION TIME " << tstop-tstart
                  << " (secs)" << std::endl;
        std::cout << "TOTAL TIME " << tstop-tstartInit
                  << " (secs)" << std::endl;
    }
    // Output for counting the number of failures and success
    //std::stringstream out;
    //out << params->oErrName << params->mpi_env->rank() ;
    //ecdata->output(out.str());
    int p = mpi_env->size();
    if(mpi_env->rank() == 0){
        std::cout << "rank" << "\t" << "phase" << "\t"
                  << "start" << "\t" << "stop" << "\t"
                  << "duration"
                  << std::endl;
        std::cout << "0" << "\t" << "read global" << "\t"
                  << read_sync_start << "\t"
                  << read_sync_stop << "\t"
                  << read_sync_stop - read_sync_start
                  << std::endl;
        std::cout << "0" << "\t" << "kmer global" << "\t"
                  << kmer_sync_start << "\t"
                  << kmer_sync_stop << "\t"
                  << kmer_sync_stop - kmer_sync_start
                  << std::endl;
        std::cout << "0" << "\t" << "ec global" << "\t"
                  << ec_sync_start << "\t"
                  << ec_sync_stop << "\t"
                  << ec_sync_stop - ec_sync_start
                  << std::endl;
    }
    for(int i = 0; i < p; i++){
        if(i == mpi_env->rank()){
            std::cout << i << "\t" << "read local" << "\t"
                      << difftime(tstart_read_p, tstart_read_p) << "\t"
                      << difftime(tstop_read_p, tstart_read_p) << "\t"
                      << difftime(tstop_read_p, tstart_read_p)
                      << std::endl;
            std::cout << i << "\t" << "kmer local" << "\t"
                      << difftime(tstart_kmer_p, tstart_read_p) << "\t"
                      << difftime(tstop_kmer_p, tstart_read_p) << "\t"
                      << difftime(tstop_kmer_p, tstart_kmer_p)
                      << std::endl;
            std::cout << i << "\t" << "ec local" << "\t"
                      << difftime(tstart_kmer_p, tstart_read_p) << "\t"
                      << difftime(tstop_kmer_p, tstart_read_p) << "\t"
                      << difftime(tstop_kmer_p, tstart_kmer_p)
                      << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef QUERY_COUNTS
    std::cout << "proc" << "\t" << "type" << "\t" << "query counts"  << "\t"
              << "query fails" << "\t" << "query success" << std::endl;
    for(int i = 0; i < p; i++){
        if(i == mpi_env->rank()){
            std::cout << i << "\t" << "kmer" <<  ecdata->m_kmerQueries << "\t"
                      << ecdata->m_kmerQueryFails << "\t"
                      << (ecdata->m_kmerQueries) - (ecdata->m_kmerQueryFails)
                      << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    delete params;
    delete ecdata;
    return 0;
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

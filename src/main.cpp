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

void load_reads(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart_read_p = clock();
    // If we have to store the reads, we read and store the reads
    if(ecdata.getParams().storeReads) {
        ecdata.getReadsFromFile();
    }

    ecstx.tstop_read_p = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateFileReadTime(ofs);
    }
}

void construct_spectrum(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart = MPI_Wtime();  ecstx.tstart_kmer_p = clock();

    // counts the k-mers and loads them in the ECData object
    count_kmers(ecdata);
    // sort kmers and tiles
    sort_kmers(ecdata);

    ecstx.tstop_kmer_p = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateSpectrumTime(ecdata, ofs);
    }
}

void run_reptile(ECData& ecdata, ECRunStats& ecstx, std::ostream& ofs){
    ecstx.tstart = MPI_Wtime(); ecstx.tstart_ec_p = clock();

    // Cache Optimized layout construction
    ecdata.buildCacheOptimizedLayout();

    // Run reptile
    ECDriver ecdr(ecdata, ecdata.getParams().outputFilename,
                  ecdata.getParams());
    ecdr.ec();

    ecstx.tstop_ec_p = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    ecstx.tstop = MPI_Wtime();
    if (ecdata.getParams().m_rank == 0) {
        ecstx.updateECTime(ofs);
    }
}

int parallelEC( char *inputFile){
    Para params(inputFile);

    std::ostream& ofs = std::cout;
    // Validations on the parameters given in config file
    if(params.validate() == false) {
        if(params.m_rank == 0)
            std::cout << "Validation Failed" << std::endl;
        MPI::COMM_WORLD.Abort(1);
    }

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

    parallelEC(argv[1]);

    MPI::Finalize();
    return 0;
}

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

void run_reptile(ECData *ecdata,Para *params, double &tstart, double &tstop){

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
    // If we have to store the reads, we read and store the reads 
    if(params->storeReads) {
        getReadsFromFile(ecdata);
    }
     MPI_Barrier(MPI_COMM_WORLD);
    double tstop = MPI_Wtime();  
    
    if (mpi_env->rank() == 0) {
        std::cout << "READING FILE " << tstop-tstart 
                  << " (secs)" << std::endl;
    }

     MPI_Barrier(MPI_COMM_WORLD);

    tstart = MPI_Wtime();
    // counts the k-mers and loads them in the ECData object
    kmer_count(ecdata); 
    // sort kmers and tiles
    kmer_sort(ecdata);

    MPI_Barrier(MPI_COMM_WORLD);
    tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        std::cout << "K-Spectrum COnstruction Time " << tstop-tstart 
                  << " (secs)" << std::endl;
    }
    tstart = tstop;

    // run reptile
    run_reptile(ecdata,params,tstart,tstop);

    MPI_Barrier(MPI_COMM_WORLD);
    tstop = MPI_Wtime();
    if (mpi_env->rank() == 0) {
        std::cout << "ERR CORRECTION TIME " << tstop-tstart 
                  << " (secs)" << std::endl;
        std::cout << "TOTAL TIME " << tstop-tstartInit
                  << " (secs)" << std::endl;
    }

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




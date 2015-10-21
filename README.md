# PSbEC : Parallel Spectrum-based Error Correction

PSbEC is a parallel spectrum-based error correction software. It is
developed as a generic framework to parallelize any spectrum-based error
correction and includes three major optimizations (1) dynamic load
balancing , (2) cache optimized array layout, and (3) Use of shared
memory threads to facilitate processing of large spectrum.

Though any spectrum-based error correction can be parallelized with our
framework, currently we have demonstrated it by a parallel
implementation of the error correction software
[Reptile](http://aluru-sun.ece.iastate.edu/doku.php?id=reptile).

# Installation
## Dependencies

* A modern, C++11 ready compiler such as `g++` version 4.8 or higher or `clang` version 3.2 or higher.
* C++11 threads.
* The [cmake](www.cmake.org) build system (Version >= 2.8.11).
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* MPI installation. Tested with MPICH2 and Intel MPI.

## Compilation

Create a build directory. For example,

     mkdir build
     cd build

Configure and build the executable 'cao-preptile'. Continuing the example,

     cmake ..
     make

If successfully built, cao-reptile  should be available in the build directory.

# Usage

Input config file should be prepared similar to Reptile. Documentation for
setting parameters for reptile in the config file is given at the
Reptile page [here](http://aluru-sun.ece.iastate.edu/doku.php?id=reptile) .

We follow Reptile's conventions except for the read input format. Unlike Reptile, we require the read input to be in the fastq format. We also provide the script fastq-filter-v1.0.pl to filter reads with 'N' character.

We provide the following additional options related to the parallel algorithm and the optimizations.

1. RunType : Possible options 0, 1, and 2. Set to 0 to construct
   spectrum and run error correction, Set it to 1 to run only
   distributed spectrum construction, Set it to 2 to run only error
   correction.
2. KmerSpectrumInFile : Input kmer spectrum file when RunType is set
   to 2. Input must be generated using the WriteSpectrum option.
3. TileSpectrumInFile : Input tile spectrum file when RunType is set
   to 2. Input must be generated using the WriteSpectrum option.
4. CacheOptimizedSearch : Default 0, no cache optimized search. Set it
   to 1 for using cache aware layout for kmer/tile spectrums; 2 for
   cache oblivious layout of kmer/tile spectrums.
5. WriteSpectrum : Default 0. If set to 1, then the spectrum is written
   to an output file.
6. Threads : Default 1. No. of shared memory threads per node.
7. WorkDistribution : Default 0. Set to 0 for static distribution, 1 for
   dynamic distribution.


# Citations

1. Nagakishore Jammula, Sriram Chockalingam, and Srinivas Aluru. _Parallel Read Error Correction for Big Genomic Datasets._ 22nd IEEE International Conference on High Performance Computing (HiPC 2015) (Accepted).
2. Ankit Shah, Sriram Chockalingam, and Srinivas Aluru. _A Parallel Algorithm for Spectrum-based Short Read Error Correction._ 26th IEEE International Symposium on Parallel and Distributed Processing, IPDPS 2012.


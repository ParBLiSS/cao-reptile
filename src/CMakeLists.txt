#
#
cmake_minimum_required(VERSION 2.8.12)
project(cao-reptile)
#
#
add_compile_options(-std=c++11 -Wall -Wextra -O3 -ffast-math -funroll-loops)
include_directories(../mpiext)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
message(STATUS "Found MPI:")
message(STATUS "    headers: ${MPI_INCLUDE_PATH}")
message(STATUS "    libs:    ${MPI_LIBRARIES}")
message(STATUS "    link flags: ${MPI_LINK_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_VERBOSE_MAKEFILE on)
add_executable (cao-reptile caware_layout.cpp  ECData.cpp  find_neighbors.cpp  kmer_count.cpp  main.cpp  Parser.cpp  sort_kmers.cpp)
target_link_libraries(cao-reptile "${MPI_LIBRARIES}")

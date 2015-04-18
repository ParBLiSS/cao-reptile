#
#
cmake_minimum_required(VERSION 2.8.11)
project(cao-reptile)
#
# set of c++ 11 and add warnings
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra  -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -std=c++11 -DDEBUG_EMIT")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -std=c++11")
include_directories(../mpiext)
#
# add MPI includes and libraries
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
message(STATUS "Found MPI:")
message(STATUS "    headers: ${MPI_INCLUDE_PATH}")
message(STATUS "    libs:    ${MPI_LIBRARIES}")
message(STATUS "    link flags: ${MPI_LINK_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${MPI_LINK_FLAGS}")
#
# verbose for listing the command
# set(CMAKE_VERBOSE_MAKEFILE on)
#
# main executable
add_executable (cao-reptile util.cpp ECData.cpp ECDriver.cpp find_neighbors.cpp  count_kmers.cpp   sort_kmers.cpp caware_layout.cpp  coblivious_layout.cpp ECImpl.cpp main.cpp)
target_link_libraries(cao-reptile "${MPI_LIBRARIES}")
#
# test executable
add_executable (cao-reptile-test util.cpp ECData.cpp  ECDriver.cpp find_neighbors.cpp count_kmers.cpp sort_kmers.cpp caware_layout.cpp coblivious_layout.cpp ECImpl.cpp test_main.cpp)
target_link_libraries(cao-reptile-test "${MPI_LIBRARIES}")

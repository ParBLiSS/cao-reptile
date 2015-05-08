#
#
cmake_minimum_required(VERSION 2.8.11)
set(PROJECT_NAME_STR cao-reptile)
PROJECT(${PROJECT_NAME_STR} C CXX)
#
# verbose for listing the command
set(CMAKE_VERBOSE_MAKEFILE on)
#
# set of c++ 11 and add warnings
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra  -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -std=c++11 -DDEBUG_EMIT")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -std=c++11")
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
# adding threads
#
find_package(Threads REQUIRED)
#
# main executable
add_executable (cao-reptile src/util.cpp src/mpi_util.cpp src/ECData.cpp src/ECDriver.cpp src/find_neighbors.cpp  src/count_kmers.cpp src/sort_kmers.cpp src/ECImpl.cpp src/ECWDist.cpp src/ECRunStats.cpp src/main.cpp)
target_link_libraries(cao-reptile "${MPI_LIBRARIES}" stdc++)
#
# test executable
add_executable (cao-reptile-test src/caware_layout.cpp src/coblivious_layout.cpp src/test_main.cpp)
target_link_libraries(cao-reptile-test "${MPI_LIBRARIES}" stdc++)
#
#
#
add_executable(tsample src/wdist_test.cpp src/util.cpp src/mpi_util.cpp)
target_link_libraries(tsample "${MPI_LIBRARIES}" stdc++)

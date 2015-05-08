#include <mpi.h>
#include <stdexcept>
#include "mpi_util.hpp"

template <typename T>
MPI_Datatype get_mpi_dt()
{
    throw std::runtime_error("Unsupported MPI datatype");
    // default to int
    return MPI_INT;
}

template <>
MPI_Datatype get_mpi_dt<double>(){
    return MPI_DOUBLE;
}
template <>
MPI_Datatype get_mpi_dt<int>(){
    return MPI_INT;
}

template <>
MPI_Datatype get_mpi_dt<unsigned>(){
    return MPI_UNSIGNED;
}

template <>
MPI_Datatype get_mpi_dt<long>(){
    return MPI_UNSIGNED;
}

template <>
MPI_Datatype get_mpi_dt<unsigned long>(){
    return MPI_UNSIGNED_LONG;
}

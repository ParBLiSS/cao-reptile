#ifndef MPI_UTIL_H
#define MPI_UTIL_H

template <typename T>
MPI_Datatype get_mpi_dt();

template <>
MPI_Datatype get_mpi_dt<double>();

template <>
MPI_Datatype get_mpi_dt<int>();

template <>
MPI_Datatype get_mpi_dt<unsigned>();

template <>
MPI_Datatype get_mpi_dt<unsigned long>();

template <>
MPI_Datatype get_mpi_dt<long>();

#endif /* MPI_UTIL_H */

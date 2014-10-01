/***
 *  $Id: mpi_hostnames.hpp 715 2007-02-05 15:56:38Z zola $
 **
 *  File: mpi_hostnames.hpp
 *  Developed: Jan 20, 2005
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2005-2007 Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  For copyright details please see attached LICENSE
 */

#ifndef EMPI_MPI_HOSTNAMES_HPP
#define EMPI_MPI_HOSTNAMES_HPP

#include <mpi.h>
#include <string>


/** extMPI main namespace.
 */
namespace empi {

  /** mpi_hostnames performs collective communication to distribute
   *  names of all processors in the communicator @a comm.
   *  @param root is rank of the processor which initiates communication.
   *  @param rank is rank of calling process.
   *  @param size is size of the communicator.
   *  @param comm is MPI communicator.
   *  @param outR is Output Iterator where ranks of the hosts should be written.
   *  @param outN is Output Iterator where names (related to ranks)
   *  should be written.
   */
  template <typename IterR, typename IterN>
  inline void mpi_hostnames(int root, int rank, int size, MPI_Comm comm, 
			    IterR outR, IterN outN) {
      // get my name
      int len;
      char name[MPI_MAX_PROCESSOR_NAME];
      ::MPI_Get_processor_name(name, &len);

      // size of data: we assume that to pack one int 
      // we need sizeof(long int)^2
      const long int SIZE = MPI_MAX_PROCESSOR_NAME
	  + (sizeof(long int) * sizeof(long int));

      // allocate buffer for packing
      int pack_pos = 0;
      int pack_buf[SIZE];

      // pack rank and name
      ::MPI_Pack(&rank, 1, MPI_INT, pack_buf, SIZE, &pack_pos, comm);
      ::MPI_Pack(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 
		 pack_buf, SIZE, &pack_pos, comm);
      
      // allocate receive storage and send data
      char* all_names = new char[size * SIZE];

      ::MPI_Allgather(pack_buf, SIZE, MPI_PACKED,
		      all_names, SIZE, MPI_PACKED, comm);
      
      int my_rank;

      // unpack data and copy to user space
      for (int pos = 0; pos < size; ++pos) {
	  pack_pos = 0;

	  ::MPI_Unpack(all_names + (pos * SIZE), SIZE, &pack_pos, 
		       &my_rank, 1, MPI_INT, comm);
	  ::MPI_Unpack(all_names + (pos * SIZE), SIZE, &pack_pos, 
		       name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, comm);

	  *outR = my_rank;
	  *outN = name;

	  ++outR;
	  ++outN;
      }
    
      delete[] all_names;
  } // mpi_hostnames

} // namespace empi

#endif // EMPI_MPI_HOSTNAMES_HPP

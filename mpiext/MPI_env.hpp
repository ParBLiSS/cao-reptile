/***
 *  $Id: MPI_env.hpp 715 2007-02-05 15:56:38Z zola $
 **
 *  File: MPI_env.hpp
 *  Developed: Nov 03, 2005
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2005-2007 Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  For copyright details please see attached LICENSE
 */

#ifndef EMPI_MPI_ENV_HPP
#define EMPI_MPI_ENV_HPP

#include <map>
#include <vector>
#include <mpi_hostnames.hpp>


/** extMPI main namespace.
 */
namespace empi {
  
  /** MPI_env stores basic information about process and its communicatior.
   */
  class MPI_env {
  public:
      /** Constructs MPI_env object and sets its current communicator
       *  to @a comm, then calls reset().
       *  @param comm is a valid MPI communicator.
       */
      explicit MPI_env(MPI_Comm comm = MPI_COMM_WORLD) : root_(-1) {
	  comm_ = comm;
	  reset();
      } // MPI_env

      /**
       */
      virtual ~MPI_env() { }


      /** Takes current size of the communicator and rank of the process.
       */
      void reset() {
	  ::MPI_Comm_size(comm_, &size_);
	  ::MPI_Comm_rank(comm_, &rank_);
	  
	  int stat;
	  char buf[MPI_MAX_PROCESSOR_NAME];

	  ::MPI_Get_processor_name(buf, &stat);
	  name_ = buf;
      } // reset


      /** This method is a collective all-to-all operation.
       *  All processes in the communicator exchange information
       *  about rank and hostname.
       *  @param root is a rank of the coordinating processor 
       *  (it can be any valid rank in the communicator).
       */
      void publish_names(int root = 0) {
	  std::vector<int> ranks(size_);
	  std::vector<std::string> names(size_);

	  mpi_hostnames(root, rank_, size_, comm_, 
			ranks.begin(), names.begin());
	  
	  for (unsigned int i = 0; i < size_; ++i) {
	      rank2name_.insert(is_pair(ranks[i], names[i]));
	      name2rank_.insert(si_pair(names[i], ranks[i]));
	  }
      } // publish_names

      
      /** This is collective one-to-all operation.
       *  Process @a who broadcasts to the communicator
       *  information about rank of the root process,
       *  and all other processes set this information locally.
       *  @param who is rank of the process which will broadcast info.
       *  @return rank of the root process sent by @a who.
       */
      int ask4root(int who = 0) {
	  ::MPI_Bcast(&root_, 1, MPI_INT, who, comm_);
	  return root_;
      } // ask4root


      /** @return current MPI communicator.
       */
      MPI_Comm comm() const { return comm_; }

      /** @return size of the current MPI communicator.
       */
      int size() const { return size_; }


      /** @return rank of the process in the current communicator.
       */
      int rank() const { return rank_; }

      /** @return rank of the first process running on the host @a name
       *  or -1 if such host is not available.
       */
      int rank(const std::string& name) const {
	  std::multimap<std::string, int>::const_iterator iter
	      = name2rank_.find(name);
	  return iter != name2rank_.end() ? iter->second : -1;
      } // rank

      /** @return MPI name of the host running given process.
       */
      const std::string& name() const { return name_; }

      /** @return name of the host with given @a rank.
       *  If such information is not available the method
       *  returns empty string.
       */
      std::string name(int rank) const {
	  std::map<int, std::string>::const_iterator iter 
	      = rank2name_.find(rank);
	  return iter != rank2name_.end() ? iter->second : "";
      } // name
      

      /** Set local info about root to @a rank.
       *  This method has only local effect!
       *  To set this information for the whole communicator 
       *  the method must be followed by @a ask4root.
       *  @param rank is rank of the root.
       */
      void set_root(int rank) {
	  root_ = rank;
      } // set_root

      /** Get local information about root.
       *  @return rank of the root process. If -1 is returned
       *  it may mean that rank of the root has been not broadcasted,
       *  i.e. @a ask4root has been not called.
       */
      int get_root() const { return root_; }

      /** @return true if calling process is root, false otherwise.
       */
      bool am_I_root() const { return rank_ == root_; }


  protected:
      MPI_Comm comm_; // MPI communicator
      int size_; // size of the communicator
      int rank_; // rank of the process in the communicator
      std::string name_; // MPI name of the host

      int root_; // rank of the root process

      std::map<int, std::string> rank2name_;
      std::multimap<std::string, int> name2rank_;

      // to avoid problems with Sun Studio
      typedef std::map<int, std::string>::value_type is_pair;
      typedef std::multimap<std::string, int>::value_type si_pair;

  }; // class MPI_env

} // namespace empi

#endif // EMPI_MPI_ENV_HPP

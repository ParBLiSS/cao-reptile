/***
 *  $Id: MPI_obuffer.hpp 715 2007-02-05 15:56:38Z zola $
 **
 *  File: MPI_obuffer.hpp
 *  Developed: Feb 09, 2006
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2005-2007 Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  For copyright details please see attached LICENSE
 */

#ifndef EMPI_MPI_OBUFFER_HPP
#define EMPI_MPI_OBUFFER_HPP

#include <mpi.h>
#include "iterator_range.hpp"


/** extMPI main namespace.
 */
namespace empi {

  /** MPI_obuffer provides a convenient way to use MPI_Unpack routine.
   */
  class MPI_obuffer {
  public:
      /** Constructs MPI_obuffer which refers to @a buf
       *  of size @a s. MPI_obuffer does not agregate @a buf!
       *  @param s is size of the buffer.
       *  @param buf is buffer with packed data.
       *  @param comm is communicator that was used for receiving data.
       */
      explicit MPI_obuffer(int s = 0, const void* buf = 0, 
			   MPI_Comm comm = MPI_COMM_WORLD)
	  : comm_(comm), ob_(0), size_(s), buffer_(buf), pos_(0) { }


      /** Resets buffer.
       */
      void reset() { ob_ = 0; pos_ = 0; }

      /** @return size of the buffer.
       */
      int size() const { return size_; }

      /** @return pointer to the buffer with packed data.
       */
      const void* data() const { return buffer_; }

      /** Sets new source of data.
       */
      void data(int s, const void* buf) {
	  size_ = s;
	  buffer_ = buf;
      } // data
       
      /** @return current position in the buffer.
       */
      int position() const { return pos_; }


      template <typename T>
      friend MPI_obuffer& operator,(MPI_obuffer&, T&);

      template <typename T>
      friend MPI_obuffer& operator,(MPI_obuffer&, T*);

      template <typename T>
      MPI_obuffer& operator>>(T& t) { return operator,(*this, t); }


  private:
      MPI_obuffer(const MPI_obuffer&);
      void operator=(const MPI_obuffer&);

      MPI_Comm comm_;

      void* ob_;

      int size_;
      const void* buffer_;
      int pos_;

  }; // class MPI_out_buffer


  // *** HERE GO SPECIALIZATIONS ***
  
  // iterator_range<T*>
  template <typename T> 
  inline MPI_obuffer& operator,(MPI_obuffer& obuf, iterator_range<T>& iter) {
      T i = iter.begin;
      T n = iter.end;
      for (; i != n; ++i) operator,(obuf, *i);
      return obuf;
  } // operator,


  // T*
  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf, char* obj) {

      if (obuf.ob_ == 0) obuf.ob_ = (void*)obj;
      else {
	  int s = (obj - (char*)obuf.ob_);
	  ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		       obuf.ob_, s, MPI_CHAR, obuf.comm_);
	  obuf.ob_ = 0;
      }
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    signed int* obj) {

      if (obuf.ob_ == 0) obuf.ob_ = (void*)obj;
      else {
	  int s = (obj - (signed int*)obuf.ob_);
	  ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		       obuf.ob_, s, MPI_INT, obuf.comm_);
	  obuf.ob_ = 0;
      }
      return obuf;
  } // operator,



  // T&
  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    char& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_CHAR, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    signed char& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_CHAR, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    signed short int& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_SHORT, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    signed int& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_INT, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    signed long int& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_LONG, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    unsigned char& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_UNSIGNED_CHAR, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    unsigned short int& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_UNSIGNED_SHORT, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    unsigned int& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_UNSIGNED, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    unsigned long int& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_UNSIGNED_LONG, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    float& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_FLOAT, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    double& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_, 
		   (void*)&obj, 1, MPI_DOUBLE, obuf.comm_);
      return obuf;
  } // operator,

  template <> inline MPI_obuffer& operator,(MPI_obuffer& obuf,
					    long double& obj) {
      ::MPI_Unpack((void*)obuf.buffer_, obuf.size_, &obuf.pos_,
		   (void*)&obj, 1, MPI_LONG_DOUBLE, obuf.comm_);
      return obuf;
  } // operator,

} // namespace empi

#endif // EMPI_MPI_OBUFFER_HPP

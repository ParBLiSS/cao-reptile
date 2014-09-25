/***
 *  $Id: MPI_ibuffer.hpp 715 2007-02-05 15:56:38Z zola $
 **
 *  File: MPI_ibuffer.hpp
 *  Developed: Feb 09, 2006
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2005-2007 Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  For copyright details please see attached LICENSE
 */

#ifndef EMPI_MPI_IBUFFER_HPP
#define EMPI_MPI_IBUFFER_HPP

#include <mpi.h>
#include "iterator_range.hpp"


/** extMPI main namespace.
 */
namespace empi {

  /** MPI_ibuffer provides a convenient way to use MPI_Pack routine.
   */
  class MPI_ibuffer {
  public:
      /** Constructs MPI_ibuffer of given size.
       *  @param s is size of the buffer to be created.
       *  @param comm is communicator that will be used 
       *  for sending packed data.
       */
      explicit MPI_ibuffer(int s, MPI_Comm comm = MPI_COMM_WORLD) 
	  : comm_(comm), ib_(0), size_(s), pos_(0), del_(true) {
	  buffer_ = size_ > 0 ? new char[size_] : 0;
      } // MPI_ibuffer
      
      /** Constructs MPI_ibuffer on top of buffer @a buf of size @a s.
       *  @param s is size of the buffer to be created.
       *  @param buf is buffer to use.
       *  @param comm is communicator that will be used 
       *  for sending packed data.
       */
      MPI_ibuffer(int s, char* buf, MPI_Comm comm = MPI_COMM_WORLD)
	  : comm_(comm), ib_(0), size_(s), buffer_(buf), pos_(0), del_(false) {
      } // MPI_ibuffer

      /**
       */
      ~MPI_ibuffer() { if (del_ == true) delete[] buffer_; }
      

      /** Resets buffer.
       */
      void reset() { ib_ = 0; pos_ = 0; }
      
      /** @return size of the buffer.
       */
      int size() const { return size_; }

      /** @return pointer to the buffer with packed data.
       */
      const void* data() const { return buffer_; }
       
      /** @return current position in the buffer.
       */
      int position() const { return pos_; }


      /** Operator utilised to insert single data into buffer.
       */
      template <typename T>
      friend MPI_ibuffer& operator,(MPI_ibuffer&, const T&);

      /** Operator utilised to insert array into buffer.
       */
      template <typename T>
      friend MPI_ibuffer& operator,(MPI_ibuffer&, T*);
      
      /** Operator utilised to insert array into buffer.
       */
      template <typename T>
      friend MPI_ibuffer& operator,(MPI_ibuffer&, const T*);

      /** Operator utilised to insert range into buffer.
       */
      template <typename T>
      friend MPI_ibuffer& operator,(MPI_ibuffer&, const iterator_range<T>&);

      /** Operator utilised to insert first element into buffer.
       *  @param t can be either object or beginning of the range.
       */
      template <typename T>
      MPI_ibuffer& operator<<(const T& t) { return operator,(*this,t); }


  private:
      MPI_ibuffer(const MPI_ibuffer&);
      void operator=(const MPI_ibuffer&);

      MPI_Comm comm_;

      void* ib_;

      int size_;
      char* buffer_;
      int pos_;

      bool del_;

  }; // class MPI_ibuffer


  // *** HERE GO SPECIALIZATIONS ***

  // iterator_range<T*>

  template <typename T> 
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<T>& iter) {
      T i = iter.begin;
      T n = iter.end;
      for (; i != n; ++i) operator,(ibuf, *i);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<char*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<signed char*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<signed short int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_SHORT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<signed int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_INT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<signed long int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_LONG,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<unsigned char*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_UNSIGNED_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<unsigned short int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_UNSIGNED_SHORT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<unsigned int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_UNSIGNED,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<unsigned long int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_UNSIGNED_LONG,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<float*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_FLOAT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<double*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_DOUBLE,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<long double*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack(iter.begin, s, MPI_LONG_DOUBLE,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,



  // iterator_range<const T*>
  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<const char*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const signed char*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const signed short int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_SHORT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const signed int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_INT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
			 const iterator_range<const signed long int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_LONG,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const unsigned char*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_UNSIGNED_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const unsigned short int*>& iter)
  {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_UNSIGNED_SHORT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const unsigned int*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_UNSIGNED,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const unsigned long int*>& iter)
  {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_UNSIGNED_LONG,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<const float*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_FLOAT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <>
  inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
				const iterator_range<const double*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_DOUBLE,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline 
  MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
			 const iterator_range<const long double*>& iter) {

      int s = (iter.end - iter.begin);
      ::MPI_Pack((void*)iter.begin, s, MPI_LONG_DOUBLE,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,



  // T*
  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, char* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (char*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_CHAR,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    signed char* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (signed char*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_CHAR, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    signed short int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = 
	      (obj - (signed short int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_SHORT,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    signed int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (signed int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_INT, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    signed long int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (signed long int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_LONG, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    unsigned char* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned char*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED_CHAR,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    unsigned short int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned short int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED_SHORT, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    unsigned int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    unsigned long int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned long int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED_LONG, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    float* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (float*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_FLOAT, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,
  
  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    double* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (double*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_DOUBLE, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    long double* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (long double*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_LONG_DOUBLE,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,



  // const T*
  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf, 
					    const char* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (char*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_CHAR, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed char* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (signed char*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_CHAR, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed short int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = 
	      (obj - (signed short int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_SHORT,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (signed int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_INT, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed long int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (signed long int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_LONG, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned char* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned char*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED_CHAR,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned short int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned short int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED_SHORT, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned long int* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (unsigned long int*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_UNSIGNED_LONG, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const float* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (float*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_FLOAT, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,
  
  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const double* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (double*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_DOUBLE, 
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const long double* obj) {

      if (ibuf.ib_ == 0) ibuf.ib_ = (void*)obj;
      else {
	  int s = (obj - (long double*)ibuf.ib_);
	  ::MPI_Pack(ibuf.ib_, s, MPI_LONG_DOUBLE,
		     ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
	  ibuf.ib_ = 0;
      }

      return ibuf;
  } // operator,



  // const T&
  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const char& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed char& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed short int& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_SHORT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed int& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_INT, 
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const signed long int& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_LONG,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned char& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_UNSIGNED_CHAR,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned short int& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_UNSIGNED_SHORT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned int& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_UNSIGNED,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const unsigned long int& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_UNSIGNED_LONG,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const float& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_FLOAT,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,
  
  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const double& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_DOUBLE,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

  template <> inline MPI_ibuffer& operator,(MPI_ibuffer& ibuf,
					    const long double& obj) {
      ::MPI_Pack((void*)&obj, 1, MPI_LONG_DOUBLE,
		 ibuf.buffer_, ibuf.size_, &ibuf.pos_, ibuf.comm_);
      return ibuf;
  } // operator,

} // namespace empi

#endif // EMPI_MPI_IBUFFER_HPP

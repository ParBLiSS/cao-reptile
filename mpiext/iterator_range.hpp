/***
 *  $Id: iterator_range.hpp 715 2007-02-05 15:56:38Z zola $
 **
 *  File: iterator_range.hpp
 *  Developed: Feb 09, 2006
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2005-2007 Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  For copyright details please see attached LICENSE
 */

#ifndef EMPI_ITERATOR_RANGE_HPP
#define EMPI_ITERATOR_RANGE_HPP


/** extMPI main namespace.
 */
namespace empi {

  /** Structure to describe range of iterators.
   */
  template <typename Iter>
  struct iterator_range {
      /**
       */
      typedef Iter iterator_type;

      /** Constructs iterator_range [@a first, @a last).
       *  @param first is beginning of the range.
       *  @param last is end of the range.
       */
      iterator_range(iterator_type first, iterator_type last)
	  : begin(first), end(last) { }

      /** Beginning of the range.
       */
      iterator_type begin;
      
      /** End of the range.
       */
      iterator_type end;
      
  }; // struct iterator_range

} // namespace empi

#endif // EMPI_ITERATOR_RANGE_HPP

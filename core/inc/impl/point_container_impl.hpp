/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POINT_CONTAINER_IMPL_HPP_INCLUDED
#define POINT_CONTAINER_IMPL_HPP_INCLUDED

#include <type_traits>

namespace cellsim {

     namespace impl_{
	  typedef Point_Container P_C;
	  typedef P_C::iterator iterator_t;
	  typedef P_C::const_iterator const_iterator_t;

	  GCC_DIAG_OFF(return-type)

	  template <TYPE type, TIME time>
	  const_iterator_t cbegin(const P_C& /*c*/)
	  { 
	       static_assert(type == CA && time == NEXT, 
			     "Instantiation of cbegin not allowed !"); 
	  }
	  template <TYPE type, TIME time>
	  iterator_t begin(P_C& /*c*/)
	  { 
	       static_assert(type == CA && time == NEXT, 
			     "Instantiation of begin not allowed !"); 
	  }
	  template <TYPE type, TIME time>
	  const_iterator_t cend(const P_C& /*c*/)
	  { 
	       static_assert(type == CA && time == NEXT, 
			     "Instantiation of cend not allowed !"); 
	  }
	  template <TYPE type, TIME time>
	  iterator_t end(P_C& /*c*/)
	  { 
	       static_assert(type == CA && time == NEXT,
			     "Instantiation of cend not allowed !"); 
	  }
	  GCC_DIAG_ON(return-type)

	  template<>
	  inline const_iterator_t cbegin<CA, NOW>(const P_C& c) {return c.ca_begin_now();}
	  template<>
	  inline const_iterator_t cend<CA, NOW>(const P_C& c) {return c.ca_end_now();}
	  template<>
	  inline const_iterator_t cbegin<CA, NEXT>(const P_C& c) {return c.ca_begin_next();}
	  template<>
	  inline const_iterator_t cend<CA, NEXT>(const P_C& c) {return c.ca_end_next();}
	  template<>
	  inline iterator_t begin<CA, NEXT>(P_C& c) {return c.ca_begin_next();}
	  template<>
	  inline iterator_t end<CA, NEXT>(P_C& c) {return c.ca_end_next();}

	  template<>
	  inline const_iterator_t cbegin<IP3, NOW>(const P_C& c) {return c.ip3_begin_now();}
	  template<>
	  inline const_iterator_t cend<IP3, NOW>(const P_C& c) {return c.ip3_end_now();}
	  template<>
	  inline const_iterator_t cbegin<IP3, NEXT>(const P_C& c) {return c.ip3_begin_next();}
	  template<>
	  inline const_iterator_t cend<IP3, NEXT>(const P_C& c) {return c.ip3_end_next();}
	  template<>
	  inline iterator_t begin<IP3, NEXT>(P_C& c) {return c.ip3_begin_next();}
	  template<>
	  inline iterator_t end<IP3, NEXT>(P_C& c) {return c.ip3_end_next();}

	  template<>
	  inline const_iterator_t cbegin<S, NOW>(const P_C& c) {return c.s_begin_now();}
	  template<>
	  inline const_iterator_t cend<S, NOW>(const P_C& c) {return c.s_end_now();}
	  template<>
	  inline iterator_t begin<S, NEXT>(P_C& c) {return c.s_begin_next();}
	  template<>
	  inline iterator_t end<S, NEXT>(P_C& c) {return c.s_end_next();}
     } // namespace impl_

     using impl_::begin;
     using impl_::end;
     using impl_::cbegin;
     using impl_::cend;

} // namespace cellsim

#endif //POINT_CONTAINER_IMPL_HPP_INCLUDED

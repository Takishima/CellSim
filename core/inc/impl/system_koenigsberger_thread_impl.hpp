/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef SYSTEM_KOENIGSBERGER_THREAD_IMPL_HPP_INCLUDED
#define SYSTEM_KOENIGSBERGER_THREAD_IMPL_HPP_INCLUDED

namespace cellsim {
namespace impl_ {

     template <typename time_t,
	       typename function_t, 
	       typename activator_t,
	       typename iterator_t>
     void execute_thread(time_t t, 
			 time_t dt, 
			 function_t f, 
			 activator_t activator,
			 iterator_t it,
			 const iterator_t last)
     {
	  for(; it != last ; ++it) {
	       // f(*it, t, dt, activator_impl::no_op<double>());
	       f(*it, t, dt, activator);
	  }
     }

} // namespace impl_
} // namespace cellsim


#endif //SYSTEM_KOENIGSBERGER_THREAD_IMPL_HPP_INCLUDED

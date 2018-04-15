/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef NO_OP_HPP_INCLUDED
#define NO_OP_HPP_INCLUDED

#include "equation_type.hpp"
#include "definitions.hpp"

namespace cellsim {
namespace activator_impl {

     /*!
      * No-operation activation functor
      *
      * The purpose of this class is to have an activator that does nothing
      * when calling the do_step method of the system.
      *
      * Typically, the function operator of this class is never even called,
      * but required (need to fix that...)
      */
     template <typename return_type>
     struct no_op
     {
	  /*!
	   * This activator should not affect \e any equation.
	   */
	  typedef equation_type::NO_EQN equation_type;

	  constexpr bool do_activation(size_type /*i*/) {return false;}

	  /*!
	   * Simple no-op function operator.
	   * \return return_type()
	   */
	  template <typename T1, typename T2>
	  constexpr return_type operator()(const T1 /*t1*/, const T2 /*t2*/) const
	       {
	  	    return return_type();
	       }
	  /*!
	   * Simple no-op function operator.
	   * \return return_type()
	   */
	  template <typename T1>
	  constexpr return_type operator()(const T1 /*t1*/) const
	       {
	  	    return return_type();
	       }
     };

} // namespace activator_impl
} // namespace cellsim

#endif //NO_OP_HPP_INCLUDED

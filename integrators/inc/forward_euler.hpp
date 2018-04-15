/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef FORWARD_EULER_HPP_INCLUDED
#define FORWARD_EULER_HPP_INCLUDED

#include "equation_type.hpp"
#include "impl/rk4_impl.hpp"

#include <iomanip>		
#include <iostream>

namespace cellsim {
namespace integrators {

     /*!
      * \brief Forward-Euler integrator
      *
      * We are basically solving the following equation : \f[\frac{{\rm d}y}{{\rm d}\alpha} = f(y) \qquad \text{where $\alpha \in \lbrace x, t\rbrace$ typically}\f]
      *
      * \tparam equation_t The type of equation this integrator will integrate
      */
     template <typename equation_t>
     class Forward_Euler
     {
     public:
	  /*!
	   * Main method to do a Forward Euler step from a to a + da.
	   *
	   * \param y_now Value to integrate at a (can be a vector or a number)
	   * \param y_new Value to integrate at a+da (can be a vector or a number)
	   * \param a Value of integration variable at current step
	   * \param da Step in a.
	   * \param f Function to call
	   * \param mem_func Member function pointer to extract the value of
	   *                 interest from \c y_now.
	   * \param act_func Activator to use while integrating.\n
	   *                 Requirements on \c activator_t: \n
	   *                 \c return_type \c operator(const input_type);
	   *                 where \c return_type is the type of \c f(y_now, a)
	   */
	  template <typename input_type,
		    typename number_type,
		    typename variable_type,
		    typename incr_type,
		    typename functor,
		    typename mem_func_t,
		    typename activator_t = activator_impl::no_op<double> >
	  static void do_step(input_type y_now, 
			      number_type& y_new,
			      const variable_type a, 
			      incr_type& da, 
			      const functor f,
			      mem_func_t mem_func,
			      activator_t& act_func)
	       {
		    // using rk4_impl_::eval_if;
		    typedef typename activator_t::equation_type act_type;
		    using if_same_eq = 
			 eval_if<is_same_eq<equation_t, act_type>::value>;

		    y_new = mem_func(y_now) 
#ifdef NO_UNITS
			 + double(da * if_same_eq::apply_activator(f, 
								   y_now, a,
								   act_func));
#else
		    + da * if_same_eq::apply_activator(f, 
						       y_now, a,
						       act_func);
#endif
	       }

	  /*!
	   * Secondary method to do a Forward Euler step from a to a + da.
	   *
	   * \param y_now Value to integrate at a (can be a vector or a number)
	   * \param y_new Value to integrate at a+da (can be a vector or a number)
	   * \param a Value of integration variable at current step
	   * \param da Step in a.
	   * \param f Function to call
	   * \param act_func Activator to use while integrating.\n
	   *                 Requirements on \c activator_t: \n
	   *                 \c return_type \c operator(const input_type);
	   *                 where \c return_type is the type of \c f(y_now, a)
	   */
	  template <typename number_type,
		    typename variable_type,
		    typename incr_type, 
		    typename functor,
		    typename activator_t = activator_impl::no_op<double> >
	  static void do_step(const number_type y_now, 
			      number_type& y_new,
			      const variable_type a, 
			      incr_type& da, 
			      const functor f,
			      activator_t act_func =
			      activator_impl::no_op<double>())
	       {
		    // using rk4_impl_::eval_if;
		    // typedef typename activator_t::equation_type act_type;
		    // using if_same_eq = 
		    // 	 eval_if<is_same_eq<equation_t, act_type>::value>;
		    do_step(y_now, y_new, a, da, f,
			    &identity<number_type>, act_func);
	       }

     };
     
} // namespace integrators
} // namespace cellsim

#endif //FORWARD_EULER_HPP_INCLUDED

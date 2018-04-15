/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef RK4_IMPL_HPP_INCLUDED
#define RK4_IMPL_HPP_INCLUDED

#include <cmath>
#include <functional>

#include "equation_type.hpp"
#include "no_op.hpp"
#include "identity.hpp"
#include "apply_activator.hpp"

namespace cellsim {
namespace integrators {
namespace rk4_impl_ {
     /*!
      * Template metafunction to do 2 half RK4 steps and compare the precision
      * thus obtained to the full RK4 step.
      *
      * Whether to do the half steps or not is decided at \b compile \b time.
      *
      * \tparam adapt_da Whether we are using the adaptative scheme or not
      *                  (ie. do the half steps or not)
      * \tparam equation_t Equation type of the integrator
      */
     template <bool adapt_da, typename equation_t>
     struct is_adaptative
     {
	  template <typename input_type,
		    typename number_type, 
		    typename variable_type,
		    typename incr_type, 
		    typename functor,
		    typename mem_func_t,
		    typename activator_t>
	  static double apply(const input_type y_now, 
			      number_type& y_new,
			      const variable_type a, 
			      incr_type& da, 
			      const functor f,
			      const mem_func_t mem_func,
			      activator_t& act_func)
	       {
		    typedef typename activator_t::equation_type act_type;
		    using if_same_eq = 
			 eval_if<is_same_eq<equation_t, act_type>::value>;

		    const double da_2(da/2.), da_4(da/4.), half_a(a + da_2);

		    number_type y_half_a(0.), y_full_a(0.);
		    input_type tmp(y_now);

		    // Two half-step integration for the adaptative da
		    // k1 = f(y_now, a);
		    // k2 = f(y_now + .5*k1, a + da_4);
		    number_type k1(if_same_eq::apply_activator(
					f, y_now, a, act_func) * da_2);

		    mem_func(tmp) += 0.5 * k1;
		    number_type k2(if_same_eq::apply_activator(
					f, tmp, a + da_4, act_func) * da_2);

		    mem_func(tmp) = mem_func(y_now) + 0.5 * k2;
		    number_type k3(if_same_eq::apply_activator(
					f, tmp, a + da_4, act_func) * da_2);

		    mem_func(tmp) = mem_func(y_now) + k1;
		    number_type k4(if_same_eq::apply_activator(
					f, tmp, half_a, act_func) * da_2);

		    y_half_a = mem_func(y_now) + 1./6.*(k1 + 2.*k2 + 2.*k3 + k4);

		    mem_func(tmp) = mem_func(y_half_a);
		    k1 = if_same_eq::apply_activator(
			 f, tmp, half_a, act_func) * da_2;

		    mem_func(tmp) += 0.5 * k1;
		    k2 = if_same_eq::apply_activator(
			 f, tmp, half_a + da_4, act_func) * da_2;

		    mem_func(tmp) = mem_func(y_half_a) + 0.5 * k2;
		    k3 = if_same_eq::apply_activator(
			 f, tmp, half_a + da_4, act_func) * da_2;
		    
		    mem_func(tmp) = mem_func(y_half_a) + k1;
		    k4 = if_same_eq::apply_activator(
			 f, tmp, a + da, act_func) * da_2;

		    y_full_a = y_half_a + 1./12. * (k1 + 2.*k2 + 2*k3 + k4);
	  
		    return std::fabs(y_new - y_full_a);
	       }
     };

     /*!
      * Template specialisation in the case we are \b not using the adaptative
      * da scheme.
      */
     template <typename equation_t>
     struct is_adaptative<false, equation_t>
     {
	  template <typename input_type,
		    typename number_type,
		    typename variable_type,
		    typename incr_type, 
		    typename functor,
		    typename mem_func_t,
		    typename activator_t>
	  static constexpr double apply(const input_type /*y_now*/, 
					const number_type /*y_new*/,
					const variable_type /*a*/, 
					const incr_type /*da*/, 
					const functor /*f*/,
					const mem_func_t /*mem_func*/,
					activator_t /*act_func*/)
	       {
		    return 0.;
	       }
     };

} // namespace rk4_impl_
} // namespace integrators
} // namespace cellsim

#endif //RK4_IMPL_HPP_INCLUDED

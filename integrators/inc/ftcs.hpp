/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef FTCS_HPP_INCLUDED
#define FTCS_HPP_INCLUDED

CLANG_DIAG_OFF(sign-conversion)

#include "impl/ftcs_impl.hpp"

#include "point_index.hpp"

namespace cellsim {
     namespace integrators {

	  /*!
	   * \brief Integrator for the diffusion equation using a FTCS scheme.
	   *
	   * The equation we want to integrate here is of the form :  \f[\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2} + F(u, t)\f]
	   * where \f$ F(u, t) \f$ is an arbitrary function of \f$u\f$ at time \f$t\f$.
	   *
	   * Once discretized, we defined the following coefficients: \f[\alpha = \frac{D \Delta t}{\Delta x} \qquad \beta = \Delta t\f]
	   * 
	   * Thus we get the following expression :\f[u_i^{n+1} = \alpha\,u_{i-1}^n + (1 - 2\alpha)\,u_i^n + \alpha\,u_{i+1}^n + F\left(u_i^n\right)\f]
	   *
	   * \tparam equation_t The type of equation this integrator will integrate
	   */
	  template <typename equation_t, 
		    typename type_select = ftcs_impl_::default_t>
	  class FTCS
	  {
	       typedef typename 
	       ftcs_impl_::type_helper<type_select>::time_t time_t;
	       typedef typename 
	       ftcs_impl_::type_helper<type_select>::length_t length_t;
	       typedef typename 
	       ftcs_impl_::type_helper<type_select>::diff_coef_t diff_coef_t;
	  public:
	       /*!
		* Default Constructor
		*/
	       FTCS(time_t dt, length_t dx);

	       FTCS(diff_coef_t D, time_t dt, length_t dx);

	       /*!
		* Setup integrator
		*
		* \param D Diffusion coefficient
		*/
	       void setup_system(diff_coef_t D);
	       /*!
		* Setup integrator
		*
		* \param D Diffusion coefficient
		* \param dt Timestep
		* \param dx Spatial step
		*/
	       void setup_system(diff_coef_t D, time_t dt, length_t dx);


	       /*!
		* \brief Integrate the diffusion equation one time step using an FTCS scheme.
		*
		* \param t Current time
		* \param y_now Array of values at time \f$ t \f$.
		* \param out_it Output iterator (pointing to the output array)
		* \param F Function of the values of \c y_now as defined above.
		* \param act_func Activation functor.\n
		*                 For each element of \c y_now, when evaluating the F term,
		*                 the algorithm will add the return value of the function
		*                 operator of \c act_func applied at that element.\n
		*                 Thus this functor will dictate which point is activated and
		*                 how it will be activated.\n
		*                 The function operator () should have the following signature :
		*                 \c ret_type \c operator(const v_type v, size_type idx) const;  
		*                 where \c ret_type is the same as the return value of F.
		*                 v_type is the type of element inside \c y_now.\n
		*                 By default, the no_op functor is used, which does nothing 
		*                 (ie return 0).
		*
		* \return Array of values at time \f$ t + \Delta t\f$.
		*/
	       template <typename array_type, 
			 typename out_iterator_t,
			 typename function,
			 typename act_functor = activator_impl::no_op<double> >
	       void do_step(const time_t t,
			    const array_type y_now,
			    out_iterator_t out_it,
			    function F,
			    act_functor& act_func) const;

	       /*!
		* \deprecated This function is not up-to-date and should'nt be used as is.
		*
		* \brief Convenience "overload"
		*
		* Stores the result directly inside y_now by using \c mem_func
		* function pointer (typically a pointer to a member function of 
		* array_type::value_type).
		*
		* \param t Current time
		* \param y_now Array of values at time \f$ t \f$.
		* \param F Function of the values of \c y_now as defined above.
		* \param mem_func Function pointer to a member of array_type::value_type
		*                 to get the value of interest for each value of \c y_now.
		* \param act_func Activation functor.\n
		*                 See do_step() for more info on this last parameter
		*/
	       template <typename mem_func_t,
	       		 typename function,
	       		 typename array_type, 
	       		 typename act_functor = activator_impl::no_op<double> >
	       void do_step_in_place(const time_t t,
	       			     array_type& y_now,
	       			     const function F,
	       			     const mem_func_t mem_func,
	       			     act_functor& act_func) const;
	       
	  private:
	       time_t dt_; //!< Timestep
	       length_t dx_; //!< Spatial step
	       diff_coef_t D_; //!< Diffusion coefficient
	  };

	  // Constructors
	  template <typename equation_t,
		    typename type_select>
	  FTCS<equation_t, type_select>::FTCS(time_t dt, length_t dx)
	       : FTCS(diff_coef_t(), dt, dx)
	  {}

	  template <typename equation_t,
		    typename type_select>
	  FTCS<equation_t,
	       type_select>::FTCS(diff_coef_t D, time_t dt, length_t dx)
	       : dt_(dt), dx_(dx), D_(D)
	  {}
     
	  // setup_system
	  template <typename equation_t,
		    typename type_select>
	  void FTCS<equation_t, type_select>::setup_system(diff_coef_t D)
	  {
	       D_ = D;
	  }
	  template <typename equation_t,
		    typename type_select>
	  void FTCS<equation_t, 
		    type_select>::setup_system(diff_coef_t D, 
					       time_t dt, 
					       length_t dx)
	  {
	       D_ = D;
	       dt_ = dt;
	       dx_ = dx;
	  }


	  // do_step
	  template <typename equation_t,
		    typename type_select>
	  template <typename array_type, 
		    typename out_iterator_t,
		    typename function,
		    typename act_functor>
	  void FTCS<equation_t, 
		    type_select>::do_step(time_t t,
					  const array_type y_now,
					  out_iterator_t out_it,
					  function F,
					  act_functor& act_func) const
	  {
	       // using ftcs_impl_::eval_if;
	       typedef typename act_functor::equation_type act_type;
	       typedef eval_if<is_same_eq<equation_t, act_type>::value> if_same_eq;

	       const auto alpha(dt_ * D_ / (dx_ * dx_));
	       const auto& beta(dt_);
	       const auto factor(1. - 2. * alpha);

	       Point_Index index(&y_now);
	  
	       // Leave first and last point intact
	       const auto BEGIN(cbegin<equation_t::value, NOW>(y_now)), 
		    END(cend<equation_t::value, NOW>(y_now));

	       for (auto prev(BEGIN-1), it(BEGIN), next(BEGIN+1) ; it != END ; 
	       	    ++prev, ++it, ++next, ++out_it, ++index) {
	       	    using namespace std::placeholders;

		    // Standard FTCS scheme
		    *out_it = 
			 alpha * (*prev)     // u(i-1,n)
			 + factor * (*it)     // u(i,n)
			 + alpha * (*next);  // u(i+1,n)

		    GCC_DIAG_OFF(sign-compare)
		    assert(index.index() == it-BEGIN);
		    GCC_DIAG_ON(sign-compare)

	       	    /*
		     * Additionnal F term
		     *
	       	     * Only use activator if we are integrating the equation
	       	     * corresponding to the activator functor.
	       	     * This is done at COMPILE TIME
	       	     */
	       	    *out_it += beta * if_same_eq::apply_activator(F,
								  index,
								  t,
								  act_func);
	       	    // *out_it += beta * if_same_eq::apply_activator(std::bind(F, _1, t),
		    // 						  index,
		    // 						  act_func);
		    assert(out_it->value() < 1e30 && out_it->value() > -1e30);
	       }
	  }
     } // namespace integrators
} // namespace cellsim

CLANG_DIAG_ON(sign-conversion)

#endif //FTCS_HPP_INCLUDED

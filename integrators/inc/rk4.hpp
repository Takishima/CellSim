/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef RK4_HPP_INCLUDED
#define RK4_HPP_INCLUDED

#include "impl/rk4_impl.hpp"

namespace cellsim {
namespace integrators {
     /*!
      * \brief Runge-Kutta 4th order integrator.
      *
      * We are basically solving the following equation : \f[\frac{{\rm d}y}{{\rm d}\alpha} = f(y) \qquad \text{where $\alpha \in \lbrace x, t\rbrace$ typically}\f]
      *
      * \tparam use_adaptative_ Template parameter controlling whether we are using
      *                         an adaptative step scheme or not.
      * \tparam equation_t The type of equation this integrator will integrate
      */
     template <bool use_adaptative_, typename equation_t>
     class RK4_Integrator
     {
	  template <bool val>
	  using is_adaptative = rk4_impl_::is_adaptative<val, equation_t>;
     public:
	  /*!
	   * Main method to do an RK4 step from a to a + da.
	   *
	   * \param y_now Value to integrate at a (can be a vector or a number)
	   * \param y_new Value to integrate at a+da (can be a vector or a number)
	   * \param a Value of integration variable at current step
	   * \param da Step in a.
	   * \param f Function to call
	   * \param mem_func Member function to call on \c y_now to read the value
	   *                 of the quantity to integrate.
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
	  static void do_step(const input_type y_now, 
			      number_type& y_new,
			      const variable_type a, 
			      incr_type& da, 
			      const functor f,
			      const mem_func_t mem_func,
			      activator_t& act_func)
	       {
		    // if (use_adaptative_) {
		    // 	 double old_da(da.value());
		    // 	 double deviation(0.);
		    // 	 bool da_precise_enough(false);
            
		    // 	 for(unsigned int iter_num(0) ;
		    // 	     iter_num < MAX_ITER_ && !da_precise_enough ;
		    // 	     ++iter_num){
                
		    // 	      // Do one RK4 step and compute deviation
		    // 	      deviation = do_rk4_step_compare_(y_now, 
		    // 					       y_new, 
		    // 					       a,
		    // 					       da, 
		    // 					       f, 
		    // 					       mem_func,
		    // 					       act_func);

		    // 	      da_precise_enough = (deviation < precision_);

		    // 	      // Adapt dt for next timestep or next try
		    // 	      adapte_da(da.value(), old_da, deviation, da_precise_enough);
		    // 	 }
		    // }
		    // else {
			 do_rk4_step_compare_(y_now, y_new, a, da, f, mem_func, act_func);
		    // }
	       }

	  /*!
	   * Overload to use with plain numbers (as opposed to the other one
	   * that assumes that \c y_now is in fact a class or struct.
	   * \see do_step
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
			      activator_t& act_func)
	       {
		    typedef std::function<number_type&(number_type&)> func_t;
		    func_t dummy_func(&identity<number_type>);

		    // if (use_adaptative_) {
		    // 	 double old_da(da.value());
		    // 	 double deviation(0.);
		    // 	 bool da_precise_enough(false);
            
		    // 	 for(unsigned int iter_num(0) ;
		    // 	     iter_num < MAX_ITER_ && !da_precise_enough ;
		    // 	     ++iter_num){
                
		    // 	      // Do one RK4 step and compute deviation
		    // 	      deviation = do_rk4_step_compare_(y_now, 
		    // 					       y_new, 
		    // 					       a,
		    // 					       da, 
		    // 					       f,
		    // 					       dummy_func,
		    // 					       act_func);

		    // 	      da_precise_enough = (deviation < precision_);

		    // 	      // Adapt dt for next timestep or next try
		    // 	      adapte_da(quantity_cast<double>(da), old_da, deviation, da_precise_enough);
		    // 	 }
		    // }
		    // else {
			 do_rk4_step_compare_(y_now, y_new, a, da, f,
					      dummy_func, act_func);
		    // }
	       }

	  /*!
	   * Set precision required based on component value of vec
	   *
	   * \param vec Array of value to consider
	   */
	  template <typename vector_type>
	  void set_da_precision(const vector_type vec) 
	       {
		    precision_ = pow( 10, floor(log10(vec.max())) - 6 );
	       }

     private:

	  /*!
	   * Do an RK4 step for arrays.
	   *
	   * If using the adaptative scheme, then do one full RK4 step followed by
	   * two half RK4 steps and compare the deviation between the two results.
	   *
	   * \param y_now Value to integrate at a
	   * \param y_new Value to integrate at a+da
	   * \param a Value of integration variable at current step
	   * \param da Step in a.
	   * \param f Function to call
	   * \param mem_func Member function to call on \c y_now to read the value
	   *                 of the quantity to integrate.
	   * \param act_func Activator to use while integrating.\n
	   *                 Requirements on \c activator_t: \n
	   *                 \c return_type \c operator(const input_type);
	   *                 where \c return_type is the type of \c f(y_now, a)
	   * \return Deviation between full and two half steps if using adaptive da
	   *         0 otherwise.
	   */
	  template <typename input_type,
		    typename number_type, 
		    typename variable_type,
		    typename incr_type, 
		    typename functor,
		    typename mem_func_t,
		    typename activator_t = activator_impl::no_op<double> >
	  static double do_rk4_step_compare_(input_type y_now, 
					     number_type& y_new,
					     const variable_type a, 
					     incr_type& da, 
					     const functor f,
					     const mem_func_t mem_func,
					     activator_t& act_func)
	       {
	  	    // using rk4_impl_::eval_if;
	  	    typedef typename activator_t::equation_type act_type;
	  	    using if_same_eq = 
	  		 eval_if<is_same_eq<equation_t, act_type>::value>;

	  	    const variable_type half_a(a + 0.5 * da);
		    input_type tmp(y_now);

	  	    number_type k1(if_same_eq::apply_activator(f, tmp, a, 
	  	    					       act_func) * da);
	  	    mem_func(tmp) += 0.5 * k1;
	  	    number_type k2(if_same_eq::apply_activator(
					f, tmp, half_a, act_func) * da);

	  	    mem_func(tmp) = mem_func(y_now) + 0.5 * k2;
	  	    number_type k3(if_same_eq::apply_activator(f, tmp, half_a,
							       act_func) * da);
	  	    mem_func(tmp) = mem_func(y_now) + k3;
	  	    number_type k4(if_same_eq::apply_activator(f, tmp, a + da, 
							       act_func) * da);

	  	    y_new = mem_func(y_now) + 1./6.*(k1 + 2.*k2 + 2.*k3 + k4);

	  	    /*
	  	     * Adaptative da algorithm is only executed if required
	  	     * Note that the choice is made at compile time so the
	  	     * compiler can optimize the code
	  	     */
	  	    return is_adaptative<use_adaptative_>::apply(
	  	    	 y_now, y_new, a, da, f, mem_func, act_func);
	       }

	  /*!
	   * \deprecated
	   * Calculate new value for da based on the deviation and precision.
	   * 
	   * We try not to change da too abruptly. 
	   * Also avoid too big or too small da.
	   *
	   * \param da Current da
	   * \param old_da Old value of da
	   * \param deviation Deviation between one full RK4 step and two half RK4 steps
	   * \param da_precise_enough Boolean indicating whether da gives a
	   *                          good enough precision.
	   */
	  // static void adapte_da(double& da, const double old_da,
	  // 			const double deviation,
	  // 			bool& da_precise_enough)
	  //      {
	  // 	    const double da_new( CORRECTION_FACTOR_ * da * 
	  // 				 pow((precision_/deviation), 1./5.)
	  // 		 );

	  // 	    if (da_new < precision_) {
	  // 		 // Avoid too small da
	  // 		 da = precision_;
	  // 	    }
	  // 	    else if (da_new > 1e5 * precision_) {
	  // 		 // Avoid too big da
	  // 		 da = 1e5*precision_;
	  // 	    }
	  // 	    else if (da_new < 0.01 * old_da) {
	  // 		 // Avoid too sudden changes of da
	  // 		 da = 0.1 * old_da;
	  // 		 da_precise_enough = true;
	  // 	    }
	  // 	    else if (da_new > 100 * old_da) {
	  // 		 // Avoid too sudden changes of da
	  // 		 da = 100 * old_da;
	  // 		 da_precise_enough = true;
	  // 	    }
	  // 	    else {
	  // 		 da = da_new;
	  // 	    }
	  //      }

	  //! Maximum number of iteration when using adaptative da
	  static constexpr unsigned int MAX_ITER_ = 10;
	  //! RK4 adaptative da correction factor
	  static constexpr double CORRECTION_FACTOR_ = 0.92;
	  //! Precision required on solution at a+da
	  static double precision_;
     };

     //! Convenience typedef
     template <typename equation_t>
     using Adaptative_RK4 = RK4_Integrator<true, equation_t>;
     //! Convenience typedef
     template <typename equation_t>
     using Standard_RK4 = RK4_Integrator<false, equation_t>;

     template <bool b, typename equation_t> 
     double RK4_Integrator<b, equation_t>::precision_ = 1e-5;
     
} // namespace integrators
} // namespace cellsim

#endif // RK4_HPP_INCLUDED

/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTENTIAL_ACTIVATOR_BASE_HPP_INCLUDED
#define POTENTIAL_ACTIVATOR_BASE_HPP_INCLUDED

#include "equation_type.hpp"
#include "units.hpp"

namespace cellsim {
     GCC_DIAG_OFF(effc++)

     /*!
      * \brief Base class for potential activators
      *
      * This class implements the basic interface of membrane potential 
      * activators.
      */
     class Potential_Activator_Base
     {
     public:
	  /*!
	   * This activator is designed to enforce the value of one of the 
	   * dynamical variables. Therefore, it does not apply to any equation
	   */
	  typedef equation_type::NO_EQN equation_type;

	  /*!
	   * \brief Simple constructor
	   *
	   * \param v_act Activation value of the potential
	   * \param t0 Time during which to hold maximum value before decreasing
	   * \param alpha Decay constant
	   */
	  Potential_Activator_Base(electric_potential_ut v_act,
				   time_ut t0,
				   double alpha);
	  
	  //! Simple accessor
	  electric_potential_ut value() const { return v_act_; }

	  /*!
	   * \brief Function operator
	   *
	   * This operator expects the following to be valid:
	   * \li Cell must have a electric_potential_ut v() const method
	   * \li Cell must have a electric_potential_ut& v_new() method
	   *
	   * \param c Cell to activate
	   * \param dt Timestep
	   */
	  template <typename Cell>
	  void operator() (Cell& c, time_ut dt)
	       {
		    if (v0_ == 0_mV) {
			 v0_ = c.v();
		    }

		    c.v_new() = v0_ + v_act_;
		    if (t0_ < t_) {
			 v_act_ *= (1. - alpha_ * dt.value());
		    }
		    t_ += dt.value();	 
	       }

     protected:
	  electric_potential_ut v0_; //!< Initial value of the eletric potential
	  electric_potential_ut v_act_; //!< Activation value for the membrane potential
	  const double alpha_; //! Decay constant
	  double t_; //! Current time of activator (from first use)
	  double t0_; //! Time to hold activation constant
     };
     
     GCC_DIAG_ON(effc++)

} // namespace cellsim

#endif //POTENTIAL_ACTIVATOR_BASE_HPP_INCLUDED

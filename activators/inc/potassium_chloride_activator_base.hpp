/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTASSIUM_CHLORIDE_ACTIVATOR_BASE_HPP_INCLUDED
#define POTASSIUM_CHLORIDE_ACTIVATOR_BASE_HPP_INCLUDED

#include "equation_type.hpp"
#include "units.hpp"

#include <iostream>

namespace cellsim {
namespace activator {
     GCC_DIAG_OFF(effc++)
     GCC_DIAG_OFF(padded)
     
     template <typename model_constants_t>
     class K_Cl_Activator_Base
     {
     public:
	  /*!
	   * This activator is designed to modify some values of the model
	   * parameters.
	   */
	  typedef equation_type::MODEL_PARAM_EQN equation_type;

	  K_Cl_Activator_Base(model_constants_t constants,
			      time_ut t0,
			      double alpha_rise,
			      double alpha_decay,
			      electric_potential_ut vK,
			      electric_potential_ut vCl)
	       : constants_(constants), vK_(vK), vCl_(vCl),
		 old_vK_(constants.vK), old_vCl_(constants.vCl),
		 t_(-4. * alpha_rise), 
		 t0_(t0.value()),
		 t_fin_(4. * alpha_rise),
		 alpha_rise_(alpha_rise), alpha_decay_(alpha_decay)
	       {}

	  model_constants_t get_constants(time_ut dt)
	       {
		    const auto delta_K = vK_ - old_vK_;
		    const auto delta_Cl = vCl_ - old_vCl_;
		    
		    if (t_ < t_fin_) {
			 constants_.vK = old_vK_ + 0.5 * delta_K * (std::tanh(t_ / alpha_rise_)+1);
			 constants_.vCl = old_vCl_ + 0.5 * delta_Cl * (std::tanh(t_ / alpha_rise_)+1);
		    }
		    else {
			 // After Tanh increase, holding values
			 constants_.vK = vK_;
			 constants_.vCl = vCl_;
		    }

		    // Tanh decay after holding
		    if ((t_fin_ + t0_) < t_) {
			 const auto tanh_factor = 
			      (1 - std::tanh((t_- t0_ - t_fin_ - 4*alpha_decay_)
					     / alpha_decay_));
		    	 constants_.vK = old_vK_ + 0.5 * delta_K * tanh_factor;
		    	 constants_.vCl = old_vCl_ + 0.5 * delta_Cl * tanh_factor;

		    }
		    t_ += dt.value();
		    return constants_;
	       }

     private:
	  model_constants_t constants_; //!< Model constants
	  const electric_potential_ut vK_; //!< Activation value for vK
	  const electric_potential_ut vCl_; //!< Activation value for vCl
	  const electric_potential_ut old_vK_;
	  const electric_potential_ut old_vCl_;
	  double t_; //! Starting time for tanh function
	  const double t0_; //! Time to hold activation constant
	  const double t_fin_; //! End time for tanh function
	  const double alpha_rise_; //! Rise constant for tanh
	  const double alpha_decay_; //! Decay constant
     };

     GCC_DIAG_ON(padded)
     GCC_DIAG_ON(effc++)
} // namespace activators
} // namespace cellsim

#endif //POTASSIUM_CHLORIDE_ACTIVATOR_BASE_HPP_INCLUDED

/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTENTIAL_ACTIVATOR_TANH_BASE_HPP_INCLUDED
#define POTENTIAL_ACTIVATOR_TANH_BASE_HPP_INCLUDED

#include "potential_activator_base.hpp"
#include "units.hpp"

namespace cellsim {
namespace activator {

     GCC_DIAG_OFF(effc++)
     class Potential_Activator_Tanh_Base 
	  : public Potential_Activator_Base
     {
     public:
	  Potential_Activator_Tanh_Base(electric_potential_ut v_act,
					time_ut t0,
					double alpha_rise,
					double alpha_decay)
	       : Potential_Activator_Base(v_act, t0, alpha_decay),
		 t_(-4. * alpha_rise), 
		 t_fin_(4. * alpha_rise),
		 alpha_rise_(alpha_rise)
	       {}
	  
	  template <typename cell_t>
	  void operator() (cell_t& c, time_ut dt) 
	       {
		    if (v0_ == 0_mV) {
			 v0_ = c.v();
		    }
     
		    if (t_ < t_fin_) {
			 c.v_new() = v0_ + 0.5 * v_act_ * (std::tanh(t_ / alpha_rise_)+1);
		    }
		    else {
			 // After Tanh increase, holding v_act_ for t0_
			 c.v_new() = v0_ + v_act_;  
		    }
     
		    // Tanh decay after holding
		    if ((t_fin_ + t0_) < t_) {
			 const auto shifted_t = (t_- t0_ - t_fin_ - 4*alpha_);
			 c.v_new() = v0_ + 0.5 * v_act_ 
			      * (1 - std::tanh(shifted_t / alpha_));
		    }
		    t_ += dt.value();
	       }
     private:
	  double t_; //! Starting time for tanh function
	  double t_fin_; //! End time for tanh function rise
	  double alpha_rise_; //! Rise constant for tanh
     };
     GCC_DIAG_ON(effc++)

} // namespace activator
} // namespace cellsim

#endif //POTENTIAL_ACTIVATOR_TANH_BASE_HPP_INCLUDED

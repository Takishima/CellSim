/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CONSTANTS_HALIDI_HPP_INCLUDED
#define CONSTANTS_HALIDI_HPP_INCLUDED

#include "config.hpp"
#include "units.hpp"

namespace cellsim {
     namespace model_halidi {
	  GCC_DIAG_OFF(padded)
	  GCC_DIAG_OFF(effc++)

	  struct Constants_Values
	  {
	       Constants_Values() = default;
	       Constants_Values(Constants_Values&&) = default;
	       Constants_Values& operator=(Constants_Values&&) = default;
	       Constants_Values(const Constants_Values&) = default;
	       Constants_Values& operator=(const Constants_Values&) = default;
		    
	       flux_ut F;
	       concentration_ut KI;
	       conductance_ut GCa;
	       electric_potential_ut vCa1;
	       electric_potential_ut vCa2;
	       electric_potential_ut RCa;
	       conductance_ut GNaCa;
	       concentration_ut cNaCa;
	       electric_potential_ut vNaCa;
	       flux_ut B;
	       concentration_ut cb;
	       flux_ut C;
	       concentration_ut sc;
	       concentration_ut cc;
	       inverse_time_ut D;
	       electric_potential_ut vd;
	       electric_potential_ut Rd;
	       inverse_time_ut L;
	       gamma_ut gamma;
	       flux_ut FNaK;
	       conductance_ut GCl;
	       concentration_ut cCl;
	       electric_potential_ut vCl;
	       conductance_ut GK;
	       electric_potential_ut vK;
	       inverse_time_ut lambda;
	       concentration_ut cw;
	       beta_ut beta;
	       electric_potential_ut vCa3;
	       electric_potential_ut RK;
	       conductance_ut Gback;
	       electric_potential_ut vrest;
	       inverse_time_ut k;
	       flux_ut E;
	       concentration_ut KCa;
	       inverse_time_ut g;

	       diff_coef_ut Dip;
	       diff_coef_ut Dc;
	       permeablity_ut Pip;
	       permeablity_ut Pc;
	  };

	  typedef const Constants_Values& Constants_f_param_t;

	  GCC_DIAG_ON(effc++)
	  GCC_DIAG_ON(padded)
     } // namespace model_halidi
} // namespace cellsim


#endif //CONSTANTS_HALIDI_HPP_INCLUDED

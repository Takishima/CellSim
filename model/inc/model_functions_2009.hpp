/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef MODEL_FUNCTIONS_2009_HPP_INCLUDED
#define MODEL_FUNCTIONS_2009_HPP_INCLUDED

#include "constants_halidi.hpp"
#include "units.hpp"

#include <cmath>

namespace cellsim {
     namespace model_halidi {	  
	  namespace model_functions_2009 {
	       flux_ut JVOCC(electric_potential_ut v, Constants_f_param_t C);
	       flux_ut JNaCa(concentration_ut c, 
			     electric_potential_ut v,
			     Constants_f_param_t C);
	       flux_ut JSRuptake(concentration_ut c, Constants_f_param_t C);
	       flux_ut JCICR(concentration_ut c, 
			     concentration_ut s, 
			     Constants_f_param_t C);
	       flux_ut Jextrusion(concentration_ut c, 
				  electric_potential_ut v, 
				  Constants_f_param_t C);
	       flux_ut Jleak(concentration_ut s, Constants_f_param_t C);
	       flux_ut JPLCd(concentration_ut c, Constants_f_param_t C);
	       flux_ut Jdegrad(concentration_ut I, Constants_f_param_t C);
	       flux_ut JIP3(concentration_ut I, Constants_f_param_t C);
	       flux_ut JNaK(Constants_f_param_t C);
	       flux_ut JCl(electric_potential_ut v, Constants_f_param_t C);
	       flux_ut JK(electric_potential_ut v, 
			  double w, 
			  Constants_f_param_t C);
	       double Kactiv(concentration_ut c, 
			     electric_potential_ut v, 
			     Constants_f_param_t C);
	  } // namespace model_functions_2009
     } // namespace model_halidi
} // namespace cellsim



#endif //MODEL_FUNCTIONS_2009_HPP_INCLUDED

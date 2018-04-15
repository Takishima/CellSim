/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef MODEL_FUNCTIONS_KOENIGSBERGER_HPP_INCLUDED
#define MODEL_FUNCTIONS_KOENIGSBERGER_HPP_INCLUDED

#include "constants_koenigsberger.hpp"
#include "units.hpp"

namespace cellsim {
     namespace model_koenigsberger {
	  class Cell;

	  namespace model_functions {
	       flux_ut JVOCC(cellsim::electric_potential_ut v, 
			     Constants_f_param_t C);
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
	       flux_ut Jdegrad(concentration_ut I, Constants_f_param_t C);
	       flux_ut JIP3(concentration_ut I, Constants_f_param_t C);
	       flux_ut JNaK(Constants_f_param_t C);
	       flux_ut JCl(concentration_ut c, 
			   electric_potential_ut v, 
			   Constants_f_param_t C);
	       flux_ut Jback(electric_potential_ut v, Constants_f_param_t C);
	       flux_ut JK(electric_potential_ut v, 
			  double w, 
			  Constants_f_param_t C);
	       double Kactiv(concentration_ut c, 
			     electric_potential_ut v, 
			     Constants_f_param_t C);
	       potential_SR_ut Vcoupling(const Cell c, Constants_f_param_t C);
	  } // namespace model_functions
     } // namespace model_koenigsberger
} // namespace cellsim



#endif //MODEL_FUNCTIONS_KOENIGSBERGER_HPP_INCLUDED

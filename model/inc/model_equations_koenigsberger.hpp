/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef MODEL_EQUATIONS_KOENIGSBERGER_HPP_INCLUDED
#define MODEL_EQUATIONS_KOENIGSBERGER_HPP_INCLUDED

#include "constants_koenigsberger.hpp"
#include "model_functions_koenigsberger.hpp"
#include "units.hpp"

namespace cellsim {
     namespace model_koenigsberger {
	  class Cell;

	  namespace model_equations {
	       /*!
		* Differential equation for the calcium concentration
		*
		* \param c A cell
		* \param C Model constants values
		* \param t Current time
		*/
	       flux_ut dc_dt(const Cell c, 
			     Constants_f_param_t C,
			     time_ut t);
	       /*!
		* Differential equation for the IP3 concentration
		*
		* \param c A cell
		* \param C Model constants values
		* \param t Current time
		*/
	       flux_ut dI_dt(const Cell c, 
			     Constants_f_param_t C,
			     time_ut t);
	       /*!
		* Differential equation for the calcium concentration in the SR
		*
		* \param c A cell
		* \param C Model constants values
		* \param t Current time
		*/
	       flux_ut ds_dt(const Cell c, 
			     Constants_f_param_t C,
			     time_ut t);
	       /*!
		* Differential equation for the cell membrane potential
		*
		* \param c A cell
		* \param C Model constants values
		* \param t Current time
		*/
	       potential_SR_ut dv_dt(const Cell c, 
				     Constants_f_param_t C,
				     time_ut t);
	       /*!
		* Differential equation for the Open state probability
		* of activated potassium channels
		*
		* \param c A cell
		* \param C Model constants values
		* \param t Current time
		*/
	       inverse_time_ut dw_dt(const Cell c, 
				     Constants_f_param_t C,
				     time_ut t);
	  } // namespace model_equations
     } // namespace model_koenigsberger
} // namespace cellsim


#endif //MODEL_EQUATIONS_KOENIGSBERGER_HPP_INCLUDED

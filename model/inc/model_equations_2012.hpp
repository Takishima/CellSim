/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef MODEL_EQUATIONS_2012_HPP_INCLUDED
#define MODEL_EQUATIONS_2012_HPP_INCLUDED

#include "constants_halidi.hpp"
#include "model_functions_2012.hpp"
#include "point_index.hpp"

#include "point.hpp"

namespace cellsim {
     namespace model_halidi {
	  class Cell;

	  namespace model_equations_2012 {
	       /*!
		* Differential equation for the calcium concentration
		*
		* \param point Point_Index to container with current data
		* \param v Cell membrane potentialmo
		* \param C Model constants
		* \param t Current time
		*/
	       flux_ut dc_dt(const Point_Index point, 
			     electric_potential_ut v, 
			     Constants_f_param_t C,
			     time_ut t);
	       /*!
		* Differential equation for the IP3 concentration
		*
		* \param point Point_Index to container with current data
		* \param JPLCago JPLCago concentration inside the cell
		* \param C Model constants
		* \param t Current time
		*/
	       flux_ut dI_dt(const Point_Index point, 
			     flux_ut JPLCago, 
			     Constants_f_param_t C,
			     time_ut t);
	       /*!
		* Differential equation for the calcium concentration in the SR
		*
		* \param point Point_Index to container with current data
		* \param C Model constants
		* \param t Current time
		*/
	       flux_ut ds_dt(const Point_Index point, 
			     Constants_f_param_t C,
			     time_ut t);

	       /*!
		* Differential equation for the cell membrane potential
		*
		* \param c A cell
		* \param C Model constants
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
		* \param C Model constants
		* \param t Current time
		*/
	       inverse_time_ut dw_dt(const Cell c,
				     Constants_f_param_t C,
				     time_ut t);
	  } // namespace model_equations_2012
     } // namespace model_halidi
} // namespace cellsim

#endif //MODEL_EQUATIONS_2012_HPP_INCLUDED

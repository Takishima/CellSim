/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTASSIUM_CHLORIDE_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED
#define POTASSIUM_CHLORIDE_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED

#include "potassium_chloride_activator_base.hpp"
#include "koenigsberger_activator_base.hpp"
#include "constants_koenigsberger.hpp"

namespace cellsim {
namespace model_koenigsberger {
namespace activator {
     GCC_DIAG_OFF(effc++)
     GCC_DIAG_OFF(padded)

     /*!
      * \brief Potassium and Cloride activator
      *
      * This class is designed to activate a cell by modifying the values
      * of the vK and vCl model constants.
      */
     class K_Cl_Activator
	  : public Activator_Base,
	    public cellsim::activator::K_Cl_Activator_Base<Constants_Values>
     {
     public:

	  /*!
	   * \brief Simple constructor
	   *
	   * The resulting Potential_Activator will activate columns in range
	   * [start, end) (ie. end non inclusive)
	   *
	   * \param start_column Index of first column to activate
	   * \param end_column Index of last column (non-inclusive) to activate
	   * \param v_act Initial value of potential perturbation 
	   *              (see \f$A\f$ in equation above)
	   * \param t0 Time after which the perturbation has an exponential decay
	   * \param alpha_rise Rise constant
	   * \param alpha_decay Decay constant
	   */
	  K_Cl_Activator(size_type start_column,
			 size_type end_column,
			 time_ut t0,
			 double alpha_rise,
			 double alpha_decay,
			 const Constants_Values& constants,
			 electric_potential_ut vK,
			 electric_potential_ut vCl);	  
     };
     
     GCC_DIAG_ON(padded)
     GCC_DIAG_ON(effc++)
} // namespace activator
} // namespace model_koenigsberger
} // namespace cellsim

#endif //POTASSIUM_CHLORIDE_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED

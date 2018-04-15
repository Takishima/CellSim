/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTASSIUM_CHLORIDE_ACTIVATOR_HPP_INCLUDED
#define POTASSIUM_CHLORIDE_ACTIVATOR_HPP_INCLUDED

#include "halidi_activator_base.hpp"
#include "potassium_chloride_activator_base.hpp"
#include "constants_halidi.hpp"

namespace cellsim {
namespace model_halidi {
namespace activator {
     GCC_DIAG_OFF(effc++)

     /*!
      * \brief Potassium and Cloride activator
      *
      * This class is designed to activate a cell by modifying the values
      * of the vK and vCl model constants.
      *
      * \note This class' purpose is to activate a \e unique cell
      *       each time the activator is called (= each timestep)
      */
     class K_Cl_Activator 
	  : public Activator_Base,
	    public cellsim::activator::K_Cl_Activator_Base<Constants_Values>
     {
     public:
	  /*!
	   * Simple constructor
	   *
	   * \param cell_idx Index of cell to activate
	   * \param t0 Time after which the perturbation has an exponential decay
	   * \param v_act Initial value of potential poerturbation 
	   *              (see \f$A\f$ in equation above)
	   * \param t0 Time after which the perturbation decays
	   * \param alpha_rise Rise constant
	   * \param alpha_decay Decay constant
	   */
	  K_Cl_Activator(size_type cell_idx,
			 time_ut t0,
			 double alpha_rise,
			 double alpha_decay,
			 const Constants_Values& constants,
			 electric_potential_ut vK,
			 electric_potential_ut vCl);	  
     };

     GCC_DIAG_ON(effc++)
} // namespace activator
} // namespace model_halidi
} // namespace cellsim

#endif //POTASSIUM_CHLORIDE_ACTIVATOR_HPP_INCLUDED

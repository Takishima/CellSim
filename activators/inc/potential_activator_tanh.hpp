/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTENTIAL_ACTIVATOR_TANH_HPP_INCLUDED
#define POTENTIAL_ACTIVATOR_TANH_HPP_INCLUDED

#include "halidi_activator_base.hpp"
#include "potential_activator_tanh_base.hpp"

namespace cellsim {
namespace model_halidi {
namespace activator {
     GCC_DIAG_OFF(effc++)

     /*!
      * This class purpose is to have a smoother increase in membrane potential
      * than a step function.
      * The equation for the rise is given by :
      * \f[ v_0 + 0.5 \, v_{\text{act}} \,\left[\tanh\left(\frac{t}{\alpha}\right)+1\right] \f]
      *
      * Then, there is a holding given by:
      * \f[ v_0 + v_{\text{act}} \quad \text{ for }t_0\f]
      *
      * And finally, the exponential decrease:
      * \f[  A\,e^{-\alpha t} \f]
      * 
      */
     class Potential_Activator_Tanh 
	  : public Activator_Base,
	    public cellsim::activator::Potential_Activator_Tanh_Base 
     {
     public:
	  /*!
	   * Simple constructor
	   *
	   * \param cell_idx Index of cell to activate
	   * \param v_act Initial value of potential perturbation 
	   *              (see \f$A\f$ in equation above)
	   * \param t0 Time after which the perturbation has an exponential decay
	   * \param alpha_rise Rise constant of the exponential 
	   *              (\f$\alpha > 0\f$ cf equation above)
	   * \param alpha_decay Decay constant of the exponential 
	   *              (\f$\alpha > 0\f$ cf equation above)
	   */	  
	  Potential_Activator_Tanh(size_type cell_idx,
				   electric_potential_ut v_act,
				   time_ut t0,
				   double alpha_rise,
				   double alpha_decay);
     };

     GCC_DIAG_ON(effc++)
     
} // namespace activator
} // namespace model_halidi
} // namespace cellsim

#endif //POTENTIAL_ACTIVATOR_TANH_HPP_INCLUDED

/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POTENTIAL_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED
#define POTENTIAL_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED

// #include "equation_type.hpp"
// #include "units.hpp"
// #include "cell_koenigsberger.hpp"

#include "koenigsberger_activator_base.hpp"
#include "potential_activator_base.hpp"

namespace cellsim {
namespace model_koenigsberger {
namespace activator {
     GCC_DIAG_OFF(effc++)

     /*!
      * \brief Cell activator for the membrane potential
      *
      * This class is designed to activate a cell using a perturbation of
      * its membrane potential. The perturbation follows an exponential decay
      * given by:
      * \f[ u(t < t_0) = A \text{ and } u(t > t_0) = A\,e^{-\alpha (t-t_0)} \f]
      *
      * \note This class' purpose is to activate a \e unique cell
      *       each time the activator is called (= each timestep)
      */
     class Potential_Activator : public Activator_Base,
				 public Potential_Activator_Base
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
	   * \param alpha Decay constant of the exponential 
	   *              (\f$\alpha > 0\f$ cf equation above)
	   */
	  Potential_Activator(size_type start_column,
			      size_type end_column,
			      electric_potential_ut v_act,
			      time_ut t0,
			      double alpha);
     };
     
     GCC_DIAG_ON(effc++)
} // namespace activator
} // namespace model_koenigsberger
} // namespace cellsim

#endif //POTENTIAL_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED

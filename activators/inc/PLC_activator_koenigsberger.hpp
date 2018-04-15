/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef PLC_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED
#define PLC_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED

#include <cmath>

#include "equation_type.hpp"
#include "koenigsberger_activator_base.hpp"
#include "units.hpp"

namespace cellsim {
namespace model_koenigsberger {
namespace activator {
     GCC_DIAG_OFF(effc++)
     /*!
      * Cell activator for the PLC agonist concentration in the cell
      *
      * This version is tailor-made to work with the Koenigsberger model.
      */
     class PLC_Activator : public Activator_Base
     {
     public:
	  /*!
	   * \brief This activator equation type.
	   *
	   * In other term, this is the only equation which this activator
	   * will have an effect on.
	   */
	  typedef equation_type::IP3_EQN equation_type;

	  /*!
	   * \brief Simple constructor
	   *
	   * The resulting PLC_Activator will activate columns in range
	   * [start, end) (ie. end non inclusive)
	   *
	   * \param start_column Index of first column to activate
	   * \param end_column Index of last column (non-inclusive) to activate
	   * \param J_act Activation PLC agonist concentration
	   */
	  PLC_Activator(size_type start_column,
			size_type end_column, 
			flux_ut J_act);

	  //! Simple accessor
	  flux_ut Jact() const {return J_act_;}
	  
	  /*!
	   * \brief Function operator
	   *
	   * This operator expects the following to be valid:
	   * \li f must accept exactly 2 parameters: f(c, a)
	   * \li f return type must be convertible to flux_ut
	   *
	   * \param f Function to apply
	   * \param c Current cell
	   * \param a Value of integration variable at current step
	   * \return Value of derivative with activation
	   */
	  template <typename function_t,
		    typename cell_t,
		    typename variable_t>
	  flux_ut operator()(const function_t f,
			     cell_t c,
			     const variable_t a) const
	       {
		    return f(c, a) + J_act_;
	       }
	  
     private:
	  /*!
	   * Activation value for the PLCago concentration
	   *
	   * This value will be added to the background PLCago concentration.
	   */
	  const flux_ut J_act_;
     };
     GCC_DIAG_ON(effc++)

} // namespace activator
} // namespace model_koenigsberger
} // namespace cellsim


#endif //PLC_ACTIVATOR_KOENIGSBERGER_HPP_INCLUDED

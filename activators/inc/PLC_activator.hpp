/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef PLC_ACTIVATOR_HPP_INCLUDED
#define PLC_ACTIVATOR_HPP_INCLUDED

#include "halidi_activator_base.hpp"
#include "equation_type.hpp"
#include "units.hpp"

#include <cmath>

namespace cellsim {
namespace model_halidi {
namespace activator {
     GCC_DIAG_OFF(effc++)
     /*!
      * Cell activator for the PLC agonist concentration in the cell
      *
      * This activation functor should be used only to activate a
      * particular cell with a greater PLC agonist concentration.
      */
     class PLC_Activator : public Activator_Base
     {
     public:
	  /*!
	   * This activator equation type.
	   *
	   * In other term, this is the only equation which this activator
	   * will have an effect on.
	   */
	  typedef equation_type::IP3_EQN equation_type;

	  /*!
	   * \brief Simple constructor
	   *
	   * \param cell_idx Index of cell to activate
	   * \param center_idx Index of the point at the center of the 
	   *                   activation region
	   * \param region_size Size of activation region in nb of points
	   * \param J_act Activation PLC agonist concentration
	   */
	  PLC_Activator(size_type cell_idx,
			size_type center_idx,
			size_type region_size,
			flux_ut J_act);
	  	  
	  /*!
	   * \brief Function operator
	   *
	   * This operator expects the following to be valid:
	   * \li f must accept exactly 2 parameters: f(i, a)
	   * \li f return type must be convertible to potential_SR_ut
	   *
	   * \param f Function to apply
	   * \param i Index of current point
	   * \param a Value of integration variable at current step
	   * \return Value of derivative with activation
	   */
	  template <typename function_t,
		    typename point_index_t,
		    typename variable_t>
	  flux_ut operator()(const function_t f,
			     const point_index_t i,
			     const variable_t a) const
	  {
	       size_type d(0);
	       if (i.index() > center_) {
		    d = i.index() - center_;
	       }
	       else if (i.index() < center_) {
		    d = center_ - i.index() + 1;
	       }

	       // only activate if in activation region
	       if (d > region_size2_) {
		    return f(i, a);
	       }
	       else {
		    return f(i, a) + J_act_;
	       }
	  }
	  
     private:
	  const flux_ut J_act_; //!< Activation value for the PLCago concentration
	  const size_type center_; //!< Center of activation region
	  const size_type region_size2_; //!< Half the size of the activation region
     };

     GCC_DIAG_ON(effc++)
} // namespace activator
} // namespace model_halidi
} // namespace cellsim

#endif //PLC_ACTIVATOR_HPP_INCLUDED

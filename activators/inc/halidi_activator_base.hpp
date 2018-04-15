/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef HALIDI_ACTIVATOR_BASE_HPP_INCLUDED
#define HALIDI_ACTIVATOR_BASE_HPP_INCLUDED

#include "definitions.hpp"

namespace cellsim {
namespace model_halidi {
namespace activator {
     GCC_DIAG_OFF(effc++)
     /*!
      * \brief Base class for activators within the Koenigsberger model
      */
     class Activator_Base
     {
     public:
	  /*! 
	   * \brief Simple constructor
	   *
	   * \param cell_idx Index of cell to activate
	   */
	  Activator_Base(size_type cell_idx);

	  /*!
	   * \brief Check method
	   *
	   * This method is used to determine whether or not to do the activation
	   * \param index Index of cell
	   * \return True if cell needs to be activated, false otherwise.
	   */
	  bool do_activation(size_type index) const;

     private:
	  const size_type cell_idx_; //!< Index of the cell to activate
     };
     GCC_DIAG_ON(effc++)
} // namespace activator
} // namespace model_halidi
} // namespace cellsim

#endif //HALIDI_ACTIVATOR_BASE_HPP_INCLUDED

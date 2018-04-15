/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef KOENIGSBERGER_ACTIVATOR_BASE_HPP_INCLUDED
#define KOENIGSBERGER_ACTIVATOR_BASE_HPP_INCLUDED

#include "definitions.hpp"

namespace cellsim {
namespace model_koenigsberger {
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
	   * \param start_column Index of first column to activate
	   * \param end_column Index after last column to activate
	   *                   (ie. activation until end_column - 1)
	   */
	  Activator_Base(size_type start_column,
			 size_type end_column);

	  /*!
	   * \brief Check method
	   *
	   * This method is used to determine whether or not to do the activation
	   * \param index Index of cell
	   * \return True if cell needs to be activated, false otherwise.
	   */
	  bool do_activation(size_type index) const;

     private:
	  const size_type start_column_; //!< Index of first column to activate
	  const size_type end_column_; //!< Index of last column to activate
     };
     GCC_DIAG_ON(effc++)
} // namespace activator
} // namespace model_koenigsberger
} // namespace cellsim

#endif //KOENIGSBERGER_ACTIVATOR_BASE_HPP_INCLUDED

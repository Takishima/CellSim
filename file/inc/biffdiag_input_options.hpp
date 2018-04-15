/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef BIFFDIAG_INPUT_OPTIONS_HPP_INCLUDED
#define BIFFDIAG_INPUT_OPTIONS_HPP_INCLUDED

#include "definitions.hpp"
#include "units.hpp"
#include "print_input_options.hpp"

namespace boost {
     namespace program_options {
	  class options_description;
	  class variables_map;
     }
}

namespace cellsim {
namespace input {
     namespace po = boost::program_options;
     /*!
      * \brief Definition of the options related to the bifurcation diagram
      *
      * \return Option description object with all the options that have to do
      *         with the \ref BIFFDIAG_HALIDI and  \ref BIFFDIAG_KOENIGS
      *         execution modes.
      */
     po::options_description biffdiag_input_options();

     GCC_DIAG_OFF(effc++)

     struct Biffdiag_Input {
	  size_type Jstep;
	  flux_ut Jmin;
	  flux_ut Jmax;	  
     };

     GCC_DIAG_ON(effc++)

     /*!
      * \brief Read bifurcation diagram option values
      *
      * \param vm_input Variables map to read data from
      * \param error Boolean value set to \c true if anything went wrong, 
      *              unmodified otherwise.
      * \return Biffdiag_Input with data
      */
     Biffdiag_Input read_biffdiag_values(po::variables_map vm_input, bool& error);

     /*!
      * \brief Print the values of the options to stdout
      *
      * \param value Input_Value to print
      */
     void print_option(Biffdiag_Input value);

} // namespace input
} // namespace cellsim


#endif //BIFFDIAG_INPUT_OPTIONS_HPP_INCLUDED

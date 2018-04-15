/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef HALIDI_INPUT_OPTIONS_HPP_INCLUDED
#define HALIDI_INPUT_OPTIONS_HPP_INCLUDED

#include "definitions.hpp"
#include "units_io.hpp"
#include "print_input_options.hpp"

#include <vector>

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
      * \brief Definition of all the halidi options
      *
      * \return Option description object with all the options that have to do
      *         with the \ref HALIDI execution mode
      */
     po::options_description halidi_input_options();

     GCC_DIAG_OFF(effc++)

     struct Halidi_Input {
	  // General options
     	  size_type npoints;
     	  size_type ncells;
     	  length_ut dx;
     	  flux_ut JPLCago;
	  
	  // Activation options
	  size_type cell_index;

#ifdef PLC_ACTIVATION
	  flux_ut Jact;
	  size_type center_index;
	  size_type region_size;
#elif defined POTASSIUM_CHORIDE_ACTIVATION
	  time_ut t0;
	  double v_rise;
	  double v_decay;
	  electric_potential_ut vK;
	  electric_potential_ut vCl;
#elif defined ELECTRIC_ACTIVATION
	  // If using electrical potential activation
	  electric_potential_ut vact;
	  time_ut t0;
	  double v_rise;
	  double v_decay;
#endif // PLC_ACTIVATION

	  // Output options
	  size_type point_index;
 	  std::vector<size_type> cell_indices;
    };

     GCC_DIAG_ON(effc++)

     /*!
      * \brief Read Halidi option values
      *
      * \param vm_input Variables map to read data from
      * \param error Boolean value set to \c true if anything went wrong, 
      *              unmodified otherwise.
      * \return Halidi_Input with data
      */
     Halidi_Input read_halidi_values(po::variables_map vm_input, bool& error);

     /*!
      * \brief Print the values of the options to stdout
      *
      * \param value Input_Value to print
      */
     void print_option(Halidi_Input value);

} // namespace input
} // namespace cellsim


#endif //HALIDI_INPUT_OPTIONS_HPP_INCLUDED

/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef INPUT_OPTIONS_HPP_INCLUDED
#define INPUT_OPTIONS_HPP_INCLUDED

#include "definitions.hpp"
#include "units_io.hpp"
#include "read_parameters.hpp"
#include "print_input_options.hpp"

#include <string>

namespace boost {
     namespace program_options {
	  class options_description;
	  class variables_map;
     }
}

namespace cellsim {
//! All that is related to program input
namespace input {
     namespace po = boost::program_options;

     GCC_DIAG_OFF(effc++)

     /*!
      * \brief Basic tuple-like class
      *
      * Only useful to return all the generic parameters
      */
     struct Input_Value
     {
	  length_ut cell_length;
	  concentration_ut c_init;
	  concentration_ut I_init;
	  concentration_ut s_init;
	  electric_potential_ut v_init;
	  double w_init;
     };

     GCC_DIAG_ON(effc++)

     /*!
      * \brief Definition of all generic options
      *
      * \return Option description object with all the options common
      *         to all the execution mode.
      */
     po::options_description generic_input_options();
     /*!
      * \brief Read generic option values
      *
      * \param vm_input Variables map to read data from
      * \param error Boolean value set to \c true if anything went wrong, 
      *              unmodified otherwise.
      * \return Input_Value with data
      */
     Input_Value read_generic_values(po::variables_map vm_input, bool& error);

     /*!
      * \brief Print the values of the options to stdout
      *
      * \param value Input_Value to print
      */
     void print_option(Input_Value value);

     
     GCC_DIAG_OFF(effc++)

     struct Simulation_Input {
	  double t_equ;
	  time_ut dt_equ;
	  size_type output_step_equ;

	  double t_act;
	  time_ut dt_act;
	  size_type output_step_act;

	  double t_run;
	  time_ut dt_run;
	  size_type output_step_run;
     };

     GCC_DIAG_ON(effc++)

     /*!
      * \brief Definition of all simulation options
      *
      * \return Option description object with all the options relating
      *         to the simulation
      */
     po::options_description simulation_input_options();

     /*!
      * \brief Read simulation option values
      *
      * \param vm_input Variables map to read data from
      * \param error Boolean value set to \c true if anything went wrong, 
      *              unmodified otherwise.
      * \return Simulation_Input with data
      */
     Simulation_Input read_simulation_values(po::variables_map vm_input, bool& error);

     /*!
      * \brief Print the values of the options to stdout
      *
      * \param value Input_Value to print
      */
     void print_option(Simulation_Input value);

}
} // namespace cellsim

#endif //INPUT_OPTIONS_HPP_INCLUDED

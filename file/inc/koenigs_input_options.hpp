/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef KOENIGS_INPUT_OPTIONS_HPP_INCLUDED
#define KOENIGS_INPUT_OPTIONS_HPP_INCLUDED

#include "definitions.hpp"
#include "units_io.hpp"
#include "print_input_options.hpp"

#include <tuple>
#include <vector>

namespace boost {
     namespace program_options {
	  class options_description;
	  class variables_map;
     }
}

namespace cellsim {
namespace input {
     typedef std::tuple<cellsim::size_type, cellsim::size_type> Point;

     namespace po = boost::program_options;

     /*!
      * \brief Definition of all the koenigsberger options
      *
      * \return Option description object with all the options that have to do
      *         with the \ref KOENIGSBERGER execution mode
      */
     po::options_description koenigsberger_input_options();

     GCC_DIAG_OFF(effc++)

     struct Koenigs_Input {
     	  size_type Nx;
     	  size_type Ny;
     	  flux_ut JPLCago;

     	  size_type start_col;
     	  size_type end_col;

#ifdef PLC_ACTIVATION
     	  flux_ut Jact;
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
#else
#  error Unknown or unspecified activation type for Koenigsberger! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION

	  // Output options
 	  std::vector<Point> cell_indices;
     };

     GCC_DIAG_ON(effc++)

     Koenigs_Input read_koenigs_values(po::variables_map vm_input, bool& error);

     /*!
      * \brief Print the values of the options to stdout
      *
      * \param value Koenigs_Input value to print
      */
     void print_option(Koenigs_Input value);

     /*!
      * \brief Definition of the options related to the Koenigsberger
      *        bifurcation diagram
      *
      * \return Option description object with all the options that have to do
      *         with the \ref BIFFDIAG_KOENIGS execution mode
      */
     po::options_description biffdiag_koenigs_input_options();

} // namespace input
} // namespace cellsim

#endif //KOENIGS_INPUT_OPTIONS_HPP_INCLUDED

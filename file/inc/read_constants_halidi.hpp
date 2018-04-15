/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef READ_CONSTANTS_HALIDI_HPP_INCLUDED
#define READ_CONSTANTS_HALIDI_HPP_INCLUDED

#include "definitions.hpp"
#include "read_constants.hpp"
 
#include <string>

namespace boost {
     namespace program_options {
	  class variables_map;
     }
}

namespace cellsim {
     namespace model_halidi {
	  struct Constants_Values;
     } // namespace model_halidi

namespace input {
     namespace po = boost::program_options;

     /*!
      * \brief Print all the constants value
      * \param C Model constants values
      */
     void print_constants_values(const model_halidi::Constants_Values& C);

     /*!
      * \brief Set all the constants value based on data from \c vm_constants
      *
      * Halidi model overload.
      *
      * \param vm_constants Input data
      * \param model_constants Objects with constants values to modify
      * \return True or false depending if successful or not.
      *         If unsuccessful for any reason (e.g. if one is missing),
      *         function is no-op on current constants values.
      */
     bool set_constants_values(po::variables_map vm_constants,
			       model_halidi::Constants_Values& model_constants);

     /*!
      * \brief Read model constants values from a file
      *
      * Halidi model overload.
      *
      * \param filename File to read data from
      * \param model_constants Objects with constants values to modify
      * \return Constants value read from file.
      *         If unsuccessful for any reason, all constants values are set
      *         to 0.
      */
     bool read_constants(std::string filename, 
			 model_halidi::Constants_Values& model_constants);

} // namespace input
} // namespace cellsim

#endif //READ_CONSTANTS_HALIDI_HPP_INCLUDED

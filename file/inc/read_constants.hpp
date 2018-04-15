/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef READ_CONSTANTS_HPP_INCLUDED
#define READ_CONSTANTS_HPP_INCLUDED

#include "definitions.hpp"

#include <string>
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

     template<typename T, size_type sizeOfArray> 
     size_type get_array_size (T (&)[sizeOfArray])
     {
	  return sizeOfArray;
     }
     
     /*!
      * \brief Helper structure for constants and their units
      *
      */
     struct constant_value {
	  const char* name; //!< Self-explanatory
	  const char* units; //!< Self-explanatory
     };

     /*!
      * \brief Definition of all the constants options
      *
      * \param list Array with all the constants names and units.
      * \param size Size of \c list
      *
      * \return Option description object with all the options that have to do
      *         with model constants
      */
     po::options_description constants_input_options(
	  const constant_value* list, size_type size);

     /*!
      * Read the filename of the constants file from input data
      *
      * \param vm_input Input data
      * \return Filename
      */
     std::string get_constants_filename(po::variables_map vm_input);

} // namespace input
} // namespace cellsim

#endif //READ_CONSTANTS_HPP_INCLUDED

/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef ARGUMENT_PARSING_HPP_INCLUDED
#define ARGUMENT_PARSING_HPP_INCLUDED

#include "definitions.hpp"

#include <string>
#include <tuple>

namespace cellsim {
namespace input {
     /*!
      * \brief Parse command line arguments
      *
      * \param argc Number of arguments (including program name)
      * \param argv Array of command line arguments
      * \return Tuple with \c <input_file,\c output_file>
      * \note If \c input_file, can also be one of the following:
      *       \li \c 1 Meaning an error occurred
      *       \li \c help-more Meaning we need to print the extended help message\n
      *
      * In those cases, \c output_file is always the empty string.
      */
     std::tuple<std::string,std::string> parse_arguments(int argc, char** argv);
}
} // namespac cellsim

#endif //ARGUMENT_PARSING_HPP_INCLUDED

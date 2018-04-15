/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef READ_PARAMETERS_HPP_INCLUDED
#define READ_PARAMETERS_HPP_INCLUDED

#include "definitions.hpp"

#include <string>
#include <vector>

namespace  cellsim {
namespace  input {
     // namespace po = boost::program_options;

     /*!
      * \brief Read an \ref EXEC_MODE value from a string
      *
      * \param str Input data
      * \return Data read and converted from input
      */
     EXEC_MODE read_exec_mode(std::string str);
    
} // namespace input
} // namespace cellsim

#endif //READ_PARAMETERS_HPP_INCLUDED

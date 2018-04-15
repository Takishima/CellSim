/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef PRINT_INPUT_OPTIONS_HPP_INCLUDED
#define PRINT_INPUT_OPTIONS_HPP_INCLUDED

#include <iostream>

namespace cellsim {
     template <typename T>
     void print_input_options_impl(T t)
     {
	  print_option(t);
     }

     template <typename T, typename... Args>
     void print_input_options_impl(T t, Args... args)
     {
	  print_option(t);
	  print_input_options_impl(args...);
     }

     template <typename... Args>
     void print_input_options(Args... args)
     {
	  std::cout << "========================================"
		    << "========================================\n"
		    << "INPUT FILE:\n";
	  print_input_options_impl(args...);
	  std::cout << "########################################"
		    << "########################################\n";
		   
     }
} // namespace cellsim

#endif //PRINT_INPUT_OPTIONS_HPP_INCLUDED

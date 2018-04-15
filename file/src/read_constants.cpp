#include "read_constants.hpp"
#include "constants_halidi.hpp"

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>

namespace po = boost::program_options;

po::options_description cellsim::input::constants_input_options(
     const constant_value* list, 
     size_type size)
{
     using boost::program_options::options_description;
     using boost::program_options::value;
     
     options_description constants_options("Model constants and their units");

     // Add all constants options definitions from data in list
     std::for_each(list,
		   list + size,
		   [&constants_options] (const constant_value& el) {
			constants_options.add_options()
			     (el.name, value<double>(), 
			      std::string(std::string("Value for ") + el.name 
					  + " " + el.units).c_str()
				  );
		   }
	  );

     return constants_options;
}

std::string cellsim::input::get_constants_filename(po::variables_map vm_input)
{
     if (vm_input["system.constants_file"].empty()) {
	  std::cerr << "Missing value : system.constants_file\n";
	  return "";
     }
     else {
	  return vm_input["system.constants_file"].as<std::string>();
     }
}

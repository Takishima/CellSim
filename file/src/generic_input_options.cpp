#include "input_options.hpp"
#include <iostream>
#include <boost/program_options.hpp>

using namespace boost::program_options;

options_description cellsim::input::generic_input_options()
{
     options_description generic_options;

     // Model parameters and initial values
     generic_options.add_options()
	  ("system.constants_file", value<std::string>()->required(),
	   "Path to constants file (relative path).")
	  ("system.cell_length", value<double>()->required(),
	   "Length of cell (in [uM])")
	  ("system.c", value<double>()->required(),
	   "Initial value for [Ca2+] (in [uM])")
	  ("system.I", value<double>()->required(),
	   "Initial value for [IP3] (in [uM])")
	  ("system.s", value<double>()->required(),
	   "Initial value for [Ca2+] in SR (in [uM])")
	  ("system.v", value<double>()->required(),
	   "Initial value for the membrane potential (in [mV])")
	  ("system.w", value<double>()->required(),
	   "Initial value for the open state probability (dimensionless)")
	  ;

     return generic_options;
}

cellsim::input::Input_Value 
cellsim::input::read_generic_values(variables_map vm_input, bool& error)
{
     if (vm_input["system.cell_length"].empty()) {
	  std::cerr << "Missing value : system.cell_length\n";
	  error = true;
     }
     if (vm_input["system.c"].empty()) {
	  std::cerr << "Missing value : system.c\n";
	  error = true;
     }
     if (vm_input["system.I"].empty()) {
	  std::cerr << "Missing value : system.I\n";
	  error = true;
     }
     if (vm_input["system.s"].empty()) {
	  std::cerr << "Missing value : system.s\n";
	  error = true;
     }
     if (vm_input["system.v"].empty()) {
	  std::cerr << "Missing value : system.v\n";
	  error = true;
     }
     if (vm_input["system.w"].empty()) {
	  std::cerr << "Missing value : system.w\n";
	  error = true;
     }
     
     namespace po = boost::program_options;

     if (!error) {
	  return Input_Value{
	       length_ut::from_value(vm_input["system.cell_length"].as<double>()),
		    concentration_ut::from_value(vm_input["system.c"].as<double>()),
		    concentration_ut::from_value(vm_input["system.I"].as<double>()),
		    concentration_ut::from_value(vm_input["system.s"].as<double>()),
		    electric_potential_ut::from_value(vm_input["system.v"].as<double>()),
		    vm_input["system.w"].as<double>()
		    };
     }
     else {
	  return Input_Value();
     }
}

void cellsim::input::print_option(Input_Value value)
{
     std::cout << "System:\n"
	       << "    cell_length = " << value.cell_length << std::endl
	       << "    c =           " << value.c_init << std::endl
	       << "    I =           " << value.I_init << std::endl
	       << "    s =           " << value.s_init << std::endl
	       << "    v =           " << value.v_init << std::endl
	       << "    w =           " << value.w_init << std::endl;
}

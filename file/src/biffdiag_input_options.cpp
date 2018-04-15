#include "biffdiag_input_options.hpp"

#include <iostream>
#include <stdexcept>
#include <boost/program_options.hpp>

namespace input = cellsim::input;

using namespace boost::program_options;


options_description cellsim::input::biffdiag_input_options()
{
     options_description biffdiag_options;
     biffdiag_options.add_options()
	  ("biffdiag.Jstep", value<size_type>(),
	   "Number of values of JPLCago to simulate.")
	  ("biffdiag.Jmin", value<double>(),
	   "Minimum background PLC agonist current.")
	  ("biffdiag.Jmax", value<double>(),
	   "Maximum background PLC agonist current.")
	  ;

     return biffdiag_options;
}

cellsim::input::Biffdiag_Input
cellsim::input::read_biffdiag_values(po::variables_map vm_input, bool& error)
{
     if (vm_input["biffdiag.Jstep"].empty()) {
	  std::cerr << "Missing value : biffdiag.Jstep\n";
	  error = true;
     }
     if (vm_input["biffdiag.Jmin"].empty()) {
	  std::cerr << "Missing value : biffdiag.Jmin\n";
	  error = true;
     }
     if (vm_input["biffdiag.Jmax"].empty()) {
	  std::cerr << "Missing value : biffdiag.Jmax\n";
	  error = true;
     }

     if (!error) {
	  return Biffdiag_Input{
	       vm_input["biffdiag.Jstep"].as<size_type>(),
		    flux_ut::from_value(vm_input["biffdiag.Jmin"].as<double>()),
		    flux_ut::from_value(vm_input["biffdiag.Jmax"].as<double>())
		    };
     }
     else {
	  return Biffdiag_Input();
     }

}

void cellsim::input::print_option(Biffdiag_Input value)
{
     std::cout << "Bifurcation Diagram Halidi:\n"
	       << "    Jstep =        " << value.Jstep << std::endl
	       << "    Jmin =         " << value.Jmin << std::endl
	       << "    Jmax =         " << value.Jmax << std::endl;
}


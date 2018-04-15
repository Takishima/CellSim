#include "input_options.hpp"

#include <iostream>
#include <boost/program_options.hpp>

using namespace boost::program_options;

options_description cellsim::input::simulation_input_options()
{
     options_description simulation_options;
     simulation_options.add_options()
	  ("simulation.dt", value<double>(),
	   "Integration time step (for all parts) (in [s])")
	  ("simulation.output_step", value<size_type>()->default_value(10),
	   "Output data every 'arg' time steps (for all parts).")

	  // For equlibration
	  ("simulation.t_equ", value<double>()->default_value(1),
	   "Time in [s] to equilibrate the system")
	  ("simulation.dt_equ", value<double>(),
	   "Integration time step (for equilibration, in [s])")
	  ("simulation.output_step_equ", value<size_type>(),
	   "Output data every 'arg' time steps (for equilibration).")

	  // For activation
	  ("simulation.t_act", value<double>()->default_value(1),
	   "Time in [s] to activate the system.")
	  ("simulation.dt_act", value<double>(),
	   "Integration time step (for activation, in [s])")
	  ("simulation.output_step_act", value<size_type>(),
	   "Output data every 'arg' time steps (for activation).")

	  // For run after activation
	  ("simulation.t_run", value<double>()->default_value(10),
	   "Time in [s] to run the simulation after activation.")
	  ("simulation.dt_run", value<double>(),
	   "Integration time step (for run, in [s])")
	  ("simulation.output_step_run", value<size_type>(),
	   "Output data every 'arg' time steps (for run).")
	  ;

     return simulation_options;
}


cellsim::input::Simulation_Input
cellsim::input::read_simulation_values(po::variables_map vm_input, bool& error)
{
     // ======================
     // Take care of timesteps

     if (vm_input["simulation.dt"].empty() &&
	 vm_input["simulation.dt_equ"].empty() && 
	 vm_input["simulation.dt_act"].empty() && 
	 vm_input["simulation.dt_run"].empty()) {
	  std::cerr << "Missing value : simulation.dt (and/or)"
		    << " simulation.dt_[equ,act,run]\n";
	  error = true;
     }
     else if (vm_input["simulation.dt"].empty() &&
	 (vm_input["simulation.dt_equ"].empty() ||
	  vm_input["simulation.dt_act"].empty() ||
	  vm_input["simulation.dt_run"].empty())) {
	  std::cerr << "Missing value : got only part of simulation.dt_[equ,act,run]"
		    << " and no simulation.dt\n";
	  error = true;
     }

     // ============

     double dt_equ(0.), dt_act(0.), dt_run(0.);
     if (!vm_input["simulation.dt"].empty()) {
	  dt_equ = vm_input["simulation.dt"].as<double>();
	  dt_act = dt_equ;
	  dt_run = dt_equ;
     }
     if (!vm_input["simulation.dt_equ"].empty()) {
	  dt_equ = vm_input["simulation.dt_equ"].as<double>();
     }
     if (!vm_input["simulation.dt_act"].empty()) {
	  dt_act = vm_input["simulation.dt_act"].as<double>();
     }
     if (!vm_input["simulation.dt_run"].empty()) {
	  dt_run = vm_input["simulation.dt_run"].as<double>();
     }
     // ======================
     // Take care of output steps

     if (vm_input["simulation.output_step"].empty() &&
	 vm_input["simulation.output_step_equ"].empty() && 
	 vm_input["simulation.output_step_act"].empty() && 
	 vm_input["simulation.output_step_run"].empty()) {
	  std::cerr << "Missing value : simulation.output_step (and/or)"
		    << " simulation.output_step_[equ,act,run]\n";
	  error = true;
     }
     else if (vm_input["simulation.output_step"].empty() &&
	 (vm_input["simulation.output_step_equ"].empty() ||
	  vm_input["simulation.output_step_act"].empty() ||
	  vm_input["simulation.output_step_run"].empty())) {
	  std::cerr << "Missing value : got only part of simulation.output_step_[equ,act,run]"
		    << " and no simulation.output_step\n";
	  error = true;
     }

     // ============

     size_type output_step_equ(0), output_step_act(0), output_step_run(0);
     if (!vm_input["simulation.output_step"].empty()) {
	  output_step_equ = vm_input["simulation.output_step"].as<size_type>();
	  output_step_act = output_step_equ;
	  output_step_run = output_step_equ;
     }
     if (!vm_input["simulation.output_step_equ"].empty()) {
	  output_step_equ = vm_input["simulation.output_step_equ"].as<size_type>();
     }
     if (!vm_input["simulation.output_step_act"].empty()) {
	  output_step_act = vm_input["simulation.output_step_act"].as<size_type>();
     }
     if (!vm_input["simulation.output_step_run"].empty()) {
	  output_step_run = vm_input["simulation.output_step_run"].as<size_type>();
     }
     // ======================
     
     if (vm_input["simulation.t_equ"].empty()) {
	  std::cerr << "Missing value : simulation.t_equ\n";
	  error = true;
     }
     if (vm_input["simulation.t_act"].empty()) {
	  std::cerr << "Missing value : simulation.t_act\n";
	  error = true;
     }
     if (vm_input["simulation.t_run"].empty()) {
	  std::cerr << "Missing value : simulation.t_run\n";
	  error = true;
     }

     if (!error) {
	  return Simulation_Input{
	       vm_input["simulation.t_equ"].as<double>(),
		    time_ut::from_value(dt_equ),
		    output_step_equ,

		    vm_input["simulation.t_act"].as<double>(),
		    time_ut::from_value(dt_act),
		    output_step_act,
		    
		    vm_input["simulation.t_run"].as<double>(),
		    time_ut::from_value(dt_run),
		    output_step_run,
		    };
     }
     else {
	  return Simulation_Input();
     }
}

void cellsim::input::print_option(Simulation_Input value)
{
     std::cout << "Simulation:\n"
	       << "  Equblibration:\n"
	       << "    t_equ =           " << time_ut::from_value(value.t_equ) << std::endl
	       << "    dt_equ =          " << value.dt_equ << std::endl
	       << "    nstep_equ =       " << (value.t_equ/value.dt_equ.value()) << std::endl
	       << "    output_step_equ = " << value.output_step_equ << std::endl
	       << "  Activation:\n"
	       << "    t_act =           " << time_ut::from_value(value.t_act) << std::endl
	       << "    dt_act =          " << value.dt_act << std::endl
	       << "    nstep_act =       " << (value.t_act/value.dt_act.value()) << std::endl
	       << "    output_step_act = " << value.output_step_act << std::endl
	       << "  Run:\n"
	       << "    t_run =           " << time_ut::from_value(value.t_run) << std::endl
	       << "    dt_run =          " << value.dt_run << std::endl
	       << "    nstep_run =       " << (value.t_run/value.dt_run.value()) << std::endl
	       << "    output_step_run = " << value.output_step_run << std::endl;
}

#include "halidi_bifurcation_diagram.hpp"

#include "definitions.hpp"
#include "utility.hpp"
#include "input_options.hpp"
#include "halidi_input_options.hpp"
#include "constants_halidi.hpp"
#include "biffdiag_input_options.hpp"
#include "cell.hpp"
#include "system.hpp"

#include <iostream>
#include <fstream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

int cellsim::run_halidi_biffdiag(po::variables_map vm_input,
				 std::string output_file,
				 const halidi::Constants_Values& model_constants)
{
     using cellsim::size_type;
     using namespace cellsim::model_halidi;

     bool error(false);

     // ==========

     input::Input_Value i_v = input::read_generic_values(vm_input, error);
     
     constexpr auto npoints(1), ncells(1);

     // ==========

     input::Halidi_Input h_i = input::read_halidi_values(vm_input, error);
     input::Biffdiag_Input bh_i = input::read_biffdiag_values(vm_input, error);
     input::Simulation_Input s_i = input::read_simulation_values(vm_input, error);

     const auto& dt_equ(s_i.dt_equ);
     const auto& dt_run(s_i.dt_run);

     const size_type nstep_equ(s_i.t_equ / dt_equ.value());
     const size_type nstep_run(s_i.t_run / dt_run.value());

     // ==========

     if (error) {
	  std::cerr << "Error while parsing values from config file!\n";
	  return 1;
     }

     // ==========
     
     std::cout << "Output file pattern: " << output_file << std::endl;

     print_input_options(i_v, s_i, h_i, bh_i);
     // input::print_generic_option_values(i_v);
     // input::print_simulation_option_values(s_i);
     std::cout << "** NB: ignoring simulation.t_act & simulation.dt_act        **\n";
     // input::print_halidi_option_values(h_i);
     std::cout << "** NB: ignoring halidi.npoints,halidi.ncells,halidi.JPLCago **\n"
	       << "**              & halidi.activation.*                       **\n";
     // input::print_biffdiag_option_values(bh_i);
     std::cout << "\n****************************************"
	       << "****************************************\n\n";

     // ==========

     const auto split = utility::split_filename(output_file);

     std::ofstream cfile(utility::create_filename(split, "", "ca_out"));
     std::ofstream Ifile(utility::create_filename(split, "", "ip3_out"));

     // ==========
     
     Cell single_cell(npoints, 
		      i_v.c_init, i_v.I_init, i_v.s_init, i_v.v_init, i_v.w_init,
		      model_constants); 

     System mysystem(ncells, single_cell, h_i.dx, false);

     std::cout << "Start finding values!" << std::endl;

     using namespace boost::accumulators;
     typedef accumulator_set<double, stats<tag::min, tag::max, tag::mean>> acc_type;

     const auto dJ((bh_i.Jmax - bh_i.Jmin) / (bh_i.Jstep * 1.));

     for(auto JPLCago(bh_i.Jmin) ; JPLCago <= bh_i.Jmax ; JPLCago += dJ){
     	  std::cout << "\tcalculating for JPLCago = " << JPLCago << std::endl;
	  
     	  for(size_type i(0) ; i < nstep_equ ; ++i){
     	       mysystem.do_step(dt_equ, JPLCago);
     	  }
	  acc_type ca_acc, ip3_acc;

     	  for(size_type i(0); i < nstep_run; ++i){
     	       mysystem.do_step(dt_run, JPLCago);
	       ca_acc(mysystem.get_c_at(0).value());
	       ip3_acc(mysystem.get_I_at(0).value());
     	  }
	  cfile << JPLCago << " "
		<< mean(ca_acc) << " "
		<< min(ca_acc) << " "
		<< max(ca_acc) << std::endl;
	  Ifile << JPLCago << " "
		<< mean(ip3_acc) << " "
		<< min(ip3_acc) << " "
		<< max(ip3_acc) << std::endl;
     }

     return 0;
}

#pragma GCC diagnostic pop

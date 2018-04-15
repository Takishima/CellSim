#include "halidi.hpp"

#include "utility.hpp"
#include "input_options.hpp"
#include "halidi_input_options.hpp"
#include "system.hpp"
#include "cell.hpp"
#include "PLC_activator.hpp"
#include "potassium_chloride_activator.hpp"
#include "potential_activator.hpp"
#include "potential_activator_tanh.hpp"
#include "definitions.hpp"
#include "units_io.hpp"
#include "output_streams.hpp"
#include "flux_streams.hpp"
#include "read_constants_halidi.hpp"

#include "progress_bar.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-equal"

int cellsim::run_halidi(po::variables_map vm_input, 
			std::string output_file,
			const halidi::Constants_Values& model_constants)
{
     using cellsim::size_type;
     using namespace cellsim;
     using namespace cellsim::model_halidi;

     bool error(false);

     // ==========

     input::Input_Value init_value = input::read_generic_values(vm_input, error);
     
     const auto& c(init_value.c_init);
     const auto& I(init_value.I_init);
     const auto& s(init_value.s_init);
     const auto& v(init_value.v_init);
     const auto& w(init_value.w_init);

     // ==========

     input::Halidi_Input h_i = input::read_halidi_values(vm_input, error);
     input::Simulation_Input s_i = input::read_simulation_values(vm_input, error);

     const auto& dt_equ(s_i.dt_equ);
     const auto& dt_act(s_i.dt_act);
     const auto& dt_run(s_i.dt_run);

     size_type nstep_equ(s_i.t_equ / dt_equ.value());
     size_type nstep_act(s_i.t_act / dt_act.value());
     size_type nstep_run(s_i.t_run / dt_run.value());

     // ==========

     if (error) {
	  std::cerr << "Error while parsing values from config file!\n";
	  return 1;
     }

     // ==========

     std::cout << "Output file pattern: " << output_file << std::endl << std::endl;

     print_input_options(init_value, s_i, h_i);
     input::print_constants_values(model_constants);
     std::cout << std::endl;

     // ================================

     Cell single_cell(h_i.npoints, c, I, s, v, w, model_constants);
     // System mysystem(h_i.Nx, h_i.Ny, c, I, s, v, w);
     System mysystem(h_i.ncells, single_cell, h_i.dx, false);

#ifdef PLC_ACTIVATION
     if (h_i.cell_index == 0 && 
	 h_i.center_index == 0 &&
	 h_i.region_size == 0 &&
	 h_i.Jact == flux_ut::from_value(0.)) {
	  std::cerr << "**********\n"
		    << "WARNING: using PLC activation and all "
		    << "related parameters are null!\n"
		    << "**********\n";
     }
     model_halidi::activator::PLC_Activator act(h_i.cell_index,
						h_i.center_index,
						h_i.region_size,
						h_i.Jact);
#elif defined ELECTRIC_ACTIVATION
     if (h_i.t0 == 0_s && h_i.vact == 0_mV && h_i.v_decay == 0.0) {
	  std::cerr << "**********\n"
		    << "WARNING: using electrical activation and all "
		    << "related parameters are null!\n"
		    << "**********\n";
     }
     model_halidi::activator::Potential_Activator_Tanh act(h_i.cell_index,
							   h_i.vact,
							   h_i.t0,
							   h_i.v_rise,
							   h_i.v_decay);
#elif defined POTASSIUM_CHORIDE_ACTIVATION
     if (h_i.cell_index == 0 && 
	 h_i.t0 == 0_s &&
	 h_i.v_rise == 0.0 &&
	 h_i.v_decay == 0.0 &&
	 h_i.vK == 0_mV &&
	 h_i.vCl == 0_mV) {
     	  std::cerr << "**********\n"
     		    << "WARNING: using potassium chloride activation and all "
     		    << "related parameters are null!\n"
     		    << "**********\n";
     }
     model_halidi::activator::K_Cl_Activator act(h_i.cell_index,
						 h_i.t0,
						 h_i.v_rise,
						 h_i.v_decay,
						 model_constants,
						 h_i.vK,
						 h_i.vCl);
#else
#  error Unknown or unspecified activation type for Halidi! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION

     // ================================

     // Output for the whole system (removed because takes too much space)
     output::Output_Streams system_output(output_file, "system");
     system_output << System::get_system_header();

     // ==========
     // Output selected cells
     output::Output_Streams cells_output(output_file);
     cells_output << "# t value for idx= ";
     for(auto i: h_i.cell_indices) {
	  cells_output << i << " ";
     }
     cells_output << std::endl;
     cells_output.precision(12);

     // ==========
     // Output for activation
     namespace util = cellsim::utility;
     std::ofstream act_out(util::create_filename(util::split_filename(output_file),
						 "activation", ""));

     // ==========

     output::Flux_Streams flux_out(h_i.cell_indices.size(), output_file);
     flux_out << Cell::get_fluxes_header();
     flux_out.precision(12);

     // ================================

     // mysystem.print_system(system_output);
     mysystem.print_cell_at(h_i.cell_indices, h_i.point_index, cells_output);
     mysystem.print_fluxes_at(h_i.cell_indices, h_i.point_index, flux_out);

     Progress_Bar p(nstep_equ, "Equilibration:   ");
     p.start();
     for(size_type i(0) ; i < nstep_equ ; ++i, ++p) {
     	  mysystem.do_step(dt_equ, h_i.JPLCago);
	  if (i % s_i.output_step_equ == 0) {
	       // mysystem.print_system(system_output);
	       mysystem.print_cell_at(h_i.cell_indices, h_i.point_index, cells_output);
	       mysystem.print_fluxes_at(h_i.cell_indices, h_i.point_index, flux_out);
	  }
     }
     mysystem.print_cell_at(h_i.cell_indices, h_i.point_index, cells_output);
     mysystem.print_fluxes_at(h_i.cell_indices, h_i.point_index, flux_out);
     mysystem.print_system(system_output);

     // ==========

     p.reset(nstep_act, "Activation:      ");
     p.start();
     for(size_type i(0) ; i < nstep_act ; ++i, ++p) {
     	  mysystem.activate(dt_act, h_i.JPLCago, act);
     	  if (i % s_i.output_step_act == 0) {
	       mysystem.print_system(system_output);
	       mysystem.print_cell_at(h_i.cell_indices, h_i.point_index, cells_output);
	       mysystem.print_fluxes_at(h_i.cell_indices, h_i.point_index, flux_out);
	       // act_out << i*dt_act.value() << " " << act.value().value()*dt_act.value() << std::endl;
     	  }
     }
     mysystem.print_system(system_output);
     mysystem.print_cell_at(h_i.cell_indices, h_i.point_index, cells_output);
     mysystem.print_fluxes_at(h_i.cell_indices, h_i.point_index, flux_out);

     // ==========

     p.reset(nstep_run, "Simulation:      ");
     p.start();
     for(size_type i(0) ; i < nstep_run ; ++i, ++p) {
     	  mysystem.do_step(dt_run, h_i.JPLCago);
     	  if (i % s_i.output_step_run == 0) {
	       // mysystem.print_system(system_output);
	       mysystem.print_cell_at(h_i.cell_indices, h_i.point_index, cells_output);
	       mysystem.print_fluxes_at(h_i.cell_indices, h_i.point_index, flux_out);
     	  }
     }

     // END of simulation
     // ========================================================================

     return 0;
}

#pragma GCC diagnostic pop

#include "koenigsberger.hpp"

#include "utility.hpp"
#include "input_options.hpp"
#include "koenigs_input_options.hpp"
#include "system_koenigsberger.hpp"
#include "cell_koenigsberger.hpp"
#include "PLC_activator_koenigsberger.hpp"
#include "potential_activator_koenigsberger.hpp"
#include "potential_activator_koenigsberger_tanh.hpp"
#include "potassium_chloride_activator_koenigsberger.hpp"
#include "definitions.hpp"
#include "units_io.hpp"
#include "output_streams.hpp"
#include "flux_streams.hpp"
#include "read_constants_koenigsberger.hpp"

#include "progress_bar.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-equal"

int cellsim::run_koenigsberger(po::variables_map vm_input, 
			       std::string output_file,
			       const koenigs::Constants_Values& model_constants)
{
     using cellsim::size_type;
     using namespace cellsim;
     using namespace cellsim::model_koenigsberger;

     bool error(false);

     // ==========

     input::Input_Value init_value = input::read_generic_values(vm_input, error);
     
     const auto& c(init_value.c_init);
     const auto& I(init_value.I_init);
     const auto& s(init_value.s_init);
     const auto& v(init_value.v_init);
     const auto& w(init_value.w_init);

     // ==========

     input::Koenigs_Input k_i = input::read_koenigs_values(vm_input, error);
     input::Simulation_Input s_i = input::read_simulation_values(vm_input, error);

     const auto& dt_equ(s_i.dt_equ);
     const auto& dt_run(s_i.dt_act);
     const auto& dt_act(s_i.dt_run);

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

     print_input_options(init_value, s_i, k_i);
     input::print_constants_values(model_constants);
     std::cout << std::endl;

     // ================================

     System mysystem(k_i.Nx, k_i.Ny, c, I, s, v, w, model_constants);

     // ==========

#ifdef PLC_ACTIVATION
     if (k_i.start_col == 0 && 
	 k_i.end_col == 0 &&
	 k_i.Jact == flux_ut::from_value(0.)) {
	  std::cerr << "**********\n"
		    << "WARNING: using PLC activation and all "
		    << "related parameters are null!\n"
		    << "**********\n";
     }
     activator::PLC_Activator act(k_i.start_col, k_i.end_col,
					(k_i.Jact - k_i.JPLCago));
#elif defined ELECTRIC_ACTIVATION
     if (k_i.start_col == 0 && 
	 k_i.end_col == 0 &&
	 k_i.vact == 0_mV &&
	 k_i.v_decay == 0.0) {
     	  std::cerr << "**********\n"
     		    << "WARNING: using electrical activation and all "
     		    << "related parameters are null!\n"
     		    << "**********\n";
     }
     model_koenigsberger::activator::Potential_Activator_Tanh act(k_i.start_col,
						   k_i.end_col,
						   k_i.vact,
						   k_i.t0,
						   k_i.v_rise,
						   k_i.v_decay);
#elif defined POTASSIUM_CHORIDE_ACTIVATION
     if (k_i.start_col == 0 && 
	 k_i.end_col == 0 &&
	 k_i.t0 == 0_s &&
	 k_i.v_rise == 0.0 &&
	 k_i.v_decay == 0.0 &&
	 k_i.vK == 0_mV &&
	 k_i.vCl == 0_mV) {
     	  std::cerr << "**********\n"
     		    << "WARNING: using potassium chloride activation and all "
     		    << "related parameters are null!\n"
     		    << "**********\n";
     }
     model_koenigsberger::activator::K_Cl_Activator act(k_i.start_col,
							k_i.end_col,
							k_i.t0,
							k_i.v_rise,
							k_i.v_decay,
							model_constants,
							k_i.vK,
							k_i.vCl);
#else
#  error Unknown or unspecified activation type for Koenigsberger! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION

     // ================================

     const auto split = utility::split_filename(output_file);

     output::Output_Streams output(output_file);
     output << "# t idx= ";
     for(auto i: k_i.cell_indices) {
	  output << "(" << std::get<0>(i) << ", " << std::get<1>(i) << ") ";
     }
     output << std::endl;
     output.precision(12);

     // ==========

     output::Flux_Streams flux_out(k_i.cell_indices.size(), output_file);
     flux_out << Cell::get_fluxes_header();
     flux_out.precision(12);

     // ================================

     mysystem.print_cell_at(k_i.cell_indices, output);
     mysystem.print_fluxes_at(k_i.cell_indices, flux_out);

     Progress_Bar p(nstep_equ, "Equilibration:   ");
     p.start();
     for(size_type i(0) ; i < nstep_equ ; ++i, ++p) {
     	  mysystem.do_step(dt_equ);
	  if (i % s_i.output_step_equ == 0) {
	       mysystem.print_cell_at(k_i.cell_indices, output);
	       mysystem.print_fluxes_at(k_i.cell_indices, flux_out);
	  }
     }
     mysystem.print_cell_at(k_i.cell_indices, output);
     mysystem.print_fluxes_at(k_i.cell_indices, flux_out);

     // ==========

     p.reset(nstep_act, "Activation:      ");
     p.start();
     for(size_type i(0) ; i < nstep_act ; ++i, ++p) {
     	  mysystem.activate(dt_act, act);
     	  if (i % s_i.output_step_act == 0) {
	       mysystem.print_cell_at(k_i.cell_indices, output);
     	       mysystem.print_fluxes_at(k_i.cell_indices, flux_out);
     	  }
     }
     mysystem.print_cell_at(k_i.cell_indices, output);
     mysystem.print_fluxes_at(k_i.cell_indices, flux_out);

     // ==========

     p.reset(nstep_run, "Simulation:      ");
     p.start();
     for(size_type i(0) ; i < nstep_run ; ++i, ++p) {
     	  mysystem.do_step(dt_run);
     	  if (i % s_i.output_step_run == 0) {
	       mysystem.print_cell_at(k_i.cell_indices, output);
     	       mysystem.print_fluxes_at(k_i.cell_indices, flux_out);
     	  }
     }

     // END of simulation
     // ========================================================================

     return 0;
}

#pragma GCC diagnostic pop

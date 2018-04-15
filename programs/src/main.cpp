// Input includes
#include "argument_parsing.hpp"
#include "read_parameters.hpp"
#include "input_options.hpp"
#include "halidi_input_options.hpp"
#include "biffdiag_input_options.hpp"
#include "koenigs_input_options.hpp"

// Execution mode includes
#include "halidi_bifurcation_diagram.hpp"
#include "halidi.hpp"
#include "koenigs_bifurcation_diagram.hpp"
#include "koenigsberger.hpp"

#include "utility.hpp"
#include "read_constants_halidi.hpp"
#include "read_constants_koenigsberger.hpp"
#include "constants_halidi.hpp"
#include "constants_koenigsberger.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include <boost/program_options.hpp>

int main(int argc, char** argv)
{
     using namespace cellsim;
     namespace po = boost::program_options;
     
     const auto start_program(std::chrono::steady_clock::now());
     auto end_init(start_program);

     std::string input_filename, output_filename;
     std::tie(input_filename, output_filename) = input::parse_arguments(argc, 
									argv);

     int return_value{0};

     if (!input_filename.empty()) {
	  po::variables_map vm_input;
	  EXEC_MODE exec_mode;

	  {
	       po::options_description input_file_options("Input file options");
	       input_file_options.add_options()
		    ("system.exec_mode", po::value<std::string>(),
		     "Program execution mode.");
	       input_file_options.add(input::generic_input_options());
	       input_file_options.add(input::simulation_input_options());
	       input_file_options.add(input::halidi_input_options());
	       input_file_options.add(input::koenigsberger_input_options());
	       input_file_options.add(input::biffdiag_input_options());

	       if (input_filename == "1") {
		    return 1;
	       }
	       else if (input_filename == "help" || 
			input_filename == "dumpconfig") {
		    return 0;
	       }
	       else if (input_filename == "help-more") {
		    std::cout << input_file_options << std::endl;
		    return 0;
	       }
	       std::ifstream input_file{input_filename};
	       if (!input_file) {
	       	    std::cerr << "Error while opening input file!\n";
	       	    return 1;
	       }
	  	  
	       try {
	       	    po::store(po::parse_config_file(input_file, 
	       					    input_file_options),
	       		      vm_input);
	       	    po::notify(vm_input);

	       	    exec_mode = input::read_exec_mode(
	       		 vm_input["system.exec_mode"].as<std::string>());
	       }
	       catch(std::exception& e)
	       {
	       	    std::cerr << "Config file:\n\t" << e.what() << std::endl;
	       	    return 1;
	       }
	  }

	  // Read constants file
	  using cellsim::input::get_constants_filename;
	  using cellsim::input::read_constants;

	  std::cout << "Program execution mode: " << utility::to_string(exec_mode)
		    << std::endl;

	  end_init = std::chrono::steady_clock::now();

	  // ===================================================================
	  // Execute program

	  if (cellsim::BIFFDIAG_HALIDI <= exec_mode && 
	      exec_mode <= cellsim::HALIDI) {
	       using namespace cellsim::model_halidi;

	       Constants_Values model_constants;
	       if (!read_constants(
	       		get_constants_filename(vm_input), model_constants)) {
		    return 1;
	       }

	       if (exec_mode == cellsim::BIFFDIAG_HALIDI) {
		    return_value = cellsim::run_halidi_biffdiag(vm_input, 
								output_filename,
								model_constants);
	       }
	       else { // exec_mode == cellsim::HALIDI
		    return_value = cellsim::run_halidi(vm_input, 
						       output_filename,
						       model_constants);
	       }

	  }
	  else if (cellsim::BIFFDIAG_KOENIGS <= exec_mode &&
		   exec_mode <= cellsim::KOENIGSBERGER) {

	       using namespace cellsim::model_koenigsberger;
	       Constants_Values model_constants;
	       if (!read_constants(
			get_constants_filename(vm_input), model_constants)) {
		    return 1;
	       }


	       if (exec_mode == cellsim::BIFFDIAG_KOENIGS) {
		    return_value = cellsim::run_koenigs_biffdiag(vm_input, 
								 output_filename,
								 model_constants);
	       }
	       else { // exec_mode = KOENIGSBERGER
		    return_value = cellsim::run_koenigsberger(vm_input, 
							      output_filename,
							      model_constants);
	       }
	  }
	  else {
	       std::cerr << "Unknown execution mode! (got: '" 
			 << vm_input["system.exec_mode"].as<std::string>()
			 << "')\n";
	       return 1;
	  }
     }

     const auto end_program(std::chrono::steady_clock::now());

     using namespace std::chrono;
     using std::setw;
     using std::setfill;

     std::cout << "\n****************************************"
	       << "****************************************\n";

     auto elapsed_time(end_init - start_program);
     hours   hh = duration_cast<hours>(elapsed_time);
     minutes mm = duration_cast<minutes>(elapsed_time % hours(1));
     seconds ss = duration_cast<seconds>(elapsed_time % minutes(1));
     milliseconds msec = duration_cast<milliseconds>(elapsed_time % seconds(1));

     std::cout << "*                    "
	       << "Initialisation time:     " << setfill('0')
	       << setw(2) << hh.count() << ":"
	       << setw(2) << mm.count() << ":"
	       << setw(2) << ss.count() << "."
	       << setw(3) << msec.count()
	       << "                     *"
	       << std::endl;

     elapsed_time = end_program - end_init;
     hh = duration_cast<hours>(elapsed_time);
     mm = duration_cast<minutes>(elapsed_time % hours(1));
     ss = duration_cast<seconds>(elapsed_time % minutes(1));
     msec = duration_cast<milliseconds>(elapsed_time % seconds(1));

     std::cout << "*                    "
	       << "Execution time:          " << setfill('0')
	       << setw(2) << hh.count() << ":"
	       << setw(2) << mm.count() << ":"
	       << setw(2) << ss.count() << "."
	       << setw(3) << msec.count()
	       << "                     *"
	       << std::endl;

     elapsed_time = end_program - start_program;
     hh = duration_cast<hours>(elapsed_time);
     mm = duration_cast<minutes>(elapsed_time % hours(1));
     ss = duration_cast<seconds>(elapsed_time % minutes(1));
     msec = duration_cast<milliseconds>(elapsed_time % seconds(1));

     std::cout << "*                    "
	       << "Total execution time:    " << setfill('0')
	       << setw(2) << hh.count() << ":"
	       << setw(2) << mm.count() << ":"
	       << setw(2) << ss.count() << "."
	       << setw(3) << msec.count()
	       << "                     *"
	       << std::endl;

     std::cout << "****************************************"
	       << "****************************************\n";

     return return_value;
}

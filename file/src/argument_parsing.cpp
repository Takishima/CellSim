#include "argument_parsing.hpp"

#ifdef MULTI_THREAD
#  include "system_koenigsberger.hpp"
#endif // MULTI_THREAD

#include <iostream>
#include <boost/program_options.hpp>
using namespace boost::program_options;

/* Auxiliary functions for checking input for validity.
 * Function used to check that 'opt1' and 'opt2' are not specified
 * at the same time. 
 */
void conflicting_options(const boost::program_options::variables_map& vm,
                         const char* opt1, const char* opt2)
{
     if (vm.count(opt1) && !vm[opt1].defaulted() 
	 && vm.count(opt2) && !vm[opt2].defaulted())
	  throw std::logic_error(std::string("Conflicting options '") 
				 + opt1 + "' and '" + opt2 + "'.");
}

std::tuple<std::string,std::string> 
cellsim::input::parse_arguments(int argc, char** argv)
{
     std::string input_filename;
     std::string output_filename;

     options_description cmdline_options(
	  "This version of Cellsim has been compiled on "
	  __DATE__ ":\n"
	  "Command-line options");

     cmdline_options.add_options()
	  ("dumpconfig", "Print values of defines used at compilation")
	  ("help,h", "Print help message and exit.")
	  ("help-more", "Print full help message and exit.")
	  ("input-file,i", 
	   value(&input_filename)->default_value("input"),
	   "Input file.")
	  ("output-file,o",value(&output_filename)->default_value("output.dat"),
	   "Base name for output files.")
#ifdef MULTI_THREAD
	  ("n-threads,t", 
	   value<size_type>(&n_threads())->default_value(4),
	   "Number of threads to launch.")
#endif //MULTI_THREAD
	  ;

     variables_map vm;
     try {
	  store(parse_command_line(argc, argv, cmdline_options), vm);
	  notify(vm);
	  conflicting_options(vm, "help", "dumpconfig");
	  conflicting_options(vm, "help", "more-help");
     }
     catch(std::exception& e)
     {
	  std::cerr << "Program error:" << std::endl << e.what() << std::endl;
	  return std::make_tuple("1", "");
     }


     if (vm.count("help")) {  
	  std::cout << cmdline_options
		    << std::endl;
	  return std::make_tuple("help", "");
     }
     else if (vm.count("help-more")) {  
	  std::cout << cmdline_options
		    << std::endl;
	  return std::make_tuple("help-more", "");
     }
     else if (vm.count("dumpconfig")) {
	  std::cout << "Program compiled on " __DATE__ " at " __TIME__ "\n"
#ifdef __clang__
	       "\tClang version: " 
	       + std::to_string(__clang_major__) + "."
	       + std::to_string(__clang_minor__) + "\n"
#else
	       "\tGCC version: "
	       + std::to_string(__GNUC__) + "."
	       + std::to_string(__GNUC_MINOR__) + "\n"
#endif //__clang__
	       ;

	  std::cout << "\tCompiled in "
#ifndef NDEBUG
		    << "*DEBUG* mode.\n";
#else
		    << "*RELEASE* mode.\n";
#endif // NDEBUG

	  std::cout << "\tDefines used to configure Cellsim:\n";

#ifdef HALIDI_2009
	  std::cout << "\t\tHALIDI_2009\n";
#elif defined HALIDI_2012
	  std::cout << "\t\tHALIDI_2012\n";
#endif // HALIDI_20XX

#ifdef MULTI_THREAD
	  std::cout << "\t\tMULTI_THREAD\n";
#endif // MULTI_THREAD

#ifdef PLC_ACTIVATION
	  std::cout << "\t\tPLC_ACTIVATION\n";
#elif defined ELECTRIC_ACTIVATION
	  std::cout << "\t\tELECTRIC_ACTIVATION\n";
#elif defined POTASSIUM_CHORIDE_ACTIVATION
	  std::cout << "\t\tPOTASSIUM_CHORIDE_ACTIVATION\n";
#endif // PLC_ACTIVATION

#ifdef WITH_RANGE_CHECK
	  std::cout << "\t\tWITH_RANGE_CHECK\n";
#endif // WITH_RANGE_CHECK

	  return std::make_tuple("dumpconfig", "");
     }

#ifdef MULTI_THREAD
     if (n_threads() == 0) {
	  throw std::logic_error("Cannot have '0' as -t [--n-threads] value!");
     }
#endif // MULTI_THREAD

     return std::make_tuple(input_filename, output_filename);
}

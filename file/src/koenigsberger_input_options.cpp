#include "koenigs_input_options.hpp"
#include "units.hpp"

#include <boost/program_options.hpp>
#include <iostream>

namespace input = cellsim::input;

using namespace boost::program_options;

// =============================================================================
// Local helper function for parsing list of tuple of numbers

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_alternative.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

namespace {
     using cellsim::input::Point;

     namespace qi = boost::spirit::qi;
     namespace ascii = boost::spirit::ascii;
     namespace phoenix = boost::phoenix;

     GCC_DIAG_OFF(effc++)

     template <typename Iterator>
     struct point_grammar : qi::grammar<Iterator, std::vector<Point>(), ascii::space_type>
     {
	  point_grammar() : point_grammar::base_type(start)
	       {
		    using qi::uint_;
		    using qi::repeat;

		    point_rule %= '(' >> uint_ >> repeat(0,1)[','] >> uint_ >> ')';
		    start %= point_rule >> *(repeat(0,1)[','] >> point_rule);
		    ;
	       }

	  qi::rule<Iterator, Point(), ascii::space_type> point_rule;
	  qi::rule<Iterator, std::vector<Point>(), ascii::space_type> start;
     };

     GCC_DIAG_ON(effc++)

     template <typename Iterator>
     bool parse_numbers(Iterator first,
			Iterator last, 
			std::vector<Point>& v)
     {
	  using ascii::space;

	  point_grammar<Iterator> grammar;
	  bool r = phrase_parse(
	       first, 
	       last,
	       grammar,
	       space,
	       v);

	  // Fail if we did not get a full match
	  if (first != last) { 
	       v = std::vector<Point>{};
	       return false;
	  }

	  return r;
     }


     std::vector<Point> read_vector(std::string str, 
						 const char* option_name)
     {
	  namespace po = boost::program_options;
     
	  std::vector<Point> v;

	  if (!parse_numbers(begin(str), end(str), v)) {
	       throw po::validation_error(
		    po::validation_error::invalid_option_value,
		    option_name);
	  }

	  return v;
     }
}

// =============================================================================

options_description cellsim::input::koenigsberger_input_options()
{
     options_description koenigs_options;
     // General options
     koenigs_options.add_options()
	  ("koenigs.Nx", value<size_type>(),
	   "Number of cells in the x direction")
	  ("koenigs.Ny", value<size_type>(),
	   "Number of cells in the y direction")
	  ("koenigs.JPLCago", value<double>(),
	   "Background PLC agonist current (if 0 use model value of 0.05).")
	  ;

     // Activation options
     koenigs_options.add_options()
	  ("koenigs.activation.start_column", value<size_type>(),
	   "Index of first column to activate.")
	  ("koenigs.activation.end_column", value<size_type>(),
	   "Index past the last column to activate "
	   "(ie. will activate from start to end-1)")
#ifdef PLC_ACTIVATION
	  ("koenigs.activation.Jact", value<double>(), 
	   "Activation flux to add to stimulated cells"
	   " (if 0, use default value of 0.4")
#elif defined ELECTRIC_ACTIVATION
	  ("koenigs.activation.vact", value<double>(), 
	   "Activation electric potential to add to stimulated cells")
	  ("koenigs.activation.v_rise", value<double>(), 
	   "Rise constant for electrical activation")
	  ("koenigs.activation.v_decay", value<double>(), 
	   "Decay constant for electrical activation")
	  ("koenigs.activation.t0", value<double>(), 
	   "Time after which electric perturbation starts decaying")
#elif defined POTASSIUM_CHORIDE_ACTIVATION
	  ("koenigs.activation.t0", value<double>(), 
	   "Time after which K-Cl perturbation starts decaying after reaching max")
	  ("koenigs.activation.v_rise", value<double>(), 
	   "Rise constant for K-Cl activation")
	  ("koenigs.activation.v_decay", value<double>(), 
	   "Decay constant for K-Cl activation")
	  ("koenigs.activation.vK", value<double>(), 
	   "Activation value for vK")
	  ("koenigs.activation.vCl", value<double>(), 
	   "Activation value for vCl")
#else
#  error Unknown or unspecified activation type for Koenigsberger! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION
	  ;

     // Output options
     koenigs_options.add_options()
	  ("koenigs.output.cell_indices", value<std::string>()->multitoken(),
	   "Indices of the cells to use for flux output\n"
		"(list of (#,#) & cannot be empty)")
	  ;

     return koenigs_options;
}

cellsim::input::Koenigs_Input 
cellsim::input::read_koenigs_values(po::variables_map vm_input, bool& error)
{     
     if (vm_input["koenigs.Nx"].empty()) {
	  std::cerr << "Missing value : koenigs.Nx\n";
	  error = true;
     }
     if (vm_input["koenigs.Ny"].empty()) {
	  std::cerr << "Missing value : koenigs.Ny\n";
	  error = true;
     }
     if (vm_input["koenigs.JPLCago"].empty()) {
	  std::cerr << "Missing value : koenigs.JPLCago\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.start_column"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.start_column\n";
	  error = true;
     }
     if (vm_input["koenigs.activation.end_column"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.end_column\n";
	  error = true;
     }

#ifdef PLC_ACTIVATION
     if (!vm_input["koenigs.activation.vact"].empty() ||
	 !vm_input["koenigs.activation.v_rise"].empty() ||
	 !vm_input["koenigs.activation.v_decay"].empty() ||
	 !vm_input["koenigs.activation.t0"].empty()) {
	  std::cerr << "Program compiled to use PLC activation. Got one of:\n"
		    << "  koenigs.activation.vact\n"
		    << "  koenigs.activation.v_rise\n"
		    << "  koenigs.activation.v_decay\n"
		    << "  koenigs.activation.t0\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.Jact"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.Jact\n";
	  error = true;
     }
#elif defined ELECTRIC_ACTIVATION
     if (!vm_input["koenigs.activation.Jact"].empty() ||
	 !vm_input["koenigs.activation.center_index"].empty() ||
	 !vm_input["koenigs.activation.region_size"].empty()) {
	  std::cerr << "Program compiled to use PLC activation. Got one of:\n"
		    << "  koenigs.activation.Jact\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.vact"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.vact\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.v_rise"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.v_rise\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.v_decay"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.v_decay\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.t0"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.t0\n";
	  error = true;
     }
#elif defined POTASSIUM_CHORIDE_ACTIVATION
     if (vm_input["koenigs.activation.t0"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.t0\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.v_rise"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.v_rise\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.v_decay"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.v_decay\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.vK"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.vK\n";
	  error = true;
     }

     if (vm_input["koenigs.activation.vCl"].empty()) {
	  std::cerr << "Missing value : koenigs.activation.vCl\n";
	  error = true;
     }
#else
#  error Unknown or unspecified activation type for Koenigsberger! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION
       
     // ================================

     if (!error) {
	  try {
	       return Koenigs_Input{
		    vm_input["koenigs.Nx"].as<size_type>(),
			 vm_input["koenigs.Ny"].as<size_type>(),
			 flux_ut::from_value(vm_input["koenigs.JPLCago"].as<double>()),
			 vm_input["koenigs.activation.start_column"].as<size_type>(),
			 vm_input["koenigs.activation.end_column"].as<size_type>(),
#ifdef PLC_ACTIVATION
			 flux_ut::from_value(vm_input["koenigs.activation.Jact"]
					     .as<double>()),
#elif defined ELECTRIC_ACTIVATION
			 electric_potential_ut::from_value(vm_input["koenigs.activation.vact"]
							   .as<double>()),
			 time_ut::from_value(vm_input["koenigs.activation.t0"].as<double>()),
			 vm_input["koenigs.activation.v_rise"].as<double>(),
			 vm_input["koenigs.activation.v_decay"].as<double>(),
#elif defined POTASSIUM_CHORIDE_ACTIVATION
			 time_ut::from_value(vm_input["koenigs.activation.t0"].as<double>()),
			 vm_input["koenigs.activation.v_rise"].as<double>(),
			 vm_input["koenigs.activation.v_decay"].as<double>(),
			 electric_potential_ut::from_value(vm_input["koenigs.activation.vK"]
							   .as<double>()),
			 electric_potential_ut::from_value(vm_input["koenigs.activation.vCl"]
							   .as<double>()),
#else
#  error Unknown or unspecified activation type for Koenigsberger! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION
			 read_vector(vm_input["koenigs.output.cell_indices"]
				     .as<std::string>(),
				     "koenigs.output.cell_indices")
			 };
	  }
	  catch(po::validation_error& e) {
	       std::cerr << "Error: " << e.what() << "\n"
			 << "    Got: '" 
			 << vm_input["koenigs.output.cell_indices"].as<std::string>()
			 << "'\n";
	       
	       error = true;
	       return Koenigs_Input();
	  }
	  catch(boost::bad_any_cast& e) {
	       std::cerr << "Casting error: " << e.what() << std::endl;
	       error = true;
	       return Koenigs_Input();
	  }
	  catch(std::exception& e) {
	       std::cerr << "UNEXPECTED ERROR: " << e.what() << std::endl;
	       error = true;
	       return Koenigs_Input();
	  }

     }
     else {
	  return Koenigs_Input();
     }
}

void cellsim::input::print_option(Koenigs_Input value)
{
     std::cout << "Koenigsberger:\n"
	       << "    Nx =           " << value.Nx << std::endl
	       << "    Ny =           " << value.Ny << std::endl
	       << "    JPLCago =      " << value.JPLCago << std::endl
	       << "  activation:\n"
#ifdef PLC_ACTIVATION
	       << "      Jact =         " << value.Jact << std::endl
	       << "      start_column = " << value.start_col << std::endl
	       << "      end_column =   " << value.end_col << std::endl;
#elif defined ELECTRIC_ACTIVATION
	       << "      vact =         " << value.vact << std::endl
	       << "      t0 =           " << value.t0 << std::endl
	       << "      v_rise =       " << value.v_rise << std::endl
	       << "      v_decay =      " << value.v_decay << std::endl
#elif defined POTASSIUM_CHORIDE_ACTIVATION
	       << "      vK =           " << value.vK << std::endl
	       << "      vCl =          " << value.vCl << std::endl
	       << "      t0 =           " << value.t0 << std::endl
	       << "      v_rise =       " << value.v_rise << std::endl
	       << "      v_decay =      " << value.v_decay << std::endl
#else
#  error Unknown or unspecified activation type for Koenigsberger! (look at *_ACTIVATION defines)
#endif // PLC_ACTIVATION
	       << "  output:\n"
	       << "      cell_indices = ";
     for (const auto el : value.cell_indices) {
	  std::cout << '(' << std::get<0>(el) << ", " << std::get<1>(el) << ") ";
     }
     std::cout << std::endl;
}


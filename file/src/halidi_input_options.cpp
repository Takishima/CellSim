#include "halidi_input_options.hpp"

#include <iostream>
#include <boost/program_options.hpp>

namespace input = cellsim::input;

using namespace boost::program_options;


// =============================================================================
// Local helper function for parsing list of numbers

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_alternative.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

namespace {
     namespace qi = boost::spirit::qi;
     namespace ascii = boost::spirit::ascii;
     namespace phoenix = boost::phoenix;

     template <typename Iterator>
     bool parse_numbers(Iterator first,
			Iterator last, 
			std::vector<cellsim::size_type>& v)
     {
	  using qi::uint_;
	  using qi::double_;
	  using qi::phrase_parse;
	  using qi::_1;
	  using ascii::space;
	  using phoenix::push_back;
	  using qi::repeat;

	  bool r = phrase_parse(
	       first, 
	       last,

	       //  Begin grammar
	       (
		    *uint_[push_back(phoenix::ref(v), _1)]
		    >> *(repeat(0,1)[','] >> uint_[push_back(phoenix::ref(v), _1)])
		    )
	       ,
	       //  End grammar

	       space);

	  // Fail if we did not get a full match
	  if (first != last) { 
	       v = std::vector<cellsim::size_type>{};
	       return false;
	  }

	  return r;
     }


     std::vector<cellsim::size_type> read_vector(std::string str, 
						 const char* option_name)
     {
	  namespace po = boost::program_options;
     
	  std::vector<cellsim::size_type> v;

	  if (!parse_numbers(begin(str), end(str), v)) {
	       throw po::validation_error(
		    po::validation_error::invalid_option_value,
		    option_name);
	  }

	  return v;
     }
}

// =============================================================================


options_description cellsim::input::halidi_input_options()
{
     options_description halidi_options;
     // General options
     halidi_options.add_options()
	  ("halidi.npoints", value<size_type>(),
	   "Number of points per cells.")
	  ("halidi.ncells", value<size_type>(),
	   "Number of cells.")
	  ("halidi.dx", value<double>(), 
	   "Spatial step (will be deprecated in a future revision)")
	  ("halidi.JPLCago", value<double>(),
	   "Background PLC agonist current.")
	  ;

     // Activation options
     halidi_options.add_options()
	  ("halidi.activation.cell_index", value<size_type>(),
	   "Index of cell to activate.")
#ifdef PLC_ACTIVATION
	  ("halidi.activation.Jact", value<double>(), 
	   "Activation flux to add to stimulated cells")
	  ("halidi.activation.center_index", value<size_type>(),
	   "Index of the center point of the activation region")
	  ("halidi.activation.region_size", value<size_type>(),
	   "Size of the activation region (# of points).")
#elif defined ELECTRIC_ACTIVATION
	  ("halidi.activation.vact", value<double>(), 
	   "Activation electric potential to add to stimulated cells")
	  ("halidi.activation.v_rise", value<double>(), 
	   "Grow constant for electrical activation (tanh only)")
	  ("halidi.activation.v_decay", value<double>(), 
	   "Decay constant for electrical activation")
	  ("halidi.activation.t0", value<double>(), 
	   "Time after which electric perturbation starts decaying")
#elif defined POTASSIUM_CHORIDE_ACTIVATION
	  ("halidi.activation.t0", value<double>(), 
	   "Time after which K-Cl perturbation starts decaying after reaching max")
	  ("halidi.activation.v_rise", value<double>(), 
	   "Grow constant for K-Cl activation")
	  ("halidi.activation.v_decay", value<double>(), 
	   "Decay constant for K-Cl activation")
	  ("halidi.activation.vK", value<double>(), 
	   "Activation value for vK")
	  ("halidi.activation.vCl", value<double>(), 
	   "Activation value for vCl")
#endif // PLC_ACTIVATION
	  ;

     // Output options
     halidi_options.add_options()
	  ("halidi.output.point_index", value<size_type>()->default_value(0), 
	   "Which point inside each cell to output")
	  ("halidi.output.cell_indices", value<std::string>()->multitoken(),
	   "Indices of the cells to use for flux output (can be empty)")
	  ;

     return halidi_options;
}

cellsim::input::Halidi_Input 
cellsim::input::read_halidi_values(po::variables_map vm_input, bool& error)
{
     if (vm_input["halidi.npoints"].empty()) {
	  std::cerr << "Missing value : halidi.npoints\n";
	  error = true;
     }
     if (vm_input["halidi.ncells"].empty()) {
	  std::cerr << "Missing value : halidi.ncells\n";
	  error = true;
     }
     if (vm_input["halidi.dx"].empty()) {
	  std::cerr << "Missing value : halidi.dx\n";
	  error = true;
     }
     if (vm_input["halidi.JPLCago"].empty()) {
	  std::cerr << "Missing value : halidi.JPLCago\n";
	  error = true;
     }

     // ================================
     // Handle activation options

     if (vm_input["halidi.activation.cell_index"].empty()) {
	  std::cerr << "Missing value : halidi.activation.cell_index\n";
	  error = true;
     }

#ifdef PLC_ACTIVATION
     if (!vm_input["halidi.activation.vact"].empty() ||
	 !vm_input["halidi.activation.v_rise"].empty() ||
	 !vm_input["halidi.activation.v_decay"].empty() ||
	 !vm_input["halidi.activation.t0"].empty()) {
	  std::cerr << "Program compiled to use PLC activation. Got one of:"
		    << "  halidi.activation.vact"
		    << "  halidi.activation.v_decay"
		    << "  halidi.activation.t0\n";
	  error = true;
     }

     if (vm_input["halidi.activation.Jact"].empty()) {
	  std::cerr << "Missing value : halidi.activation.Jact\n";
	  error = true;
     }

     if (vm_input["halidi.activation.center_index"].empty()) {
	  std::cerr << "Missing value : halidi.activation.center_index\n";
	  error = true;
     }

     if (vm_input["halidi.activation.region_size"].empty()) {
	  std::cerr << "Missing value : halidi.activation.region_size\n";
	  error = true;
     }
#elif defined ELECTRIC_ACTIVATION
     if (!vm_input["halidi.activation.Jact"].empty() ||
	 !vm_input["halidi.activation.center_index"].empty() ||
	 !vm_input["halidi.activation.region_size"].empty()) {
	  std::cerr << "Program compiled to use PLC activation. Got one of:"
		    << "  halidi.activation.Jact"
		    << "  halidi.activation.center_index"
		    << "  halidi.activation.region_size\n";
	  error = true;
     }

     if (vm_input["halidi.activation.vact"].empty()) {
	  std::cerr << "Missing value : halidi.activation.vact\n";
	  error = true;
     }

     if (vm_input["halidi.activation.v_rise"].empty()) {
	  std::cerr << "Missing value : halidi.activation.v_rise\n";
	  error = true;
     }

     if (vm_input["halidi.activation.v_decay"].empty()) {
	  std::cerr << "Missing value : halidi.activation.v_decay\n";
	  error = true;
     }

     if (vm_input["halidi.activation.t0"].empty()) {
	  std::cerr << "Missing value : halidi.activation.t0\n";
	  error = true;
     }
#elif defined POTASSIUM_CHORIDE_ACTIVATION
     if (vm_input["halidi.activation.t0"].empty()) {
	  std::cerr << "Missing value : halidi.activation.t0\n";
	  error = true;
     }

     if (vm_input["halidi.activation.v_decay"].empty()) {
	  std::cerr << "Missing value : halidi.activation.v_decay\n";
	  error = true;
     }

     if (vm_input["halidi.activation.v_rise"].empty()) {
	  std::cerr << "Missing value : halidi.activation.v_rise\n";
	  error = true;
     }

     if (vm_input["halidi.activation.vK"].empty()) {
	  std::cerr << "Missing value : halidi.activation.vK\n";
	  error = true;
     }

     if (vm_input["halidi.activation.vCl"].empty()) {
	  std::cerr << "Missing value : halidi.activation.vCl\n";
	  error = true;
     }
#endif // PLC_ACTIVATION
       
     // ================================

     if (!error) {
	  try {
	       Halidi_Input h_i = {
		    vm_input["halidi.npoints"].as<size_type>(),
		    vm_input["halidi.ncells"].as<size_type>(),
		    length_ut::from_value(vm_input["halidi.dx"].as<double>()),
		    flux_ut::from_value(vm_input["halidi.JPLCago"].as<double>()),
		    vm_input["halidi.activation.cell_index"].as<size_type>(),
#ifdef PLC_ACTIVATION
		    flux_ut::from_value(vm_input["halidi.activation.Jact"]
					.as<double>()),
		    vm_input["halidi.activation.center_index"].as<size_type>(),
		    vm_input["halidi.activation.region_size"].as<size_type>(),
#elif defined ELECTRIC_ACTIVATION
		    electric_potential_ut::from_value(vm_input["halidi.activation.vact"]
						      .as<double>()),
		    time_ut::from_value(vm_input["halidi.activation.t0"].as<double>()),
		    vm_input["halidi.activation.v_rise"].as<double>(),
		    vm_input["halidi.activation.v_decay"].as<double>(),
#elif defined POTASSIUM_CHORIDE_ACTIVATION
		    time_ut::from_value(vm_input["halidi.activation.t0"].as<double>()),
		    vm_input["halidi.activation.v_rise"].as<double>(),
		    vm_input["halidi.activation.v_decay"].as<double>(),
		    electric_potential_ut::from_value(vm_input["halidi.activation.vK"]
						      .as<double>()),
		    electric_potential_ut::from_value(vm_input["halidi.activation.vCl"]
						      .as<double>()),
#endif // PLC_ACTIVATION
		    vm_input["halidi.output.point_index"].as<size_type>(),
		    read_vector(vm_input["halidi.output.cell_indices"]
				.as<std::string>(),
				"halidi.output.cell_indices")
	       };

	       // Check for out-of-range for the output
	       if (h_i.point_index > h_i.npoints) {
		    error = true;
		    return Halidi_Input();
	       }
	       return h_i;
	  }
	  catch(po::validation_error& e) {
	       std::cerr << "Error: " << e.what() << "\n"
			 << "    Got: '" 
			 << vm_input["halidi.output.cell_indices"].as<std::string>()
			 << "'\n";
	       
	       error = true;
	       return Halidi_Input();
	  }
	  catch(std::exception& e) {
	       std::cerr << "UNEXPECTED ERROR: " << e.what() << std::endl;
	       error = true;
	       return Halidi_Input();
	  }
     }
     else {
	  return Halidi_Input();
     }

}

void cellsim::input::print_option(Halidi_Input value)
{
     std::cout << "Halidi:\n"
	       << "    npoints =      " << value.npoints << std::endl
	       << "    ncells =       " << value.ncells << std::endl
	       << "    dx =           " << value.dx << std::endl
	       << "    JPLCago =      " << value.JPLCago << std::endl
	       << "  activation:\n"
	       << "      cell_index =   " << value.cell_index << std::endl
#ifdef PLC_ACTIVATION
	       << "      Jact =         " << value.Jact << std::endl
	       << "      center_index = " << value.center_index << std::endl
	       << "      region_size =  " << value.region_size << std::endl
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
#endif // PLC_ACTIVATION
	       << "  output:\n"
	       << "      point_index =  " << value.point_index << std::endl
	       << "      cell_indices = ";
     for (const auto el : value.cell_indices) {
	  std::cout << el << " ";
     }
     std::cout << std::endl;
}


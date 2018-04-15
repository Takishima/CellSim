#include "read_constants_koenigsberger.hpp"

#include "constants_koenigsberger.hpp"
#include "units_io.hpp"

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace koenigs = cellsim::model_koenigsberger;

constexpr cellsim::input::constant_value list[] = {
     {"F", "[uM / s]"},
     {"KI", "[uM]"},
     {"GCa", "[uM / mV / s]"},
     {"vCa1", "[mV]"},
     {"vCa2", "[mV]"},
     {"RCa", "[mV]"},
     {"GNaCa", "[uM]"},
     {"cNaCa", "[uM]"},
     {"vNaCa", "[uM]"},
     {"B", "[uM / s]"},
     {"cb", "[uM]"},
     {"C", "[uM / s]"},
     {"sc", "[uM]"},
     {"cc", "[uM]"},
     {"D", "[1 / s]"},
     {"vd", "[mV]"},
     {"Rd", "[mV]"},
     {"L", "[1 / s]"},
     {"gamma", "[mV / uM]"},
     {"FNaK", "[uM / s]"},
     {"GCl", "[uM / mV / s]"},
     {"cCl", "[uM]"},
     {"vCl", "[mV]"},
     {"GK", "[uM / mV / s]"},
     {"vK", "[mV]"},
     {"lambda", "[1 / s]"},
     {"cw", "[uM]"},
     {"beta", "[uM * uM]"},
     {"vCa3", "[mV]"},
     {"RK", "[mV]"},
     {"Gback", "[uM / mV / s]"},
     {"vrest", "[mV]"},
     {"k", "[1 / s]"},
     {"E", "[uM / s]"},
     {"KCa", "[uM]"},
     {"g", "[1 / s]"},
};

namespace {
     GCC_DIAG_OFF(conversion)

     template <typename quantity_type>
     quantity_type helper(po::variables_map vm, const char* name, bool& error)
     {
	  // Just over-checking here to detect typos...
	  if (vm.count(name) == 0) {
	       std::cerr << "ERROR: Unkown constant name: " << name << std::endl;
	       error = true;
	       return quantity_type::from_value(0.);
	  }
	  return quantity_type::from_value(vm[name].as<double>());
     }

     GCC_DIAG_ON(conversion)
}


void cellsim::input::print_constants_values(const koenigs::Constants_Values& C)
{
     std::cout << "========================================"
	       << "========================================\n"
	       << "MODEL CONSTANTS:\n"
	       <<  "\tF = " << C.F       << std::endl
	       <<  "\tKI = " << C.KI     << std::endl
	       <<  "\tGCa = " << C.GCa   << std::endl 
	       <<  "\tvCa1 = " << C.vCa1 << std::endl 
	       <<  "\tvCa2 = " << C.vCa2 << std::endl 
	       <<  "\tRCa = " << C.RCa	 << std::endl 
	       <<  "\tGNaCa = " << C.GNaCa << std::endl 
	       <<  "\tcNaCa = " << C.cNaCa << std::endl 
	       <<  "\tvNaCa = " << C.vNaCa << std::endl 
	       <<  "\tB = " << C.B	 << std::endl 
	       <<  "\tcb = " << C.cb	 << std::endl 
	       <<  "\tC = " << C.C	 << std::endl 
	       <<  "\tsc = " << C.sc	 << std::endl 
	       <<  "\tcc = " << C.cc	 << std::endl 
	       <<  "\tD = " << C.D	 << std::endl 
	       <<  "\tvd = " << C.vd	 << std::endl 
	       <<  "\tRd = " << C.Rd	 << std::endl 
	       <<  "\tL = " << C.L	 << std::endl 
	       <<  "\tgamma = " << C.gamma << std::endl 
	       <<  "\tFNaK = " << C.FNaK << std::endl 
	       <<  "\tGCl = " << C.GCl	 << std::endl 
	       <<  "\tcCl = " << C.cCl	 << std::endl 
	       <<  "\tvCl = " << C.vCl	 << std::endl 
	       <<  "\tGK = " << C.GK	 << std::endl 
	       <<  "\tvK = " << C.vK	 << std::endl 
	       <<  "\tlambda = " << C.lambda<< std::endl 
	       <<  "\tcw = " << C.cw	 << std::endl 
	       <<  "\tbeta = " << C.beta << std::endl 
	       <<  "\tvCa3 = " << C.vCa3 << std::endl 
	       <<  "\tGback = " << C.Gback << std::endl 
	       <<  "\tvrest = " << C.vrest << std::endl 
	       <<  "\tRK = " << C.RK	 << std::endl 
	       <<  "\tk = " << C.k	 << std::endl 
	       <<  "\tE = " << C.E	 << std::endl 
	       <<  "\tKCa = " << C.KCa	 << std::endl 
	       <<  "\tg = " << C.g	 << std::endl 
	       << "########################################"
	       << "########################################\n";
}


bool cellsim::input::set_constants_values(po::variables_map vm, 
					  koenigs::Constants_Values& m_c)
{
     // Check that all values are present
     for (const auto el : list) {
	  if (vm[el.name].empty()) {
	       std::cerr << "Missing value : " << el.name << std::endl;
	       return false;
	  }
     }

     bool error(false);
     koenigs::Constants_Values tmp;

     tmp.F = helper<flux_ut>(vm, "F", error);
     tmp.KI = helper<concentration_ut>(vm, "KI", error);

     tmp.GCa = helper<conductance_ut>(vm, "GCa", error);
     tmp.vCa1 = helper<electric_potential_ut>(vm, "vCa1", error);
     tmp.vCa2 = helper<electric_potential_ut>(vm, "vCa2", error);
     tmp.vCa3 = helper<electric_potential_ut>(vm, "vCa3", error);
     tmp.RCa = helper<electric_potential_ut>(vm, "RCa", error);

     tmp.GNaCa = helper<conductance_ut>(vm, "GNaCa", error);
     tmp.cNaCa = helper<concentration_ut>(vm, "cNaCa", error);
     tmp.vNaCa = helper<electric_potential_ut>(vm, "vNaCa", error);

     tmp.B = helper<flux_ut>(vm, "B", error);
     tmp.cb = helper<concentration_ut>(vm, "cb", error);

     tmp.C = helper<flux_ut>(vm, "C", error);

     tmp.sc = helper<concentration_ut>(vm, "sc", error);
     tmp.cc = helper<concentration_ut>(vm, "cc", error);

     tmp.D = helper<inverse_time_ut>(vm, "D", error);

     tmp.vd = helper<electric_potential_ut>(vm, "vd", error);
     tmp.Rd = helper<electric_potential_ut>(vm, "Rd", error);

     tmp.L = helper<inverse_time_ut>(vm, "L", error);
     tmp.gamma = helper<gamma_ut>(vm, "gamma", error);
     tmp.FNaK = helper<flux_ut>(vm, "FNaK", error);

     tmp.GCl = helper<conductance_ut>(vm, "GCl", error);
     tmp.cCl = helper<concentration_ut>(vm, "cCl", error);
     tmp.vCl = helper<electric_potential_ut>(vm, "vCl", error);

     tmp.GK = helper<conductance_ut>(vm, "GK", error);
     tmp.vK = helper<electric_potential_ut>(vm, "vK", error);

     tmp.lambda = helper<inverse_time_ut>(vm, "lambda", error);
     tmp.cw = helper<concentration_ut>(vm, "cw", error);
     tmp.beta = helper<beta_ut>(vm, "beta", error);

     tmp.RK = helper<electric_potential_ut>(vm, "RK", error);

     tmp.Gback = helper<conductance_ut>(vm, "Gback", error);
     tmp.vrest = helper<electric_potential_ut>(vm, "vrest", error);

     tmp.k = helper<inverse_time_ut>(vm, "k", error);
     tmp.E = helper<flux_ut>(vm, "E", error);
     tmp.KCa = helper<concentration_ut>(vm, "KCa", error);
     tmp.g = helper<inverse_time_ut>(vm, "g", error);

     if (error) {
	  return false;
     }
     m_c = std::move(tmp);
     return true;
}

bool cellsim::input::read_constants(std::string filename, 
				    koenigs::Constants_Values& model_constants)
{
     const auto options = constants_input_options(list, get_array_size(list));
     po::variables_map vm_constants;
     
     std::ifstream input_file{filename};
     if (!input_file) {
	  std::cerr << "Error while opening constants input file!\n";
	  return false;
     }
	  	  
     try {
	  po::store(po::parse_config_file(input_file, options), vm_constants);
	  po::notify(vm_constants);
	  return set_constants_values(vm_constants, model_constants);
     }
     catch(std::exception& e)
     {
	  std::cerr << "Constants file:\n\t" << e.what() << std::endl;
	  return false;
     }

}


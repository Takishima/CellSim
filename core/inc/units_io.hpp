/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef UNITS_IO_HPP_INCLUDED
#define UNITS_IO_HPP_INCLUDED

#include "units.hpp"

#ifndef NO_UNITS

#include <boost/units/io.hpp>
#include <boost/units/reduce_unit.hpp>
#include <boost/units/systems/si/io.hpp>

namespace boost { 
     namespace units {
	  namespace c_ = cellsim::units_impl;

	  template <typename T>
	  using r_u_ = typename reduce_unit<T>::type;

	  inline std::string name_string(const r_u_<c_::molar_unit>&)
	  {return "molar";}
	  inline std::string symbol_string(const r_u_<c_::molar_unit>&)
	  {return "M";}

	  inline std::string name_string(const r_u_<c_::micro_molar_unit>&) 
	  {return "micro molar";}
	  inline std::string symbol_string(const r_u_<c_::micro_molar_unit>&)
	  {return "uM";}

	  inline std::string name_string(const r_u_<c_::molar_per_second_unit>&) 
	  {return "molar per seconds";}
	  inline std::string symbol_string(const r_u_<c_::molar_per_second_unit>&)
	  {return "M/s";}

	  inline std::string name_string(const r_u_<c_::micro_molar_per_second_unit>&) 
	  {return "micro molar per seconds";}
	  inline std::string symbol_string(const r_u_<c_::micro_molar_per_second_unit>&)
	  {return "uM/s";}

	  inline std::string name_string(const r_u_<c_::potential_SR_unit>&)
	  {return "milli volts per seconds";}
	  inline std::string symbol_string(const r_u_<c_::potential_SR_unit>&) 
	  {return "mV/s";}

	  inline std::string name_string(const r_u_<c_::diffusion_unit>&)
	  {return "micro meter squared per seconds";}
	  inline std::string symbol_string(const r_u_<c_::diffusion_unit>&) 
	  {return "um^2/s";}

	  inline std::string name_string(const r_u_<c_::permeablity_unit>&)
	  {return "micro meter per seconds";}
	  inline std::string symbol_string(const r_u_<c_::permeablity_unit>&) 
	  {return "um/s";}

	  inline std::string symbol_string(const r_u_<c_::conductance_ut::unit_type>&)
	  {return "uM/(mV s)";}
	  inline std::string symbol_string(const r_u_<c_::gamma_ut::unit_type>&)
	  {return "mV/uM";}
	  inline std::string symbol_string(const r_u_<c_::beta_ut::unit_type>&)
	  {return "uM^2";}
     }    
}
#endif // NO_UNITS
#endif //UNITS_IO_HPP_INCLUDED

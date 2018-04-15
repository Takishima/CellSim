#include "read_parameters.hpp"

#include <boost/algorithm/string.hpp>

namespace input = cellsim::input;

cellsim::EXEC_MODE input::read_exec_mode(std::string str)
{
     boost::algorithm::to_lower(str);
     if (str == "biffdiag_halidi") {
	  return BIFFDIAG_HALIDI;
     }
     else if (str == "biffdiag_koenigs") {
	  return BIFFDIAG_KOENIGS;
     }
     else if (str == "halidi") {
	  return HALIDI;
     }
     else if (str == "koenigsberger") {
	  return KOENIGSBERGER;
     }
     else {
	  return INVALID;
     }     
}



/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef DEFINITIONS_HPP_INCLUDED
#define DEFINITIONS_HPP_INCLUDED

#include "config.hpp"

#include <cstddef>

namespace cellsim {
     typedef std::size_t size_type;

     enum TYPE {
	  CA,
	  IP3,
	  S,
	  V,
	  W
     };
     
     enum TIME {
	  NOW,
	  NEXT
     };

     enum EXEC_MODE {
	  INVALID = -1,
	  BIFFDIAG_HALIDI = 0, 
	  HALIDI, 
	  BIFFDIAG_KOENIGS = 10, 
	  KOENIGSBERGER
     };

     size_type& n_threads();
} // namespace cellsim


#endif //DEFINITIONS_HPP_INCLUDED

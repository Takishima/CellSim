/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef HALIDI_BIFURCATION_DIAGRAM_HPP_INCLUDED
#define HALIDI_BIFURCATION_DIAGRAM_HPP_INCLUDED

#include <string>

namespace boost { namespace program_options {
	  class variables_map;
} }

namespace cellsim {
     namespace model_halidi {
	  struct Constants_Values;
     }
     
     namespace halidi = model_halidi;

     int run_halidi_biffdiag(boost::program_options::variables_map vm_input,
			     std::string output_file,
			     const halidi::Constants_Values& model_constants);
} // namespace cellsim

#endif //HALIDI_BIFURCATION_DIAGRAM_HPP_INCLUDED

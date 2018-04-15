/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef KOENIGS_BIFURCATION_DIAGRAM_HPP_INCLUDED
#define KOENIGS_BIFURCATION_DIAGRAM_HPP_INCLUDED

#include <string>

namespace boost { namespace program_options {
	  class variables_map;
} }

namespace cellsim {
     namespace model_koenigsberger {
	  struct Constants_Values;
     }

     namespace koenigs = model_koenigsberger;

     int run_koenigs_biffdiag(boost::program_options::variables_map vm_input,
			      std::string output_file,
			      const koenigs::Constants_Values& model_constants);
} // namespace cellsim

#endif //KOENIGS_BIFURCATION_DIAGRAM_HPP_INCLUDED

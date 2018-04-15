/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef KOENIGSBERGER_HPP_INCLUDED
#define KOENIGSBERGER_HPP_INCLUDED

#include <string>

namespace boost { namespace program_options {
	  class variables_map;
} }

namespace cellsim {
     namespace model_koenigsberger {
	  struct Constants_Values;
     }

     namespace koenigs = model_koenigsberger;
     int run_koenigsberger(boost::program_options::variables_map vm_input,
			   std::string output_file,
			   const koenigs::Constants_Values& model_constants);
} // namespace cellsim

#endif //KOENIGSBERGER_HPP_INCLUDED

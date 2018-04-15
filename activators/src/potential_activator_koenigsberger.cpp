#include "potential_activator_koenigsberger.hpp"

namespace activator = cellsim::model_koenigsberger::activator;

activator::Potential_Activator::Potential_Activator(size_type start_column,
						    size_type end_column,
						    electric_potential_ut v_act,
						    time_ut t0,
						    double alpha)
     : Activator_Base(start_column, end_column),
       Potential_Activator_Base(v_act, t0, alpha)
{}

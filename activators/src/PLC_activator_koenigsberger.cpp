#include "PLC_activator_koenigsberger.hpp"

namespace activator = cellsim::model_koenigsberger::activator;

activator::PLC_Activator::PLC_Activator(size_type start_column,
					size_type end_column,
					flux_ut J_act)
     : Activator_Base(start_column, end_column), J_act_(J_act)
{}



#include "potential_activator_koenigsberger_tanh.hpp"

#include "cell_koenigsberger.hpp"

#include <cmath>

namespace activator = cellsim::model_koenigsberger::activator;

activator::Potential_Activator_Tanh::Potential_Activator_Tanh(
     size_type start_column,
     size_type end_column,
     electric_potential_ut v_act,
     time_ut t0,
     double alpha_rise,
     double alpha)
     : Activator_Base(start_column, end_column),
       Potential_Activator_Tanh_Base(v_act, t0, alpha_rise, alpha)
{}

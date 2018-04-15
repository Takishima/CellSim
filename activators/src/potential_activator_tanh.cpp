#include "potential_activator_tanh.hpp"
#include "cell.hpp"

#include <cmath>

namespace activator = cellsim::model_halidi::activator;

activator::Potential_Activator_Tanh::Potential_Activator_Tanh(
     size_type cell_idx,
     electric_potential_ut v_act,
     time_ut t0,
     double alpha_rise,
     double alpha)
     : Activator_Base(cell_idx),
       Potential_Activator_Tanh_Base(v_act, t0, alpha_rise, alpha)
{}

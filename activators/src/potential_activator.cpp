#include "potential_activator.hpp"

namespace activator = cellsim::model_halidi::activator;

activator::Potential_Activator::Potential_Activator(size_type cell_idx,
						    electric_potential_ut v_act,
						    time_ut t0,
						    double alpha)
     : Activator_Base(cell_idx),
       Potential_Activator_Base(v_act, t0, alpha)
{}
